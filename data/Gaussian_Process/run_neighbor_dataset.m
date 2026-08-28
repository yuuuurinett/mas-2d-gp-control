function run_neighbor_dataset(DatasetName, CurrentMode, train_ratio, seed)

if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1;          end
rng(seed);

fprintf('\n[邻域静态聚合 NBR] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

%% 1. 加载数据集
switch upper(DatasetName)
    case 'KIN40K'
        tr=load('KIN40K_train.mat'); te=load('KIN40K_test.mat');
        hp=load('KIN40K_Hyperparameter.mat');
        train_x=tr.x; train_y=tr.y; test_x=te.xtest; test_y=te.ytest;
    case 'POL'
        tr=load(fullfile('POL','POL_train.mat')); te=load(fullfile('POL','POL_test.mat'));
        hp=load(fullfile('POL','POL_Hyperparameter.mat'));
        train_x=tr.x; train_y=tr.y; test_x=te.xtest; test_y=te.ytest;
    case 'PUMADYN32NM'
        tr=load(fullfile('PUMADYN32NM','PUMADYN32NM_train.mat')); 
        te=load(fullfile('PUMADYN32NM','PUMADYN32NM_test.mat'));
        hp=load(fullfile('PUMADYN32NM','PUMADYN32NM_Hyperparameter.mat'));
        train_x=tr.x; train_y=tr.y; test_x=te.xtest; test_y=te.ytest;
    case 'SARCOS'
        tr=load(fullfile('SARCOS','SARCOS_train.mat'));
        te=load(fullfile('SARCOS','SARCOS_test.mat'));
        hp_raw=load(fullfile('SARCOS','SARCOS_GP_Hyperparameter.mat'));
        hp.SigmaF=mean(cell2mat(hp_raw.SigmaF_set));
        hp.SigmaN=mean(cell2mat(hp_raw.SigmaN_set));
        hp.SigmaL=mean(cell2mat(hp_raw.SigmaL_set'),2);
        train_x=tr.sarcos_inv(:,1:21); train_y=tr.sarcos_inv(:,22:28);
        test_x=te.sarcos_inv_test(:,1:21); test_y=te.sarcos_inv_test(:,22:28);
    otherwise, error('未知数据集: %s', DatasetName);
end

%% 2. 归一化 & 数据准备
if size(hp.SigmaL,1)>1 && size(hp.SigmaL,2)>1, hp.SigmaL=mean(hp.SigmaL,1); end
SigmaL=hp.SigmaL(:);
if numel(hp.SigmaF)>1, hp.SigmaF=mean(hp.SigmaF); end
if numel(hp.SigmaN)>1, hp.SigmaN=mean(hp.SigmaN); end
SigmaF=hp.SigmaF; SigmaN=hp.SigmaN;

N_all=size(train_x,1);
n_train=round(N_all*train_ratio);
idx_tr=randperm(N_all,n_train);
X_train=train_x(idx_tr,:); Y_train=train_y(idx_tr,:);

idx_te=randperm(size(test_x,1));
X_test=test_x(idx_te,:); Y_test=test_y(idx_te,:);

X_mean=mean(X_train,1); X_std=std(X_train,0,1); X_std(X_std==0)=1;
if ~(max(abs(X_mean))<1e-2 && max(abs(X_std-1))<1e-2)
    X_train=(X_train-X_mean)./X_std;
    X_test=(X_test-X_mean)./X_std;
    SigmaL=SigmaL./X_std(:);
end

Y_mean=mean(Y_train,1); Y_std=std(Y_train,0,1); Y_std(Y_std==0)=1;
if max(abs(Y_mean))<1e-2 && max(abs(Y_std-1))<1e-2
    Y_mean=zeros(1,size(Y_train,2)); Y_std=ones(1,size(Y_train,2));
else
    Y_train=(Y_train-Y_mean)./Y_std;
    SigmaF=SigmaF/mean(Y_std); SigmaN=SigmaN/mean(Y_std);
end
prior_var=SigmaF^2;

[N_train,x_dim]=size(X_train); y_dim=size(Y_train,2);
N_eval=min(3000,size(X_test,1));
X_eval=X_test(1:N_eval,:); Y_eval_raw=Y_test(1:N_eval,:);
Y_var_base=var(Y_eval_raw,0,1);

%% 3. 分布式参数 & 拓扑定义
AgentQuantity=6;
MultiAgentSystem=Manipulator_2D_2DoF_SetMASTopology(AgentQuantity,1);
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
A = double(L ~= 0); 
for i=1:AgentQuantity, A(i,i)=1; end

%% 4. 训练局部 GP 
t_start = tic;
LocalGP_set = cell(AgentQuantity, 1);

W = floor(N_train / AgentQuantity); 

fprintf('\n>>> 数据分配: N_train = %d, 无重叠划分 <<<\n', N_train);

for n = 1:AgentQuantity
    start_idx = (n - 1) * W + 1;
    
    if n == AgentQuantity
        
        end_idx = N_train; 
    else
        end_idx = start_idx + W - 1;
    end
    
    idx = start_idx : end_idx;
    fprintf('  -> Agent %d 分配索引: %6d 到 %6d  (容量: %d 点)\n', n, start_idx, end_idx, length(idx));

    % 加载模型
    LocalGP_set{n}=LocalGP_MultiOutput(x_dim,y_dim,length(idx),SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)',Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8; LocalGP_set{n}.delta=0.01;
end
t_train = toc(t_start);
t_nbr_local_gp_train = t_train;
fprintf('局部GP训练完成(已确保100%%数据覆盖): %.4fs\n', t_train);

%% 5. 预计算所有 Agent 的局部预测 (调用外部批量加速函数版)
fprintf('\n预计算 %d 个测试点的局部推断...\n', N_eval);
tic;
Mu_Local  = zeros(AgentQuantity, y_dim, N_eval);
Var_Local = zeros(AgentQuantity, y_dim, N_eval);
X_eval_matrix = X_eval'; 

chunk_size = 500; 
for n = 1:AgentQuantity
    for start_idx = 1:chunk_size:N_eval
        end_idx = min(start_idx + chunk_size - 1, N_eval);
        X_chunk = X_eval_matrix(:, start_idx:end_idx);
        
        % 【关键修复】调用外部批量加速函数
        [mn_batch, vn_batch] = batch_predict_external(LocalGP_set{n}, X_chunk, SigmaN, SigmaF);
        
        for t_idx = 1:size(X_chunk, 2)
            t_global = start_idx + t_idx - 1;
            Mu_Local(n, :, t_global)  = mn_batch(t_idx, :);
            Var_Local(n, :, t_global) = repmat(vn_batch(t_idx), 1, y_dim);
        end
    end
    fprintf('  -> Agent %d 预计算完成\n', n);
end
Precompute_Time = toc;
fprintf('预计算全部完成: %.2fs\n', Precompute_Time);
%% 6. 方法字典定义
methods_all={'moe','gpoe','poe','bcm','rbcm'};
if strcmpi(CurrentMode,'all'), method_list=methods_all;
else, method_list={lower(CurrentMode)}; end

SaveFolder=fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
tr_tag=round(train_ratio*100);

%% 7. 邻居静态聚合主循环
comm_train = 0;
comm_test  = 1; % 邻居通信仅需 1 轮

for mi=1:numel(method_list)
    cur=method_list{mi};
    fprintf('\n[%d/%d] NBR-%s\n',mi,numel(method_list),upper(cur));
    t_nbr_test_local_prediction = Precompute_Time / max(numel(method_list), 1);
    tic_aggregation = tic;
    
    agent_smse = zeros(AgentQuantity, 1);
    agent_rmse = zeros(AgentQuantity, 1);
    agent_nlpd = zeros(AgentQuantity, 1);
    
    for n = 1:AgentQuantity
        neighbors_idx = find(A(n, :) == 1);
        num_neighbors = length(neighbors_idx);
        
        Mean_Aggregated = zeros(N_eval, y_dim);
        Var_Aggregated  = zeros(N_eval, y_dim);
        
        for t = 1:N_eval
            for d = 1:y_dim
                mu_n_set  = squeeze(Mu_Local(neighbors_idx, d, t)); 
                var_n_set = squeeze(Var_Local(neighbors_idx, d, t));
                beta = max(0.5*(log(prior_var)-log(var_n_set)), eps);
                
                switch cur
                    case 'moe'
                        Mean_Aggregated(t,d) = mean(mu_n_set);
                        Var_Aggregated(t,d)  = mean(var_n_set + mu_n_set.^2) - Mean_Aggregated(t,d)^2;
                    case 'gpoe'
                        prec = sum(beta ./ var_n_set);
                        Mean_Aggregated(t,d) = sum(beta .* mu_n_set ./ var_n_set) / max(prec, eps);
                        Var_Aggregated(t,d)  = 1 / max(prec, eps);
                    case 'poe'
                        prec = sum(1 ./ var_n_set);
                        Mean_Aggregated(t,d) = sum(mu_n_set ./ var_n_set) / max(prec, eps);
                        Var_Aggregated(t,d)  = 1 / max(prec, eps);
                    case 'bcm'
                        prec = sum(1 ./ var_n_set) - (num_neighbors-1)/prior_var;
                        Mean_Aggregated(t,d) = sum(mu_n_set ./ var_n_set) / max(prec, eps);
                        Var_Aggregated(t,d)  = 1 / max(prec, eps);
                    case 'rbcm'
                        beta_sum = sum(beta);
                        if beta_sum > 1, beta = beta / beta_sum; beta_sum = 1; end
                        prec = sum(beta ./ var_n_set) + (1 - beta_sum)/prior_var;
                        Mean_Aggregated(t,d) = sum(beta .* mu_n_set ./ var_n_set) / max(prec, eps);
                        Var_Aggregated(t,d)  = 1 / max(prec, eps);
                end
                Var_Aggregated(t,d) = max(Var_Aggregated(t,d), SigmaN^2);
            end
        end
        
        Mean_Aggregated_orig = Mean_Aggregated .* repmat(Y_std, N_eval, 1) + repmat(Y_mean, N_eval, 1);
        Var_orig = Var_Aggregated .* repmat(Y_std.^2, N_eval, 1);
        err = Y_eval_raw - Mean_Aggregated_orig;
        
        agent_smse(n) = mean(mean(err.^2) ./ Y_var_base);
        agent_rmse(n) = mean(sqrt(mean(err.^2)));
        agent_nlpd(n) = mean(mean(0.5*(log(2*pi*Var_orig) + err.^2./Var_orig)));
   
        if n == 1, all_agents_err_sq = zeros(N_eval, AgentQuantity); end
        all_agents_err_sq(:, n) = mean(err.^2, 2);
    end
    
    t_nbr_aggregation = toc(tic_aggregation);
    t_test = t_nbr_test_local_prediction + t_nbr_aggregation;
    
    t_train_per_point = (t_train / N_train) * 1000;
    t_test_per_point  = (t_test / N_eval) * 1000;
    t_train_total = t_train;
    t_test_total = t_test;

    smse = mean(agent_smse);
    rmse = mean(agent_rmse);
    nlpd = mean(agent_nlpd);
    
    err_sq_mean = mean(all_agents_err_sq, 2); 
    smse_curve  = cumsum(err_sq_mean) ./ (1:N_eval)' / mean(Y_var_base);
    rmse_curve  = sqrt(cumsum(err_sq_mean) ./ (1:N_eval)');

    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train: %.2f ms/pt  Test: %.2f ms/pt\n', ...
        smse, rmse, nlpd, t_train_per_point, t_test_per_point);
    fprintf('  Timing NBR: LocalPred=%.3fs  Aggregation=%.3fs  TotalTest=%.3fs\n', ...
        t_nbr_test_local_prediction, t_nbr_aggregation, t_test_total);

    
    save(fullfile(SaveFolder, sprintf('%s_nbr_tr%d_mc%d.mat', cur, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', 't_train', 't_test', ...
        't_train_total', 't_test_total', ...
        't_train_per_point', 't_test_per_point', ...
        't_nbr_local_gp_train', 't_nbr_test_local_prediction', ...
        't_nbr_aggregation', ...
        'comm_train', 'comm_test', ...
        'cur', 'seed', 'train_ratio', 'smse_curve', 'rmse_curve');
end
fprintf('\n[%s] NBR done. seed=%d tr=%d%%\n\n', DatasetName, seed, tr_tag);
end
