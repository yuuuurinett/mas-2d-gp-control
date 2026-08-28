function run_centralized_dataset(DatasetName, CurrentMode, train_ratio, seed)

% DatasetName  : 'KIN40K' | 'POL' | 'PUMADYN32NM' | 'SARCOS'
% CurrentMode  : 'moe'|'gpoe'|'poe'|'bcm'|'rbcm'|'all'
% train_ratio  : 训练集比例，如 0.4
% seed         : 随机种子

if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1;          end
rng(seed);

fprintf('\n[集中式基线 CEN] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

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

%% 2. 超参数提取 & 归一化
if size(hp.SigmaL,1)>1 && size(hp.SigmaL,2)>1, hp.SigmaL=mean(hp.SigmaL,1); end
SigmaL=hp.SigmaL(:);
if numel(hp.SigmaF)>1, hp.SigmaF=mean(hp.SigmaF); end
if numel(hp.SigmaN)>1, hp.SigmaN=mean(hp.SigmaN); end
SigmaF=hp.SigmaF; SigmaN=hp.SigmaN;

n_train=round(size(train_x,1)*train_ratio);
idx_tr=randperm(size(train_x,1),n_train);
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
X_eval=X_test(1:N_eval,:); Y_eval=Y_test(1:N_eval,:);
Y_var_base=var(Y_eval,0,1);
fprintf('Train=%d  Test=%d  x=%d  y=%d\n',N_train,N_eval,x_dim,y_dim);

%% 3. 分布式参数
AgentQuantity=6;
MaxDataPerAgent=min(floor(N_train/AgentQuantity),3000);

%% 4. 训练局部 GP
t_start=tic;
LocalGP_set=cell(AgentQuantity,1);
for n=1:AgentQuantity
    idx=(n-1)*MaxDataPerAgent+1 : min(n*MaxDataPerAgent,N_train);
    LocalGP_set{n}=LocalGP_MultiOutput(x_dim,y_dim,MaxDataPerAgent,SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)',Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8; LocalGP_set{n}.delta=0.01;
end
t_train=toc(t_start);
t_cen_local_gp_train = t_train;
fprintf('局部GP训练完成（无通信）: %.4fs\n',t_train);

%% 5. 方法列表
methods_all={'moe','gpoe','poe','bcm','rbcm'};
if strcmpi(CurrentMode,'all'), method_list=methods_all;
else, method_list={lower(CurrentMode)}; end

SaveFolder=fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
tr_tag=round(train_ratio*100);

%% 6. 统一预计算 (调用外部批量加速函数版)
fprintf('\n预计算 %d 个测试点的局部推断...\n', N_eval);
tic;
Mu_All  = zeros(AgentQuantity, y_dim, N_eval);
Var_All = zeros(AgentQuantity, y_dim, N_eval);
X_eval_matrix = X_eval'; % 转换为列格式 [x_dim x N_eval]

chunk_size = 500; % 分块防止爆内存
for n = 1:AgentQuantity
    for start_idx = 1:chunk_size:N_eval
        end_idx = min(start_idx + chunk_size - 1, N_eval);
        X_chunk = X_eval_matrix(:, start_idx:end_idx);
        
        % 【关键修复】使用你写好的外部加速函数，完美绕过底层的 Bug
        [mn_batch, vn_batch] = batch_predict_external(LocalGP_set{n}, X_chunk, SigmaN, SigmaF);
        
        % 将批量结果塞入总矩阵
        for t_idx = 1:size(X_chunk, 2)
            t_global = start_idx + t_idx - 1;
            Mu_All(n, :, t_global)  = mn_batch(t_idx, :);
            % var_batch 是一维的，扩展到 y_dim 维度
            Var_All(n, :, t_global) = repmat(vn_batch(t_idx), 1, y_dim);
        end
    end
    fprintf('  -> Agent %d 预计算完成\n', n);
end
Precompute_Time = toc;
fprintf('预计算完成: %.2fs\n', Precompute_Time);

%% 7. 集中式聚合主循环
comm_train = 0;
comm_test  = 1;

for mi=1:numel(method_list)
    cur=method_list{mi};
    fprintf('\n[%d/%d] CEN-%s\n',mi,numel(method_list),upper(cur));
    t_cen_test_local_prediction = Precompute_Time / max(numel(method_list), 1);
    tic_aggregation = tic;

    Final_Mean=zeros(N_eval,y_dim);
    Final_Var =zeros(N_eval,y_dim);

    for t=1:N_eval
        for d=1:y_dim
            mu_a  = squeeze(Mu_All(:,d,t));    % [AgentQuantity×1]
            var_a = squeeze(Var_All(:,d,t));   % [AgentQuantity×1]
            
            beta  = max(0.5*(log(prior_var)-log(var_a)), eps);

            switch cur
                case 'moe'
                    Final_Mean(t,d) = mean(mu_a);
                   
                    Final_Var(t,d)  = mean(var_a);

                case 'poe'
                    prec = sum(1./var_a);
                    Final_Mean(t,d) = sum(mu_a./var_a) / max(prec, eps);
                    Final_Var(t,d)  = 1 / max(prec, eps);

                case 'gpoe'
                    prec = sum(beta./var_a);
                    Final_Mean(t,d) = sum(beta.*mu_a./var_a) / max(prec, eps);
                    Final_Var(t,d)  = 1 / max(prec, eps);

                case 'bcm'
                    prec = sum(1./var_a) - (AgentQuantity-1)/prior_var;
                    Final_Mean(t,d) = sum(mu_a./var_a) / max(prec, eps);
                    Final_Var(t,d)  = 1 / max(prec, eps);

                case 'rbcm'
                    
                    prec = sum(beta./var_a) + (1 - sum(beta))/prior_var;
                    Final_Mean(t,d) = sum(beta.*mu_a./var_a) / max(prec, eps);
                    Final_Var(t,d)  = 1 / max(prec, eps);
            end
            Final_Var(t,d) = max(Final_Var(t,d), SigmaN^2);
        end
    end

    t_cen_aggregation = toc(tic_aggregation);
    t_test = t_cen_test_local_prediction + t_cen_aggregation;

    %% 8. 反归一化 & 误差计算
    mu_pred  = Final_Mean .* repmat(Y_std, N_eval, 1) + repmat(Y_mean, N_eval, 1);
    var_pred = Final_Var  .* repmat(Y_std.^2, N_eval, 1);

    err  = Y_eval - mu_pred;
    smse = mean(mean(err.^2) ./ Y_var_base);
    rmse = mean(sqrt(mean(err.^2)));
    nlpd = mean(mean(0.5*(log(2*pi*var_pred) + err.^2 ./ var_pred)));

    t_train_per_point = (t_train / N_train) * 1000;
    t_test_per_point  = (t_test / N_eval) * 1000;
    t_train_total = t_train;
    t_test_total = t_test;

    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train: %.2f ms/pt  Test: %.2f ms/pt\n', ...
        smse, rmse, nlpd, t_train_per_point, t_test_per_point);
    fprintf('  Timing CEN: LocalPred=%.3fs  Aggregation=%.3fs  TotalTest=%.3fs\n', ...
        t_cen_test_local_prediction, t_cen_aggregation, t_test_total);

    err_sq_mean = mean(err.^2, 2);
    smse_curve  = cumsum(err_sq_mean) ./ (1:N_eval)' / mean(Y_var_base);
    rmse_curve  = sqrt(cumsum(err_sq_mean) ./ (1:N_eval)');

    save(fullfile(SaveFolder, sprintf('%s_cen_tr%d_mc%d.mat', cur, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', 't_train', 't_test', ...
        't_train_total', 't_test_total', ...
        't_train_per_point', 't_test_per_point', ...
        't_cen_local_gp_train', 't_cen_test_local_prediction', ...
        't_cen_aggregation', ...
        'comm_train', 'comm_test', ...
        'cur', 'seed', 'train_ratio', 'smse_curve', 'rmse_curve');
end
fprintf('\n[%s] CEN done. seed=%d tr=%d%%\n\n', DatasetName, seed, tr_tag);
end
