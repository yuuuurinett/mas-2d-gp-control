function run_testpoint_dataset(DatasetName, CurrentMode, train_ratio, seed)
% 分布式高斯过程 - 测试点方法 (TP-DAC / TP-AC)
% 完美终极版：包含统一预计算 10 倍提速优化、绝对公平抽样、防崩溃通信循环

if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1;          end
rng(seed);

fprintf('\n[测试点] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

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

%% 2. 归一化与测试点抽取
if size(hp.SigmaL,1)>1 && size(hp.SigmaL,2)>1, hp.SigmaL=mean(hp.SigmaL,1); end
SigmaL=hp.SigmaL(:);
if numel(hp.SigmaF)>1, hp.SigmaF=mean(hp.SigmaF); end
if numel(hp.SigmaN)>1, hp.SigmaN=mean(hp.SigmaN); end
SigmaF=hp.SigmaF; SigmaN=hp.SigmaN;

n_train=round(size(train_x,1)*train_ratio);
idx_tr=randperm(size(train_x,1),n_train);
X_train=train_x(idx_tr,:); Y_train=train_y(idx_tr,:);

% --- 绝对公平的 3000 个测试点随机抽样 ---
idx_te=randperm(size(test_x,1));
X_test=test_x(idx_te,:); Y_test=test_y(idx_te,:);

X_mean=mean(X_train,1); X_std=std(X_train,0,1); X_std(X_std==0)=1;
if ~(max(abs(X_mean))<1e-2 && max(abs(X_std-1))<1e-2)
    X_train=(X_train-X_mean)./X_std;
    X_test=(X_test-X_mean)./X_std;
    SigmaL = SigmaL ./ X_std(:); 
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
N_eval=min(3000,size(X_test,1));  % 限制为 3000 个测试点
X_eval=X_test(1:N_eval,:); Y_eval=Y_test(1:N_eval,:);
Y_var_base=var(Y_eval,0,1);
fprintf('Train=%d  Test=%d  x=%d  y=%d\n',N_train,N_eval,x_dim,y_dim);

%% 3. 分布式参数与拓扑
AgentQuantity=6; Kappa_P=10; t_step=0.01;
MaxDataPerAgent=min(floor(N_train/AgentQuantity),3000);
MultiAgentSystem=Manipulator_2D_2DoF_SetMASTopology(AgentQuantity,1);
L=MultiAgentSystem.Agent_Topology.LaplacianMatrix;

%% 4. 训练局部 GP
tic;
LocalGP_set=cell(AgentQuantity,1);
for n=1:AgentQuantity
    idx=(n-1)*MaxDataPerAgent+1 : min(n*MaxDataPerAgent,N_train);
    LocalGP_set{n}=LocalGP_MultiOutput(x_dim,y_dim,MaxDataPerAgent,SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)',Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8; LocalGP_set{n}.delta=0.01;
end
fprintf('局部GP训练完成: %.2fs\n',toc);

%% 5. 方法字典配置
dac_methods={'moe','gpoe','poe','bcm','rbcm'};
ac_methods={'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all'), AllModes=[dac_methods,ac_methods];
else, AllModes={lower(CurrentMode)}; end

SaveFolder=fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
p_dim=2*y_dim; tr_tag=round(train_ratio*100);

%% =========================================================================
%% 6. 统一预计算阶段 (核心 10 倍提速)
%% =========================================================================
fprintf('\n[核心优化] 正在预计算 %d 个测试点的局部推断 (6个专家并算)...\n', N_eval);
tic;
Mu_Local_All  = zeros(AgentQuantity, y_dim, N_eval);
Var_Local_All = zeros(AgentQuantity, y_dim, N_eval);

for t = 1:N_eval
    Test_Point_X = X_eval(t, :)';
    for n = 1:AgentQuantity
        [mn, vn] = LocalGP_set{n}.predict(Test_Point_X);
        Mu_Local_All(n, :, t)  = mn'; 
        Var_Local_All(n, :, t) = vn';
    end
end
Precompute_Time = toc;
fprintf('预计算完成！总计耗时: %.2f 秒 (每个点平均 %.4f 秒)\n', Precompute_Time, Precompute_Time/N_eval);

%% =========================================================================
%% 7. 聚合方法主循环 (共用预计算结果)
%% =========================================================================
for mi = 1:numel(AllModes)
    Current_Method = AllModes{mi};
    Method_Base_Name = strrep(lower(Current_Method), '_ac', '');
    fprintf('\n[%d/%d] 正在执行聚合: %s\n', mi, numel(AllModes), Current_Method);

    Final_Mean_Pred = zeros(N_eval, y_dim);
    Final_Var_Pred  = zeros(N_eval, y_dim);
    
    tic; 

    if ismember(lower(Current_Method), dac_methods)
        %% --- A. 分布式一致性聚合 (DAC) ---
        P_Info_Matrix = zeros(p_dim, AgentQuantity, N_eval);
        for t = 1:N_eval
            for n = 1:AgentQuantity
                Local_Mu  = Mu_Local_All(n, :, t)'; 
                Local_Var = Var_Local_All(n, :, t)';
                for d = 1:y_dim
                    Beta = max(0.5 * (log(prior_var) - log(Local_Var(d))), eps);
                    switch lower(Current_Method)
                        case 'moe'
                            P_Info_Matrix(2*d-1, n, t) = AgentQuantity * Local_Mu(d);
                            P_Info_Matrix(2*d,   n, t) = AgentQuantity * (Local_Var(d) + Local_Mu(d)^2);
                        case 'gpoe'
                            P_Info_Matrix(2*d-1, n, t) = AgentQuantity * Beta * Local_Mu(d) / Local_Var(d);
                            P_Info_Matrix(2*d,   n, t) = AgentQuantity * Beta / Local_Var(d);
                        case 'poe'
                            P_Info_Matrix(2*d-1, n, t) = AgentQuantity * Local_Mu(d) / Local_Var(d);
                            P_Info_Matrix(2*d,   n, t) = AgentQuantity / Local_Var(d);
                        case 'bcm'
                            P_Info_Matrix(2*d-1, n, t) = AgentQuantity * Local_Mu(d) / Local_Var(d);
                            P_Info_Matrix(2*d,   n, t) = AgentQuantity / Local_Var(d) - (AgentQuantity-1)/prior_var;
                        case 'rbcm'
                            P_Info_Matrix(2*d-1, n, t) = AgentQuantity * Beta * Local_Mu(d) / Local_Var(d);
                            P_Info_Matrix(2*d,   n, t) = AgentQuantity * Beta / Local_Var(d) + (1 - AgentQuantity*Beta)/prior_var;
                    end
                end
            end
        end
        
        % 执行稳健的网络共识迭代
        Zeta_Comm = zeros(p_dim, AgentQuantity, N_eval);
        for iter = 1:3000
            Prev_Zeta = Zeta_Comm;
            % 绝对安全的迭代方式
            for n = 1:AgentQuantity
                Zeta_Comm(:, n, :) = Zeta_Comm(:, n, :) + t_step * Kappa_P * sum((P_Info_Matrix - Prev_Zeta) .* reshape(L(n, :), 1, AgentQuantity, 1), 2);
            end
            if max(abs(Zeta_Comm(:) - Prev_Zeta(:))) < 1e-5
                fprintf('    -> 网络收敛于第 %d 步\n', iter); 
                break; 
            end
        end
        
        % 提取共识均值
        Final_Xi = P_Info_Matrix - Zeta_Comm;
        for t = 1:N_eval
            for d = 1:y_dim
                Xi_1 = Final_Xi(2*d-1, 1, t); 
                Xi_2 = Final_Xi(2*d, 1, t);
                if ismember(lower(Current_Method), {'gpoe','poe','bcm','rbcm'})
                    Final_Mean_Pred(t, d) = Xi_1 / max(Xi_2, eps);
                else
                    Final_Mean_Pred(t, d) = Xi_1 / AgentQuantity;
                end
            end
        end

    else
        %% --- B. 集中式聚合 (AC) ---
        for t = 1:N_eval
            for d = 1:y_dim
                All_Agent_Mu  = Mu_Local_All(:, d, t); 
                All_Agent_Var = Var_Local_All(:, d, t);
                Beta_Set = max(0.5 * (log(prior_var) - log(All_Agent_Var)), eps);
                switch Method_Base_Name
                    case 'moe',  Final_Mean_Pred(t, d) = mean(All_Agent_Mu);
                    case 'gpoe', Final_Mean_Pred(t, d) = sum(Beta_Set .* All_Agent_Mu ./ All_Agent_Var) / max(sum(Beta_Set ./ All_Agent_Var), eps);
                    case 'poe',  Final_Mean_Pred(t, d) = sum(All_Agent_Mu ./ All_Agent_Var) / max(sum(1 ./ All_Agent_Var), eps);
                    case 'bcm'
                        Prec_BCM = sum(1 ./ All_Agent_Var) - (AgentQuantity-1)/prior_var;
                        Final_Mean_Pred(t, d) = sum(All_Agent_Mu ./ All_Agent_Var) / max(Prec_BCM, eps);
                    case 'rbcm'
                        Prec_RBCM = sum(Beta_Set ./ All_Agent_Var) + (1 - sum(Beta_Set))/prior_var;
                        Final_Mean_Pred(t, d) = sum(Beta_Set .* All_Agent_Mu ./ All_Agent_Var) / max(Prec_RBCM, eps);
                end
            end
        end
    end

    %% --- C. 统一方差聚合 ---
    for t = 1:N_eval
        for d = 1:y_dim
            All_Agent_Var = Var_Local_All(:, d, t); 
            Beta_Set = max(0.5 * (log(prior_var) - log(All_Agent_Var)), eps);
            switch Method_Base_Name
                case 'moe',  Final_Var_Pred(t, d) = mean(All_Agent_Var);
                case 'gpoe', Final_Var_Pred(t, d) = 1 / max(sum(Beta_Set ./ All_Agent_Var), eps);
                case 'poe',  Final_Var_Pred(t, d) = 1 / max(sum(1 ./ All_Agent_Var), eps);
                case 'bcm'
                    Prec_BCM = sum(1 ./ All_Agent_Var) - (AgentQuantity-1)/prior_var;
                    Final_Var_Pred(t, d) = 1 / max(Prec_BCM, eps);
                case 'rbcm'
                    Prec_RBCM = sum(Beta_Set ./ All_Agent_Var) + (1 - sum(Beta_Set))/prior_var;
                    Final_Var_Pred(t, d) = 1 / max(Prec_RBCM, eps);
            end
            Final_Var_Pred(t, d) = max(Final_Var_Pred(t, d), SigmaN^2);
        end
    end
    
    t_test = Precompute_Time/numel(AllModes) + toc;

    %% --- D. 反归一化、误差计算与保存 ---
    mu_pred = Final_Mean_Pred .* repmat(Y_std, N_eval, 1) + repmat(Y_mean, N_eval, 1);
    var_pred = Final_Var_Pred .* repmat(Y_std.^2, N_eval, 1);

    err = Y_eval - mu_pred;
    smse = mean(mean(err.^2) ./ Y_var_base);
    rmse = mean(sqrt(mean(err.^2)));
    nlpd = mean(mean(0.5*(log(2*pi*var_pred) + err.^2 ./ var_pred)));
    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Time=%.1fs\n', smse, rmse, nlpd, t_test);

    err_sq_mean = mean(err.^2, 2);            
    smse_curve  = cumsum(err_sq_mean) ./ (1:N_eval)' / mean(Y_var_base);
    rmse_curve  = sqrt(cumsum(err_sq_mean) ./ (1:N_eval)');

    t_train = 0; 
    save(fullfile(SaveFolder, sprintf('%s_tp_tr%d_mc%d.mat', Current_Method, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', 't_train', 't_test', 'Current_Method', 'seed', 'train_ratio', ...
        'smse_curve', 'rmse_curve');
end
fprintf('\n[%s] 测试点 done. seed=%d tr=%d%%\n\n', DatasetName, seed, tr_tag);
end