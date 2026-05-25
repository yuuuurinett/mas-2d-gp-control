function run_testpoint_dataset(DatasetName, CurrentMode, train_ratio, seed)
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
t_start=tic;
LocalGP_set=cell(AgentQuantity,1);
for n=1:AgentQuantity
    idx=(n-1)*MaxDataPerAgent+1 : min(n*MaxDataPerAgent,N_train);
    LocalGP_set{n}=LocalGP_MultiOutput(x_dim,y_dim,MaxDataPerAgent,SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)',Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8; LocalGP_set{n}.delta=0.01;
end
t_train=toc(t_start);
fprintf('局部GP训练完成: %.4fs\n',t_train);

%% 5. 方法字典配置
dac_methods={'moe','gpoe','poe','bcm','rbcm'};
ac_methods={'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all'), AllModes=[dac_methods,ac_methods];
else, AllModes={lower(CurrentMode)}; end

SaveFolder=fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
p_dim=2*y_dim; tr_tag=round(train_ratio*100);

%% 6. 统一预计算阶段 (向量化极速版)
fprintf('\n[核心优化] 正在预计算 %d 个测试点的局部推断 (6个专家并算)...\n', N_eval);
tic;
Mu_Local_All  = zeros(AgentQuantity, y_dim, N_eval);
Var_Local_All = zeros(AgentQuantity, y_dim, N_eval);

% 核心优化：将测试点直接转置为矩阵 [x_dim, N_eval]
X_eval_matrix = X_eval';

for n = 1:AgentQuantity
    try
        % 尝试极速批量预测 (秒级完成)
        [mn_batch, vn_batch] = LocalGP_set{n}.predict(X_eval_matrix);
        
        % 维度自适应装填
        if size(mn_batch, 2) == N_eval
            for t = 1:N_eval
                Mu_Local_All(n, :, t)  = mn_batch(:, t)';
                Var_Local_All(n, :, t) = vn_batch(:, t)';
            end
        else
            for t = 1:N_eval
                Mu_Local_All(n, :, t)  = mn_batch(t, :);
                Var_Local_All(n, :, t) = vn_batch(t, :);
            end
        end
    catch
        % 如果类方法不支持矩阵输入，安全回退到原版循环
        for t = 1:N_eval
            [mn, vn] = LocalGP_set{n}.predict(X_eval_matrix(:, t));
            Mu_Local_All(n, :, t)  = mn';
            Var_Local_All(n, :, t) = vn';
        end
    end
end
Precompute_Time = toc;
fprintf('预计算完成！总计耗时: %.2f 秒 (每个点平均 %.4f 秒)\n', Precompute_Time, Precompute_Time/N_eval);
%% 7. 聚合方法主循环 (共用预计算结果)
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
        % =========================================================================
        % [新增] 基于 IEEE TAC 2012 的分布式事件触发机制 (Distributed ET)
        % =========================================================================
        % Zeta   : 测试阶段本地连续演化的真实状态 (对应论文 x_i(t))
        % Zeta_k : 上一次触发时保存/广播的状态快照 (对应论文 x_i(t_k))
        Zeta   = zeros(p_dim, AgentQuantity, N_eval);
        Zeta_k = zeros(p_dim, AgentQuantity, N_eval); 
        
        trigger_count_set = zeros(AgentQuantity, 1); % 记录测试阶段触发次数
        sigma_i = 0.1;  % 触发系数
        a_param = 0.2;  % 论文参数 a
        
        comm_train = 0; % TP 架构离线训练无通信
        
        for iter = 1:3000
            Prev_Zeta = Zeta;
            
            % ---------------------------------------------------------
            % 1. 计算网络一致性驱动力 (邻居间仅通过触发状态 Zeta_k 交互)
            % 注意: TP 架构这里的目标矩阵叫 P_Info_Matrix
            % ---------------------------------------------------------
            diff_k = P_Info_Matrix - Zeta_k; 
            
            % 连续动力演化 (相当于 dx = -Lx)
            for n = 1:AgentQuantity
                Zeta(:, n, :) = Zeta(:, n, :) + t_step * Kappa_P * ...
                    sum(diff_k .* reshape(L(n, :), 1, AgentQuantity, 1), 2);
            end
            
            % ---------------------------------------------------------
            % 2. Distributed Event-Trigger 判定 (遵循实验室 rho_i > rho_bar_i 风格)
            % ---------------------------------------------------------
            for n = 1:AgentQuantity
                % 测量误差 e_i
                e_i = Zeta_k(:,n,:) - Zeta(:,n,:);
                
                % 邻居间相对误差 z_i
                z_i = sum((Zeta_k(:,n,:) - Zeta_k) .* reshape(L(n,:)<0, 1, AgentQuantity, 1), 2);
                
                % 符号完全对齐导师习惯
                rho_i     = max(sum(e_i.^2, 1)); % ||e_i||^2
                norm_z_sq = max(sum(z_i.^2, 1)); % ||z_i||^2
                
                N_i = sum(L(n,:) < 0); % 邻居数量 |N_i|
                rho_bar_i = (sigma_i * a_param * (1 - a_param * N_i) / N_i) * norm_z_sq; 
                
                % 触发判定条件 
                if rho_i > rho_bar_i || iter == 1
                    Zeta_k(:,n,:) = Zeta(:,n,:);       % 触发更新并广播
                    trigger_count_set(n) = trigger_count_set(n) + 1; 
                end
            end
            
            % ---------------------------------------------------------
            % 3. 停止准则 (基于底层真实连续状态)
            % ---------------------------------------------------------
            if max(abs(Zeta(:) - Prev_Zeta(:))) < 1e-5
                fprintf('    -> [Testing ET] 网络收敛于第 %d 步\n', iter);
                break; 
            end
        end % <--- 这是 for iter = 1:3000 的结束括号
        
        % =========================================================
        % [核心修复] 将 comm_test 移至循环体外，防止未收敛时保存报错
        % =========================================================
        comm_test = mean(trigger_count_set); 
        fprintf('    -> 平均实际物理通信: %.1f 次\n', comm_test);
        
        % 提取共识均值 (用收敛后的真实状态 Zeta 替换原来的 Zeta_Comm)
        Final_Xi = P_Info_Matrix - Zeta;
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
        comm_train = 0;
        comm_test  = 1;
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

    %% --- 反归一化、误差计算与单点时间换算 ---
    % 你的原有误差计算代码保留...
    err = Y_eval - mu_pred;
    smse = mean(mean(err.^2) ./ Y_var_base);
    rmse = mean(sqrt(mean(err.^2)));
    nlpd = mean(mean(0.5*(log(2*pi*var_pred) + err.^2 ./ var_pred)));

    % [核心新增] 计算单点耗时 (ms/point)
    t_train_per_point = (t_train / N_train) * 1000;
    t_test_per_point  = (t_test / N_eval) * 1000;

    % 修改打印输出，展示 ms/pt
    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train: %.2f ms/pt  Test: %.2f ms/pt\n', ...
        smse, rmse, nlpd, t_train_per_point, t_test_per_point);

    err_sq_mean = mean(err.^2, 2);            
    smse_curve  = cumsum(err_sq_mean) ./ (1:N_eval)' / mean(Y_var_base);
    rmse_curve  = sqrt(cumsum(err_sq_mean) ./ (1:N_eval)');

    save(fullfile(SaveFolder, sprintf('%s_tp_tr%d_mc%d.mat', Current_Method, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', 't_train', 't_test', ...
        't_train_per_point', 't_test_per_point', ... 
        'comm_train', 'comm_test', ...
        'Current_Method', 'seed', 'train_ratio', 'smse_curve', 'rmse_curve');
end