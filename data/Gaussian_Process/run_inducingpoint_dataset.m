function run_inducingpoint_dataset(DatasetName, CurrentMode, train_ratio, seed)
if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1; end
rng(seed);
fprintf('\n[诱导点] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

%% 1. 加载数据集
switch upper(DatasetName)
    case 'KIN40K'
        tr=load('KIN40K_train.mat'); te=load('KIN40K_test.mat'); hp=load('KIN40K_Hyperparameter.mat');
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
        tr=load(fullfile('SARCOS','SARCOS_train.mat')); te=load(fullfile('SARCOS','SARCOS_test.mat'));
        hp_raw=load(fullfile('SARCOS','SARCOS_GP_Hyperparameter.mat'));
        hp.SigmaF=mean(cell2mat(hp_raw.SigmaF_set)); hp.SigmaN=mean(cell2mat(hp_raw.SigmaN_set));
        hp.SigmaL=mean(cell2mat(hp_raw.SigmaL_set'),2);
        train_x=tr.sarcos_inv(:,1:21); train_y=tr.sarcos_inv(:,22:28);
        test_x=te.sarcos_inv_test(:,1:21); test_y=te.sarcos_inv_test(:,22:28);
    otherwise, error('未知数据集: %s',DatasetName);
end

%% 2. 归一化
if size(hp.SigmaL,1)>1&&size(hp.SigmaL,2)>1, hp.SigmaL=mean(hp.SigmaL,1); end
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
if ~(max(abs(X_mean))<1e-2&&max(abs(X_std-1))<1e-2)
    X_train=(X_train-X_mean)./X_std; X_test=(X_test-X_mean)./X_std; SigmaL=SigmaL./X_std(:);
end
Y_mean=mean(Y_train,1); Y_std=std(Y_train,0,1); Y_std(Y_std==0)=1;
if max(abs(Y_mean))<1e-2&&max(abs(Y_std-1))<1e-2
    Y_mean=zeros(1,size(Y_train,2)); Y_std=ones(1,size(Y_train,2));
else
    Y_train=(Y_train-Y_mean)./Y_std; SigmaF=SigmaF/mean(Y_std); SigmaN=SigmaN/mean(Y_std);
end
prior_var=SigmaF^2;
[N_train,x_dim]=size(X_train); y_dim=size(Y_train,2);
N_eval=min(3000,size(X_test,1));
X_eval=X_test(1:N_eval,:); Y_eval=Y_test(1:N_eval,:);
Y_var_base=var(Y_eval,0,1);
fprintf('Train=%d Test=%d x=%d y=%d\n',N_train,N_eval,x_dim,y_dim);

%% 3. 分布式参数
AgentQuantity=6; Kappa_P=10; t_step=0.01;
MaxDataPerAgent=min(floor(N_train/AgentQuantity),3000);
switch upper(DatasetName)
    case {'SARCOS','POL'}, NumInducingPoints=2500;
    otherwise,             NumInducingPoints=2000;
end
MultiAgentSystem=Manipulator_2D_2DoF_SetMASTopology(AgentQuantity,1);
L=MultiAgentSystem.Agent_Topology.LaplacianMatrix;
N_degree=sum(L<0,2); N_max=max(N_degree);
num_directed_neighbor_links=sum(N_degree);
num_undirected_edges=num_directed_neighbor_links/2;
% ET参数：Dimarogonas 2012 distributed event-triggered consensus
% 触发条件： ||e_i||^2 > [sigma_i * a * (1-a*|N_i|)/|N_i|] * ||z_i||^2
% 约束：0 < a < 1/|N_i|, 0 < sigma_i < 1。
% 如果触发太少：减小 sigma_i / sigma_i_ac；如果触发太多：增大它们。
a_param    = 0.5 / N_max;   % 自动满足 a < 1/max|N_i|
sigma_i    = 0.5;           % IP-DAC和IP-AC统一使用，便于对比
sigma_i_ac = 0.5;           % 和sigma_i相同，确保对比公平
fprintf('ET参数: a=%.4g  sigma_DAC=%.4g  sigma_AC=%.4g  max|N_i|=%d\n', ...
    a_param, sigma_i, sigma_i_ac, N_max);
%fprintf('Topology degree N_i = '); disp(N_degree(:)');
%fprintf('L row-sum max=%.2e, col-sum max=%.2e, symmetric=%d\n', ...
    %max(abs(sum(L,2))), max(abs(sum(L,1))), isequal(L,L'));


%% 4. 训练局部GP
t_start=tic;
LocalGP_set=cell(AgentQuantity,1);
for n=1:AgentQuantity
    idx=(n-1)*MaxDataPerAgent+1:min(n*MaxDataPerAgent,N_train);
    LocalGP_set{n}=LocalGP_MultiOutput(x_dim,y_dim,MaxDataPerAgent,SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)',Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8; LocalGP_set{n}.delta=0.01;
end
t_train_gp=toc(t_start);
fprintf('局部GP: %.4fs\n',t_train_gp);
idx_ind=randperm(N_train,NumInducingPoints);
InducingPoints_Coordinates=X_train(idx_ind,:)';

%% 5. 方法列表
dac_methods={'poe','gpoe','moe','bcm','rbcm'};
ac_methods={'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all'), AllModes=[dac_methods,ac_methods];
else, AllModes={lower(CurrentMode)}; end

%% 6. 批量预计算P矩阵
p_dim=2*y_dim;
P_poe=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_gpoe=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_moe=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_bcm=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_rbcm=zeros(p_dim,AgentQuantity,NumInducingPoints);
mu_ind=zeros(AgentQuantity,y_dim,NumInducingPoints);
var_ind=zeros(AgentQuantity,y_dim,NumInducingPoints);
fprintf('[预计算] %d个诱导点...\n',NumInducingPoints); tic;
for n=1:AgentQuantity
    [mu_n_all,var_n_all]=batch_predict_external(LocalGP_set{n},InducingPoints_Coordinates,SigmaN,SigmaF);
    for m=1:NumInducingPoints
        mu_n=mu_n_all(m,:)'; var_n=repmat(var_n_all(m),y_dim,1);
        mu_ind(n,:,m)=mu_n'; var_ind(n,:,m)=var_n';
        for d=1:y_dim
            vs=max(var_n(d),SigmaN^2);
            b=max(min(0.5*(log(prior_var)-log(vs)),10),eps);
            P_poe(2*d-1,n,m)=AgentQuantity*mu_n(d)/vs;
            P_poe(2*d,n,m)=AgentQuantity/vs;
            P_gpoe(2*d-1,n,m)=AgentQuantity*b*mu_n(d)/vs;
            P_gpoe(2*d,n,m)=AgentQuantity*b/vs;
            P_moe(2*d-1,n,m)=AgentQuantity*mu_n(d);
            P_moe(2*d,n,m)=AgentQuantity*(vs+mu_n(d)^2);
            P_bcm(2*d-1,n,m)=AgentQuantity*mu_n(d)/vs;
            P_bcm(2*d,n,m)=AgentQuantity/vs-(AgentQuantity-1)/prior_var;
            P_rbcm(2*d-1,n,m)=AgentQuantity*b*mu_n(d)/vs;
            P_rbcm(2*d,n,m)=AgentQuantity*b/vs+(1-AgentQuantity*b)/prior_var;
        end
    end
end
fprintf('预计算完成: %.2fs\n',toc);
SaveFolder=fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
tr_tag=round(train_ratio*100);

%% 7. 主循环
for mi=1:numel(AllModes)
    cur=AllModes{mi};
    fprintf('\n[%d/%d] %s\n',mi,numel(AllModes),cur);
    tic;
    trigger_count_per_agent=zeros(AgentQuantity,NumInducingPoints);

    base_method=strrep(lower(cur),'_ac','');
    if     strcmpi(base_method,'gpoe'), Pi=P_gpoe;
    elseif strcmpi(base_method,'poe'),  Pi=P_poe;
    elseif strcmpi(base_method,'bcm'),  Pi=P_bcm;
    elseif strcmpi(base_method,'rbcm'), Pi=P_rbcm;
    else,                               Pi=P_moe;
    end

    iter=0; max_iter=3000;

    % ============================================================
    % Convergence diagnostics
    % We only use one representative scalar component for the agent-state
    % trajectories, but the error curves below are averaged over all
    % p-dimensions and inducing points.
    %   DAC test 1: tracking error to the average reference input
    %   DAC test 2: inter-agent disagreement
    %   AC  test 1: error to the initial average value
    % ============================================================
    conv_r = 1;      % representative information-vector dimension
    conv_m = 1;      % representative inducing point

    conv_curve_dac = [];                 % backward-compatible field
    conv_curve_ac  = [];                 % backward-compatible field

    conv_dac_tracking_curve     = [];
    conv_dac_disagreement_curve = [];
    conv_ac_avg_error_curve     = [];
    conv_ac_disagreement_curve  = [];

    dac_state_hist = [];
    dac_ref_hist   = [];
    ac_state_hist  = [];
    ac_ref_hist    = [];

    if ismember(lower(cur),dac_methods)
        %% IP-DAC：严格区分“实时状态”和“广播快照”
        % 状态定义：
        %   Zeta              : DAC内部状态
        %   Xi_now = Pi-Zeta  : 每个agent当前实时consensus输出 x_i(t)
        %   Xi_last_trigger   : 上一次广播出去的快照 \hat{x}_i(t)
        %
        % 按事件触发论文的实现逻辑：
        %   动力学/控制输入使用快照 \hat{x}（零阶保持）；
        %   触发误差 e_i = \hat{x}_i - x_i 使用快照与实时值；
        %   z_i = sum_j (x_i - x_j) 使用实时值 x_i, x_j。
        Zeta=zeros(p_dim,AgentQuantity,NumInducingPoints);
        Xi_last_trigger=Pi;  % 初始时刻广播 x_i(0)=Pi_i

        Pi_ref   = mean(Pi, 2);
        Pi_scale = max(1, mean(Pi(:).^2)); %#ok<NASGU>

        conv_dac_tracking_curve     = zeros(max_iter,1);
        conv_dac_disagreement_curve = zeros(max_iter,1);
        dac_state_hist = zeros(max_iter,AgentQuantity);
        dac_ref_hist   = zeros(max_iter,1);
        conv_curve_dac = zeros(max_iter,1);

        while iter<max_iter
            iter=iter+1; Zeta_prev=Zeta;

            % 动力学向量化：Zeta_dot = Kappa_P * L * Xi_hat
            % permute: [p_dim x N x M] → [N x p_dim x M] → [N x p_dim*M]
            Xi_hat_perm = permute(Xi_last_trigger, [2 1 3]);  % [N x p_dim x M]
            Xi_hat_2d   = reshape(Xi_hat_perm, AgentQuantity, p_dim*NumInducingPoints);  % [N x p_dim*M]
            L_Xi_hat_2d = L * Xi_hat_2d;  % [N x p_dim*M]
            L_Xi_hat    = permute(reshape(L_Xi_hat_2d, AgentQuantity, p_dim, NumInducingPoints), [2 1 3]);  % [p_dim x N x M]
            Zeta = Zeta + t_step * Kappa_P * L_Xi_hat;

            % 当前实时consensus输出
            Xi_now = Pi - Zeta;

            % 收敛诊断
            dac_disagreement = abs(Xi_now - mean(Xi_now,2));
            conv_dac_tracking_curve(iter)     = mean(max(abs(Xi_now-Pi_ref),[],2),'all');
            conv_dac_disagreement_curve(iter) = mean(max(dac_disagreement,[],2),'all');
            dac_state_hist(iter,:) = squeeze(Xi_now(conv_r,:,conv_m));
            dac_ref_hist(iter)     = mean(squeeze(Pi(conv_r,:,conv_m)));
            conv_curve_dac(iter)   = conv_dac_disagreement_curve(iter);

            % ET触发条件（point-wise向量化，消除内层for m循环）
            % e_i = x̂_i - x_i（广播快照 - 实时状态）
            % z_i = Σ_{j∈N_i}(x_i - x_j)（实时状态的邻居差）
            for agent_i = 1:AgentQuantity
                neighbor_mask_i = (L(agent_i,:) < 0);
                N_i = N_degree(agent_i);
                coeff_i = sigma_i * a_param * (1 - a_param*N_i) / N_i;

                % e_tilde：[p_dim x 1 x M]
                e_tilde_all = Xi_last_trigger(:,agent_i,:) - Xi_now(:,agent_i,:);
                % z_i：对所有邻居求和，[p_dim x 1 x M]
                z_all = sum(Xi_now(:,agent_i,:) - Xi_now(:,neighbor_mask_i,:), 2);

                % 各诱导点的平方和：[1 x 1 x M] → squeeze成[M x 1]
                e_sq_all       = squeeze(sum(e_tilde_all.^2, 1));  % [M x 1]
                z_sq_all       = squeeze(sum(z_all.^2, 1));        % [M x 1]
                threshold_all  = coeff_i * z_sq_all;               % [M x 1]

                % 触发mask：[M x 1]
                trigger_mask = (e_sq_all > threshold_all) | (iter == 1);

                % 批量更新触发的诱导点
                Xi_last_trigger(:,agent_i,trigger_mask) = Xi_now(:,agent_i,trigger_mask);
                trigger_count_per_agent(agent_i,trigger_mask) = ...
                    trigger_count_per_agent(agent_i,trigger_mask) + 1;
            end

            if max(abs(Zeta(:)-Zeta_prev(:)))<1e-5, break; end
        end
        iter_converge=iter;

        conv_dac_tracking_curve     = conv_dac_tracking_curve(1:iter_converge);
        conv_dac_disagreement_curve = conv_dac_disagreement_curve(1:iter_converge);
        dac_state_hist = dac_state_hist(1:iter_converge,:);
        dac_ref_hist   = dac_ref_hist(1:iter_converge);
        conv_curve_dac = conv_curve_dac(1:iter_converge);

        % 导师口径：通信次数 = 触发次数
        % 主表报告 per-agent / per-inducing-point 的平均触发次数，
        % 与 Dimarogonas paper 中“每个 agent 的通信/更新次数”口径一致。
        % ============================================================
        % Communication metric for IP-DAC
        % Official metric: mean event-trigger count per agent per inducing point
        % This follows the paper-style interpretation: communication count = trigger count.
        % ============================================================

        trigger_per_agent_point = mean(trigger_count_per_agent, 2);  % [AgentQuantity x 1]
        comm_train = mean(trigger_count_per_agent(:));               % scalar, official Comm_Train
        comm_test  = 0;

        % Optional diagnostic variables, not used as the main metric
        trigger_total_per_point = sum(trigger_count_per_agent, 1);   % [1 x NumInducingPoints]
        trigger_total_mean_per_point = mean(trigger_total_per_point);

        fprintf('  [IP-DAC] 收敛步数:%d  平均触发次数/agent/point:%.1f\n', ...
            iter_converge, comm_train);

        %fprintf('    mean trigger per agent/point: ');
        %fprintf('%.4f  ', trigger_per_agent_point);
        %fprintf('\n');

%fprintf('    diagnostic only: total trigger per point mean = %.1f\n', ...
    %trigger_total_mean_per_point);
        Xi_final=Pi-Zeta;
    else
        %% IP-AC：静态平均一致性，也严格区分实时状态和广播快照
        %   Xi              : 当前实时状态 x_i(t)
        %   Xi_last_trigger : 上一次广播快照 \hat{x}_i(t)
        % 动力学使用快照： xdot = -Kappa_P * L * xhat
        % 触发条件中：
        %   e_i = xhat_i - x_i 用快照与实时值；
        %   z_i = sum_j(x_i - x_j) 用实时值。
        Xi=Pi;
        Xi_last_trigger=Pi;

        % AC reference is the initial average. This is fixed by definition:
        % x_i(k) should converge to mean_i x_i(0).
        Xi_ref0 = mean(Pi, 2);  % [p_dim x 1 x NumInducingPoints]

        conv_ac_avg_error_curve    = zeros(max_iter,1);
        conv_ac_disagreement_curve = zeros(max_iter,1);
        ac_state_hist = zeros(max_iter,AgentQuantity);
        ac_ref_hist   = zeros(max_iter,1);
        conv_curve_ac = zeros(max_iter,1);

        while iter<max_iter
            iter=iter+1; Xi_prev=Xi;

            % 动力学向量化：xdot = -Kappa_P * L * xhat
            Xi_hat_perm = permute(Xi_last_trigger, [2 1 3]);
            Xi_hat_2d   = reshape(Xi_hat_perm, AgentQuantity, p_dim*NumInducingPoints);
            L_Xi_hat_2d = L * Xi_hat_2d;
            L_Xi_hat    = permute(reshape(L_Xi_hat_2d, AgentQuantity, p_dim, NumInducingPoints), [2 1 3]);
            Xi = Xi - t_step * Kappa_P * L_Xi_hat;

            % AC convergence diagnostic
            % AC should converge to the initial average value, not to zero.
            ac_avg_err      = abs(Xi - Xi_ref0);
            ac_disagreement = abs(Xi - mean(Xi,2));

            conv_ac_avg_error_curve(iter)    = mean(max(ac_avg_err, [], 2), 'all');
            conv_ac_disagreement_curve(iter) = mean(max(ac_disagreement, [], 2), 'all');

            ac_state_hist(iter,:) = squeeze(Xi(conv_r,:,conv_m));
            ac_ref_hist(iter)     = mean(squeeze(Pi(conv_r,:,conv_m)));
            conv_curve_ac(iter)   = conv_ac_avg_error_curve(iter);

            % ET触发条件（point-wise向量化）
            % e_i = x̂_i - x_i，z_i = Σ_{j∈N_i}(x_i - x_j)（实时状态）
            for agent_i = 1:AgentQuantity
                neighbor_mask_i = (L(agent_i,:) < 0);
                N_i = N_degree(agent_i);
                coeff_i = sigma_i_ac * a_param * (1 - a_param*N_i) / N_i;

                e_tilde_all = Xi_last_trigger(:,agent_i,:) - Xi(:,agent_i,:);
                z_all = sum(Xi(:,agent_i,:) - Xi(:,neighbor_mask_i,:), 2);

                e_sq_all      = squeeze(sum(e_tilde_all.^2, 1));
                z_sq_all      = squeeze(sum(z_all.^2, 1));
                threshold_all = coeff_i * z_sq_all;

                trigger_mask = (e_sq_all > threshold_all) | (iter == 1);

                Xi_last_trigger(:,agent_i,trigger_mask) = Xi(:,agent_i,trigger_mask);
                trigger_count_per_agent(agent_i,trigger_mask) = ...
                    trigger_count_per_agent(agent_i,trigger_mask) + 1;
            end

            if max(abs(Xi(:)-Xi_prev(:)))<1e-5, break; end
        end
        iter_converge=iter;

        conv_ac_avg_error_curve    = conv_ac_avg_error_curve(1:iter_converge);
        conv_ac_disagreement_curve = conv_ac_disagreement_curve(1:iter_converge);
        ac_state_hist = ac_state_hist(1:iter_converge,:);
        ac_ref_hist   = ac_ref_hist(1:iter_converge);
        conv_curve_ac = conv_curve_ac(1:iter_converge);

        conv_curve_dac = [];  % AC分支，DAC曲线为空

        % 导师口径：通信次数 = 触发次数
        % AC分支同样报告 per-agent / per-point 的平均触发次数；
        % total per point 只作为诊断信息保存和打印。
        % ============================================================
        % Communication metric for IP-AC
        % Official metric: mean update/trigger count per agent per inducing point
        % ============================================================

        trigger_per_agent_point = mean(trigger_count_per_agent, 2);  % [AgentQuantity x 1]
        comm_train = mean(trigger_count_per_agent(:));               % scalar, official Comm_Train
        comm_test  = 0;

        trigger_total_per_point = sum(trigger_count_per_agent, 1);   % diagnostic only
        trigger_total_mean_per_point = mean(trigger_total_per_point);

        fprintf('  [IP-AC] 收敛步数:%d  平均更新次数/agent/point:%.1f\n', ...
            iter_converge, comm_train);

        %fprintf('    mean update per agent/point: ');
        %fprintf('%.4f  ', trigger_per_agent_point);
        %fprintf('\n');

        %fprintf('    diagnostic only: total update per point mean = %.1f\n', ...
            %trigger_total_mean_per_point);
        Xi_final=Xi;
    end  %% if/else结束

    %% 提取phi（DAC/AC收敛后直接用Xi_final）
    phi=zeros(y_dim,NumInducingPoints);
    for d=1:y_dim
        xi1=squeeze(Xi_final(2*d-1,1,:))';
        xi2=squeeze(Xi_final(2*d,  1,:))';
        if ismember(base_method,{'gpoe','poe','bcm','rbcm'})
            % 防止分母除零，但保留 denominator 的符号。
            % 不再使用 abs(xi2)，否则 BCM/rBCM 的先验修正项可能被改变符号。
            den = xi2;
            small_den = abs(den) < 1e-4;
            den(small_den & den >= 0) = 1e-4;
            den(small_den & den <  0) = -1e-4;
            phi(d,:) = xi1 ./ den;
        else
            phi(d,:)=xi1/AgentQuantity;
        end
    end

    %% 训练MaskedGP并预测
    MaskedGP=LocalGP_MultiOutput(x_dim,y_dim,NumInducingPoints,1e-6,SigmaF,SigmaL);
    MaskedGP.add_Alldata(InducingPoints_Coordinates,phi);
    t_train=t_train_gp+toc; tic;
    Num_Inducing=MaskedGP.DataQuantity;
    Alpha_Vec=MaskedGP.alpha(1:Num_Inducing,:);
    Cholesky_L=MaskedGP.L(1:Num_Inducing,1:Num_Inducing);
    K_star=MaskedGP.kernel(MaskedGP.X(:,1:Num_Inducing),X_eval');
    mu_normalized=(Alpha_Vec'*K_star)';
    V_matrix=Cholesky_L\K_star;
    var_normalized=max(SigmaF^2-sum(V_matrix.^2,1)',SigmaN^2);
    mu_pred=mu_normalized.*repmat(Y_std,N_eval,1)+repmat(Y_mean,N_eval,1);
    var_pred=repmat(var_normalized,1,y_dim).*repmat(Y_std.^2,N_eval,1);
    t_test=toc;

    err=Y_eval-mu_pred;
    % 防弹衣3：Y_var_base为0时（静态常数测试集）用1代替
    Y_var_base_safe = max(Y_var_base, 1e-10);
    smse=mean(mean(err.^2)./Y_var_base_safe);
    rmse=mean(sqrt(mean(err.^2)));
    nlpd=mean(mean(0.5*(log(2*pi*var_pred)+err.^2./var_pred)));
    t_train_per_point=(t_train/N_train)*1000;
    t_test_per_point=(t_test/N_eval)*1000;
    fprintf('  SMSE=%.4f RMSE=%.4f NLPD=%.4f Train:%.2fms/pt Test:%.2fms/pt\n',...
        smse,rmse,nlpd,t_train_per_point,t_test_per_point);

    % Compatibility variables expected by run_mc_dataset / plotting scripts
    t_train_total = t_train;
    t_test_total  = t_test;
    current_method = cur;

    err_sq_mean=mean(err.^2,2);
    smse_curve=cumsum(err_sq_mean)./(1:N_eval)'/mean(Y_var_base_safe);
    rmse_curve=sqrt(cumsum(err_sq_mean)./(1:N_eval)');
    event_count_mean=comm_train;

    save(fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat', current_method, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', ...
        't_train_total', 't_test_total', ...
        't_train_per_point', 't_test_per_point', ...
        'comm_train', 'comm_test', 'iter_converge', ...
        'event_count_mean', ...
        'current_method', 'seed', 'train_ratio', ...
        'smse_curve', 'rmse_curve', ...
        'trigger_count_per_agent', ...
        'trigger_per_agent_point', ...
        'trigger_total_per_point', ...
        'trigger_total_mean_per_point', ...
        'conv_curve_dac', 'conv_curve_ac', ...
        'conv_dac_tracking_curve', 'conv_dac_disagreement_curve', ...
        'conv_ac_avg_error_curve', 'conv_ac_disagreement_curve', ...
        'dac_state_hist', 'dac_ref_hist', ...
        'ac_state_hist', 'ac_ref_hist');
end
end