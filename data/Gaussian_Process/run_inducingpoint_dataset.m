function run_inducingpoint_dataset(DatasetName, CurrentMode, train_ratio, seed)

if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1;          end
rng(seed);

fprintf('\n[诱导点] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

%% 1. 加载数据集
switch upper(DatasetName)
    case 'KIN40K'
        tr = load('KIN40K_train.mat');          te = load('KIN40K_test.mat');
        hp = load('KIN40K_Hyperparameter.mat');
        train_x = tr.x;  train_y = tr.y;  test_x = te.xtest;  test_y = te.ytest;
    case 'POL'
        tr = load(fullfile('POL','POL_train.mat'));  te = load(fullfile('POL','POL_test.mat'));
        hp = load(fullfile('POL','POL_Hyperparameter.mat'));
        train_x = tr.x;  train_y = tr.y;  test_x = te.xtest;  test_y = te.ytest;
    case 'PUMADYN32NM'
        tr = load(fullfile('PUMADYN32NM','PUMADYN32NM_train.mat'));
        te = load(fullfile('PUMADYN32NM','PUMADYN32NM_test.mat'));
        hp = load(fullfile('PUMADYN32NM','PUMADYN32NM_Hyperparameter.mat'));
        train_x = tr.x;  train_y = tr.y;  test_x = te.xtest;  test_y = te.ytest;
    case 'SARCOS'
        tr = load(fullfile('SARCOS','SARCOS_train.mat'));
        te = load(fullfile('SARCOS','SARCOS_test.mat'));
        hp_raw = load(fullfile('SARCOS','SARCOS_GP_Hyperparameter.mat'));
        hp.SigmaF = mean(cell2mat(hp_raw.SigmaF_set));
        hp.SigmaN = mean(cell2mat(hp_raw.SigmaN_set));
        hp.SigmaL = mean(cell2mat(hp_raw.SigmaL_set'), 2);
        train_x = tr.sarcos_inv(:,1:21);       train_y = tr.sarcos_inv(:,22:28);
        test_x  = te.sarcos_inv_test(:,1:21);  test_y  = te.sarcos_inv_test(:,22:28);
    otherwise,  error('未知数据集: %s', DatasetName);
end

%% 2. 归一化
if size(hp.SigmaL,1)>1 && size(hp.SigmaL,2)>1, hp.SigmaL = mean(hp.SigmaL,1); end
SigmaL = hp.SigmaL(:);
if numel(hp.SigmaF)>1, hp.SigmaF = mean(hp.SigmaF); end
if numel(hp.SigmaN)>1, hp.SigmaN = mean(hp.SigmaN); end
SigmaF = hp.SigmaF;  SigmaN = hp.SigmaN;

n_train = round(size(train_x,1) * train_ratio);
idx_tr  = randperm(size(train_x,1), n_train);
X_train = train_x(idx_tr,:);  Y_train = train_y(idx_tr,:);
idx_te  = randperm(size(test_x,1));
X_test  = test_x(idx_te,:);   Y_test  = test_y(idx_te,:);

X_mean = mean(X_train,1);  X_std = std(X_train,0,1);  X_std(X_std==0)=1;
if ~(max(abs(X_mean))<1e-2 && max(abs(X_std-1))<1e-2)
    X_train=(X_train-X_mean)./X_std;
    X_test=(X_test-X_mean)./X_std;
    SigmaL  = SigmaL ./ X_std(:);
end
Y_mean = mean(Y_train,1);  Y_std = std(Y_train,0,1);  Y_std(Y_std==0)=1;
if max(abs(Y_mean))<1e-2 && max(abs(Y_std-1))<1e-2
    Y_mean=zeros(1,size(Y_train,2)); Y_std=ones(1,size(Y_train,2));
else
    Y_train=(Y_train-Y_mean)./Y_std;
    SigmaF=SigmaF/mean(Y_std);  SigmaN=SigmaN/mean(Y_std);
end
prior_var = SigmaF^2;

[N_train,x_dim] = size(X_train);  y_dim = size(Y_train,2);
N_eval = min(3000,size(X_test,1));
X_eval = X_test(1:N_eval,:);  Y_eval = Y_test(1:N_eval,:);
Y_var_base = var(Y_eval,0,1);
fprintf('Train=%d  Test=%d  x=%d  y=%d\n', N_train,N_eval,x_dim,y_dim);

%% 3. 分布式参数
AgentQuantity = 6;  Kappa_P = 10;  t_step = 0.01;
MaxDataPerAgent = min(floor(N_train/AgentQuantity), 3000);
switch upper(DatasetName)
    case {'SARCOS','POL'},  NumInducingPoints = 2500;
    otherwise,              NumInducingPoints = 2000;
end
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity,1);
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

N_degree = sum(L < 0, 2);              % |N_i|: 每个 agent 的邻居数量
N_max    = max(N_degree);
num_directed_neighbor_links = sum(N_degree);  % 一轮全图有向邻居交换次数
num_undirected_edges = num_directed_neighbor_links / 2;
sigma_i  = 0.5;
a_param  = 0.5 / N_max;
fprintf('ET 参数: sigma=%.2f  a=%.4f  (N_max=%d)\n', sigma_i, a_param, N_max);

%% 4. 训练局部 GP
t_start = tic;
LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    idx = (n-1)*MaxDataPerAgent+1 : min(n*MaxDataPerAgent,N_train);
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim,y_dim,MaxDataPerAgent,SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)', Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8;  
    LocalGP_set{n}.delta=0.01;
end
t_train_gp = toc(t_start);
fprintf('局部GP: %.4fs\n', t_train_gp);

idx_ind = randperm(N_train, NumInducingPoints);
InducingPoints_Coordinates = X_train(idx_ind,:)';

%% 5. 方法
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all'), AllModes=[dac_methods,ac_methods];
else, AllModes={lower(CurrentMode)}; end

%% 6. 批量预计算 P 矩阵
p_dim = 2*y_dim;
P_poe =zeros(p_dim,AgentQuantity,NumInducingPoints);
P_gpoe=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_moe =zeros(p_dim,AgentQuantity,NumInducingPoints);
P_bcm =zeros(p_dim,AgentQuantity,NumInducingPoints);
P_rbcm=zeros(p_dim,AgentQuantity,NumInducingPoints);
mu_ind =zeros(AgentQuantity,y_dim,NumInducingPoints);
var_ind=zeros(AgentQuantity,y_dim,NumInducingPoints);

fprintf('[批量预计算] 正在预计算 %d 个诱导点...\n', NumInducingPoints);
tic;
for n = 1:AgentQuantity
    [mu_n_all, var_n_all] = batch_predict_external( ...
        LocalGP_set{n}, InducingPoints_Coordinates, SigmaN, SigmaF);
    for m = 1:NumInducingPoints
        mu_n  = mu_n_all(m,:)';
        var_n = repmat(var_n_all(m), y_dim, 1);
        mu_ind(n,:,m)  = mu_n';
        var_ind(n,:,m) = var_n';
        for d = 1:y_dim
            var_n_safe = max(var_n(d), SigmaN^2);
            b = max(min(0.5*(log(prior_var)-log(var_n_safe)), 10), eps);
            P_poe(2*d-1,n,m)  = AgentQuantity*mu_n(d)/var_n_safe;
            P_poe(2*d,  n,m)  = AgentQuantity/var_n_safe;

            P_gpoe(2*d-1,n,m) = AgentQuantity*b*mu_n(d)/var_n_safe;
            P_gpoe(2*d,  n,m) = AgentQuantity*b/var_n_safe;

            P_moe(2*d-1,n,m)  = AgentQuantity*mu_n(d);
            P_moe(2*d,  n,m)  = AgentQuantity*(var_n_safe+mu_n(d)^2);

            P_bcm(2*d-1,n,m)  = AgentQuantity*mu_n(d)/var_n_safe;
            P_bcm(2*d,  n,m)  = AgentQuantity/var_n_safe-(AgentQuantity-1)/prior_var;
            
            P_rbcm(2*d-1,n,m) = AgentQuantity*b*mu_n(d)/var_n_safe;
            P_rbcm(2*d,  n,m) = AgentQuantity*b/var_n_safe+(1-AgentQuantity*b)/prior_var;
        end
    end
end
fprintf('批量诱导点预计算完成！耗时: %.2fs\n', toc);

SaveFolder = fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
tr_tag = round(train_ratio*100);

%% 7. 各方法主循环
for mi = 1:numel(AllModes)
    cur = AllModes{mi};
    fprintf('\n[%d/%d] %s\n', mi, numel(AllModes), cur);
    tic;

    % 默认诊断变量，保证 DAC/AC 保存字段一致。
    event_count_mean = NaN;
    trigger_count_per_agent = [];

    switch lower(cur)
        case dac_methods
            if     strcmpi(cur,'gpoe'), Pi=P_gpoe;
            elseif strcmpi(cur,'poe'),  Pi=P_poe;
            elseif strcmpi(cur,'bcm'),  Pi=P_bcm;
            elseif strcmpi(cur,'rbcm'), Pi=P_rbcm;
            else,                       Pi=P_moe;
            end

            Zeta   = zeros(p_dim, AgentQuantity, NumInducingPoints);
            Zeta_k = zeros(p_dim, AgentQuantity, NumInducingPoints);

            % Point-wise 触发计数：
            % trigger_count_per_agent(n,m) = agent n 在 inducing point m 上触发/广播了多少次。
            % 这样最后可以统计“平均每个 inducing point 的邻居交换次数”。
            trigger_count_per_agent = zeros(AgentQuantity, NumInducingPoints);
            iter = 0;
            max_iter = 3000;  % 安全上限，防止死循环

            while iter < max_iter
                iter = iter + 1;
                Zeta_prev = Zeta;

                % 动力学更新（Pi - Zeta），保证收敛性
                diff = Pi - Zeta;
                for n = 1:AgentQuantity
                    Zeta(:,n,:) = Zeta(:,n,:) + t_step * Kappa_P * ...
                        sum(diff .* reshape(L(n,:), 1, AgentQuantity, 1), 2);
                end

                % ET 触发判定：point-wise 触发，每个 inducing point 独立判断。
                % 这样通信次数统计的是“平均每个 inducing point 上发生了多少次邻居交换”，
                % 更接近导师说的 neighbor exchange count，也避免把 2000 个点全量合成一个标量。
                for n = 1:AgentQuantity
                    % error_vec: [p_dim x 1 x NumInducingPoints]
                    % 表示 agent n 当前真实状态 Zeta 与上一次广播状态 Zeta_k 的差。
                    error_vec = Zeta_k(:,n,:) - Zeta(:,n,:);

                    % z_i: [p_dim x 1 x NumInducingPoints]
                    % 局部不一致度，z_i = sum_{j in N_i}(Zeta_i - Zeta_j)。
                    z_i = sum((Zeta(:,n,:) - Zeta) .* reshape(L(n,:)<0, 1, AgentQuantity, 1), 2);

                    % 对每个 inducing point，仅在 p_dim 维度上求平方和。
                    % squeeze 后得到 [NumInducingPoints x 1]。
                    norm_e_sq = squeeze(sum(error_vec.^2, 1));
                    norm_z_sq = squeeze(sum(z_i.^2, 1));

                    % 经典分布式 ET 阈值：e_i^2 <= coeff_i * z_i^2 + delta。
                    N_i = N_degree(n);
                    delta_bound = 1e-6;
                    coeff_i = sigma_i * a_param * (1 - a_param * N_i) / N_i;
                    rho_bar_vec = coeff_i * norm_z_sq + delta_bound;

                    % 每个 inducing point 独立触发。
                    trigger_idx = find(norm_e_sq > rho_bar_vec);

                    % 第一次迭代强制广播所有 inducing points，完成初始化。
                    if iter == 1
                        trigger_idx = 1:NumInducingPoints;
                    end

                    if ~isempty(trigger_idx)
                        Zeta_k(:, n, trigger_idx) = Zeta(:, n, trigger_idx);
                        trigger_count_per_agent(n, trigger_idx) = trigger_count_per_agent(n, trigger_idx) + 1;
                    end
                end           

                % 收敛判定
                if max(abs(Zeta(:) - Zeta_prev(:))) < 1e-5
                    break;
                end
            end

            iter_converge = iter;  % while 结束后直接就是收敛步数

            % 通信统计口径：导师说的“和邻居交换的次数”。
            % 对每个 inducing point m：sum_n trigger_count(n,m) * |N_n|
            % 最后对所有 inducing points 取平均，得到平均每个 inducing point 的邻居交换次数。
            event_count_mean = mean(sum(trigger_count_per_agent, 1));
            comm_train = mean(sum(trigger_count_per_agent .* N_degree(:), 1));
            comm_test = 0;
            fprintf('  [Training ET] 收敛步数: %d  平均Event: %.1f次/point  平均邻居交换: %.1f次/point\n', ...
                iter_converge, event_count_mean, comm_train);

            Xi = Pi - Zeta;
            phi=zeros(y_dim,NumInducingPoints);
            for d=1:y_dim
                xi1=squeeze(Xi(2*d-1,1,:))'; xi2=squeeze(Xi(2*d,1,:))';
                if ismember(lower(cur),{'gpoe','poe','bcm','rbcm'})
                    phi(d,:)=xi1./max(xi2,eps);
                else
                    phi(d,:)=xi1/AgentQuantity;
                end
            end

            case ac_methods
            % IP-AC：在诱导点上直接一次性全局聚合，不需要迭代
            if     strcmpi(cur,'gpoe_ac'), Pi=P_gpoe;
            elseif strcmpi(cur,'poe_ac'),  Pi=P_poe;
            elseif strcmpi(cur,'bcm_ac'),  Pi=P_bcm;
            elseif strcmpi(cur,'rbcm_ac'), Pi=P_rbcm;
            else,                          Pi=P_moe;
            end

            base_method = strrep(lower(cur),'_ac','');
            phi = zeros(y_dim, NumInducingPoints);
            for m = 1:NumInducingPoints
                for d = 1:y_dim
                    mu_all  = squeeze(mu_ind(:,d,m));
                    var_all = squeeze(var_ind(:,d,m));
                    beta_all = max(0.5*(log(prior_var)-log(var_all)), eps);
                    switch base_method
                        case 'moe'
                            phi(d,m) = mean(mu_all);
                        case 'gpoe'
                            prec = sum(beta_all./var_all);
                            phi(d,m) = sum(beta_all.*mu_all./var_all)/max(prec,eps);
                        case 'poe'
                            prec = sum(1./var_all);
                            phi(d,m) = sum(mu_all./var_all)/max(prec,eps);
                        case 'bcm'
                            prec = sum(1./var_all)-(AgentQuantity-1)/prior_var;
                            phi(d,m) = sum(mu_all./var_all)/max(prec,eps);
                        case 'rbcm'
                            prec = sum(beta_all./var_all)+(1-sum(beta_all))/prior_var;
                            phi(d,m) = sum(beta_all.*mu_all./var_all)/max(prec,eps);
                    end
                end
            end

            comm_train    = 1;   % 一次性全局聚合，通信次数=1
            comm_test     = 0;
            iter_converge = 1;
            event_count_mean = 1;
            trigger_count_per_agent = zeros(AgentQuantity, NumInducingPoints);
            fprintf('  [IP-AC] 直接聚合完成，comm_train=1\n');

    
    end

    MaskedGP=LocalGP_MultiOutput(x_dim,y_dim,NumInducingPoints,1e-6,SigmaF,SigmaL);
    MaskedGP.add_Alldata(InducingPoints_Coordinates,phi);
    t_train = t_train_gp + toc;

    tic;
    Num_Inducing = MaskedGP.DataQuantity;
    Alpha_Vec    = MaskedGP.alpha(1:Num_Inducing, :);
    Cholesky_L   = MaskedGP.L(1:Num_Inducing, 1:Num_Inducing);
    K_star = MaskedGP.kernel(MaskedGP.X(:, 1:Num_Inducing), X_eval');
    mu_normalized = (Alpha_Vec' * K_star)';
    V_matrix = Cholesky_L \ K_star;
    var_normalized = max(SigmaF^2 - sum(V_matrix.^2, 1)', SigmaN^2);
    mu_pred  = mu_normalized .* repmat(Y_std, N_eval, 1) + repmat(Y_mean, N_eval, 1);
    var_pred = repmat(var_normalized, 1, y_dim) .* repmat(Y_std.^2, N_eval, 1);
    t_test = toc;

    err = Y_eval - mu_pred;
    smse = mean(mean(err.^2) ./ Y_var_base);
    rmse = mean(sqrt(mean(err.^2)));
    nlpd = mean(mean(0.5*(log(2*pi*var_pred) + err.^2 ./ var_pred)));
    t_train_per_point = (t_train / N_train) * 1000;
    t_test_per_point  = (t_test / N_eval) * 1000;

    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train:%.2fms/pt  Test:%.2fms/pt\n', ...
        smse, rmse, nlpd, t_train_per_point, t_test_per_point);

    err_sq_mean = mean(err.^2, 2);
    smse_curve  = cumsum(err_sq_mean) ./ (1:N_eval)' / mean(Y_var_base);
    rmse_curve  = sqrt(cumsum(err_sq_mean) ./ (1:N_eval)');

    save(fullfile(SaveFolder, sprintf('%s_tr%d_mc%d.mat', cur, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', 't_train', 't_test', ...
        't_train_per_point', 't_test_per_point', ...
        'comm_train', 'comm_test', 'iter_converge', ...
        'event_count_mean', 'trigger_count_per_agent', 'N_degree', ...
        'num_directed_neighbor_links', 'num_undirected_edges', ...
        'cur', 'seed', 'train_ratio', 'smse_curve', 'rmse_curve');
end
end