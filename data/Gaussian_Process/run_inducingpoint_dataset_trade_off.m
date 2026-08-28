function run_inducingpoint_dataset_trade_off(DatasetName, CurrentMode, train_ratio, seed, NumInducingPoints_override)

if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1; end
if nargin < 5, NumInducingPoints_override = []; end
rng(seed);

% true  : compute predictive variance and NLPD/MSLL
% false : skip predictive variance; NLPD/MSLL = NaN; mean-only test time
compute_variance = true;

fprintf('\n[诱导点] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);
if compute_variance
    fprintf('IP test mode: full prediction with variance/NLPD\n');
else
    fprintf('IP test mode: mean-only prediction, NLPD disabled\n');
end

%% 1. 加载数据集
switch upper(DatasetName)
    case 'KIN40K'
        tr = load('KIN40K_train.mat');
        te = load('KIN40K_test.mat');
        hp = load('KIN40K_Hyperparameter.mat');
        train_x = tr.x; train_y = tr.y;
        test_x  = te.xtest; test_y = te.ytest;

    case 'POL'
        tr = load(fullfile('POL','POL_train.mat'));
        te = load(fullfile('POL','POL_test.mat'));
        hp = load(fullfile('POL','POL_Hyperparameter.mat'));
        train_x = tr.x; train_y = tr.y;
        test_x  = te.xtest; test_y = te.ytest;

    case 'PUMADYN32NM'
        tr = load(fullfile('PUMADYN32NM','PUMADYN32NM_train.mat'));
        te = load(fullfile('PUMADYN32NM','PUMADYN32NM_test.mat'));
        hp = load(fullfile('PUMADYN32NM','PUMADYN32NM_Hyperparameter.mat'));
        train_x = tr.x; train_y = tr.y;
        test_x  = te.xtest; test_y = te.ytest;

    case 'SARCOS'
        tr = load(fullfile('SARCOS','SARCOS_train.mat'));
        te = load(fullfile('SARCOS','SARCOS_test.mat'));
        hp_raw = load(fullfile('SARCOS','SARCOS_GP_Hyperparameter.mat'));
        hp.SigmaF = mean(cell2mat(hp_raw.SigmaF_set));
        hp.SigmaN = mean(cell2mat(hp_raw.SigmaN_set));
        hp.SigmaL = mean(cell2mat(hp_raw.SigmaL_set'), 2);
        train_x = tr.sarcos_inv(:,1:21);
        train_y = tr.sarcos_inv(:,22:28);
        test_x  = te.sarcos_inv_test(:,1:21);
        test_y  = te.sarcos_inv_test(:,22:28);

    otherwise
        error('未知数据集: %s', DatasetName);
end

%% 2. 归一化
if size(hp.SigmaL,1) > 1 && size(hp.SigmaL,2) > 1
    hp.SigmaL = mean(hp.SigmaL, 1);
end
SigmaL = hp.SigmaL(:);
if numel(hp.SigmaF) > 1, hp.SigmaF = mean(hp.SigmaF); end
if numel(hp.SigmaN) > 1, hp.SigmaN = mean(hp.SigmaN); end
SigmaF = hp.SigmaF;
SigmaN = hp.SigmaN;

n_train = round(size(train_x,1) * train_ratio);
idx_tr = randperm(size(train_x,1), n_train);
X_train = train_x(idx_tr,:);
Y_train = train_y(idx_tr,:);

idx_te = randperm(size(test_x,1));
X_test = test_x(idx_te,:);
Y_test = test_y(idx_te,:);

X_mean = mean(X_train,1);
X_std  = std(X_train,0,1);
X_std(X_std == 0) = 1;
if ~(max(abs(X_mean)) < 1e-2 && max(abs(X_std-1)) < 1e-2)
    X_train = (X_train - X_mean) ./ X_std;
    X_test  = (X_test  - X_mean) ./ X_std;
    SigmaL  = SigmaL ./ X_std(:);
end

Y_mean = mean(Y_train,1);
Y_std  = std(Y_train,0,1);
Y_std(Y_std == 0) = 1;
if max(abs(Y_mean)) < 1e-2 && max(abs(Y_std-1)) < 1e-2
    Y_mean = zeros(1,size(Y_train,2));
    Y_std  = ones(1,size(Y_train,2));
else
    Y_train = (Y_train - Y_mean) ./ Y_std;
    SigmaF  = SigmaF / mean(Y_std);
    SigmaN  = SigmaN / mean(Y_std);
end
prior_var = SigmaF^2;

[N_train, x_dim] = size(X_train);
y_dim = size(Y_train,2);
N_eval = min(3000, size(X_test,1));
X_eval = X_test(1:N_eval,:);
Y_eval = Y_test(1:N_eval,:);
Y_var_base = var(Y_eval,0,1);
fprintf('Train=%d Test=%d x=%d y=%d\n', N_train, N_eval, x_dim, y_dim);

%% 3. 分布式参数
AgentQuantity = 6;
Kappa_P = 10;
t_step = 0.01;
MaxDataPerAgent = min(floor(N_train/AgentQuantity), 3000);

if ~isempty(NumInducingPoints_override)
    NumInducingPoints = NumInducingPoints_override;
else
    switch upper(DatasetName)
        case {'KIN40K'}
            NumInducingPoints = 2500;
        case {'POL'}
            NumInducingPoints = 2000;
        case {'PUMADYN32NM'}
            NumInducingPoints = 500;
        otherwise
            NumInducingPoints = 2500;
    end
end
fprintf('NumInducingPoints M = %d\n', NumInducingPoints);

MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity,1);
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
N_degree = sum(L < 0, 2);
N_max = max(N_degree);

a_param = 0.5 / N_max;
sigma_i = 0.5;
sigma_i_ac = 0.5;
fprintf('ET参数: a=%.4g  sigma_DAC=%.4g  sigma_AC=%.4g  max|N_i|=%d\n', ...
    a_param, sigma_i, sigma_i_ac, N_max);

%% 4. 训练局部GP
tic_local_gp = tic;
LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    idx = (n-1)*MaxDataPerAgent+1 : min(n*MaxDataPerAgent, N_train);
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim, y_dim, MaxDataPerAgent, SigmaN, SigmaF, SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)', Y_train(idx,:)');
    LocalGP_set{n}.tau = 1e-8;
    LocalGP_set{n}.delta = 0.01;
end
t_ip_local_gp_train = toc(tic_local_gp);
t_train_gp = t_ip_local_gp_train; 
fprintf('局部GP: %.4fs\n', t_ip_local_gp_train);

idx_ind = randperm(N_train, NumInducingPoints);
InducingPoints_Coordinates = X_train(idx_ind,:)';

%% 5. 方法列表
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all')
    AllModes = [dac_methods, ac_methods];
else
    AllModes = {lower(CurrentMode)};
end

%% 6. 批量预计算P矩阵
p_dim = 2*y_dim;
P_poe  = zeros(p_dim,AgentQuantity,NumInducingPoints);
P_gpoe = zeros(p_dim,AgentQuantity,NumInducingPoints);
P_moe  = zeros(p_dim,AgentQuantity,NumInducingPoints);
P_bcm  = zeros(p_dim,AgentQuantity,NumInducingPoints);
P_rbcm = zeros(p_dim,AgentQuantity,NumInducingPoints);
mu_ind  = zeros(AgentQuantity,y_dim,NumInducingPoints);
var_ind = zeros(AgentQuantity,y_dim,NumInducingPoints);

fprintf('[预计算] %d个诱导点...\n', NumInducingPoints);
tic_inducing_prediction = tic;
for n = 1:AgentQuantity
    [mu_n_all, var_n_all] = batch_predict_external(LocalGP_set{n}, InducingPoints_Coordinates, SigmaN, SigmaF);
    for m = 1:NumInducingPoints
        mu_n = mu_n_all(m,:)';
        var_n = repmat(var_n_all(m), y_dim, 1);
        mu_ind(n,:,m) = mu_n';
        var_ind(n,:,m) = var_n';
        for d = 1:y_dim
            vs = max(var_n(d), SigmaN^2);
            b = max(min(0.5*(log(prior_var)-log(vs)),10),eps);

            P_poe(2*d-1,n,m) = AgentQuantity*mu_n(d)/vs;
            P_poe(2*d,n,m)   = AgentQuantity/vs;

            P_gpoe(2*d-1,n,m) = AgentQuantity*b*mu_n(d)/vs;
            P_gpoe(2*d,n,m)   = AgentQuantity*b/vs;

            P_moe(2*d-1,n,m) = AgentQuantity*mu_n(d);
            P_moe(2*d,n,m)   = AgentQuantity*(vs + mu_n(d)^2);

            P_bcm(2*d-1,n,m) = AgentQuantity*mu_n(d)/vs;
            P_bcm(2*d,n,m)   = AgentQuantity/vs - (AgentQuantity-1)/prior_var;

            P_rbcm(2*d-1,n,m) = AgentQuantity*b*mu_n(d)/vs;
            P_rbcm(2*d,n,m)   = AgentQuantity*b/vs + (1-AgentQuantity*b)/prior_var;
        end
    end
end
t_ip_inducing_prediction = toc(tic_inducing_prediction);
fprintf('预计算完成: %.2fs\n', t_ip_inducing_prediction);


ProjectRoot = fileparts(mfilename('fullpath'));
SaveFolder = fullfile(ProjectRoot, 'Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
tr_tag = round(train_ratio*100);

%% 7. 主循环
for mi = 1:numel(AllModes)
    cur = AllModes{mi};
    fprintf('\n[%d/%d] %s\n', mi, numel(AllModes), cur);
    trigger_count_per_agent = zeros(AgentQuantity,NumInducingPoints);
    broadcast_count_per_agent = zeros(AgentQuantity,1);

    base_method = strrep(lower(cur),'_ac','');
    if strcmpi(base_method,'gpoe')
        Pi = P_gpoe;
    elseif strcmpi(base_method,'poe')
        Pi = P_poe;
    elseif strcmpi(base_method,'bcm')
        Pi = P_bcm;
    elseif strcmpi(base_method,'rbcm')
        Pi = P_rbcm;
    else
        Pi = P_moe;
    end

    iter = 0;
    max_iter = 3000;
    ConsensusBaselineIter = 500;

    % Convergence diagnostics
    conv_r = 1;
    conv_m = 1;

    conv_curve_dac = [];
    conv_curve_ac  = [];
    conv_dac_tracking_curve     = [];
    conv_dac_disagreement_curve = [];
    conv_ac_avg_error_curve     = [];
    conv_ac_disagreement_curve  = [];
    dac_state_hist = [];
    dac_ref_hist   = [];
    ac_state_hist  = [];
    ac_ref_hist    = [];

    if ismember(lower(cur), dac_methods)
        %% IP-DAC
        tic_consensus = tic;

        Zeta = zeros(p_dim,AgentQuantity,NumInducingPoints);
        Xi_last_trigger = Pi;
        Pi_ref = mean(Pi, 2);

        conv_dac_tracking_curve     = zeros(max_iter,1);
        conv_dac_disagreement_curve = zeros(max_iter,1);
        dac_state_hist = zeros(max_iter,AgentQuantity);
        dac_ref_hist   = zeros(max_iter,1);
        conv_curve_dac = zeros(max_iter,1);

        while iter < max_iter
            iter = iter + 1;
            Zeta_prev = Zeta;

            L_Xi_hat = laplacian_multiply_agent_dim(Xi_last_trigger, L);
            Zeta = Zeta + t_step*Kappa_P*L_Xi_hat;
            Xi_now = Pi - Zeta;

            dac_disagreement = abs(Xi_now - mean(Xi_now,2));
            conv_dac_tracking_curve(iter)     = mean(max(abs(Xi_now-Pi_ref),[],2),'all');
            conv_dac_disagreement_curve(iter) = mean(max(dac_disagreement,[],2),'all');
            dac_state_hist(iter,:) = squeeze(Xi_now(conv_r,:,conv_m));
            dac_ref_hist(iter)     = mean(squeeze(Pi(conv_r,:,conv_m)));
            conv_curve_dac(iter)   = conv_dac_disagreement_curve(iter);

            for agent_i = 1:AgentQuantity
                neighbor_i = (L(agent_i,:) < 0);
                N_i = N_degree(agent_i);
                coeff_i = sigma_i * a_param * (1 - a_param*N_i) / N_i;

                E_i = Xi_last_trigger(:,agent_i,:) - Xi_now(:,agent_i,:);
                e_norm_sq = squeeze(sum(E_i.^2, 1));

                Z_i = N_i*Xi_now(:,agent_i,:) - sum(Xi_now(:,neighbor_i,:), 2);
                z_norm_sq = squeeze(sum(Z_i.^2, 1));

                trigger_idx = (e_norm_sq(:).' > coeff_i * z_norm_sq(:).');
                if iter == 1
                    trigger_idx(:) = true;
                end

                if any(trigger_idx)
                    Xi_last_trigger(:,agent_i,trigger_idx) = Xi_now(:,agent_i,trigger_idx);
                    trigger_count_per_agent(agent_i,trigger_idx) = ...
                        trigger_count_per_agent(agent_i,trigger_idx) + 1;
                    broadcast_count_per_agent(agent_i) = broadcast_count_per_agent(agent_i) + 1;
                end
            end

            if max(abs(Zeta(:)-Zeta_prev(:))) < 1e-5
                break;
            end
        end

        t_ip_consensus = toc(tic_consensus);
        iter_converge = iter;

        conv_dac_tracking_curve     = conv_dac_tracking_curve(1:iter_converge);
        conv_dac_disagreement_curve = conv_dac_disagreement_curve(1:iter_converge);
        dac_state_hist = dac_state_hist(1:iter_converge,:);
        dac_ref_hist   = dac_ref_hist(1:iter_converge);
        conv_curve_dac = conv_curve_dac(1:iter_converge);

        trigger_per_agent_point = mean(trigger_count_per_agent, 2);
        point_trigger_count_mean = mean(trigger_count_per_agent(:));
        broadcast_count_mean = mean(broadcast_count_per_agent);
        comm_train = point_trigger_count_mean;
        comm_test  = 0;
        trigger_ratio_train = comm_train / ConsensusBaselineIter;
        fprintf('  [IP-DAC] 收敛步数:%d  平均触发次数/agent/point:%.1f\n', ...
            iter_converge, comm_train);
        fprintf('  [IP-DAC broadcast stats] broadcast/agent:%.1f  point-trigger/agent/point:%.1f\n', ...
            broadcast_count_mean, point_trigger_count_mean);

        Xi_final = Pi - Zeta;

    else
        %% IP-AC
        tic_consensus = tic;

        Xi = Pi;
        Xi_last_trigger = Pi;
        Xi_initial_average = mean(Pi, 2);

        conv_ac_avg_error_curve    = zeros(max_iter,1);
        conv_ac_disagreement_curve = zeros(max_iter,1);
        ac_state_hist = zeros(max_iter,AgentQuantity);
        ac_ref_hist   = zeros(max_iter,1);
        conv_curve_ac = zeros(max_iter,1);

        while iter < max_iter
            iter = iter + 1;
            Xi_prev = Xi;

            L_Xi_hat = laplacian_multiply_agent_dim(Xi_last_trigger, L);
            Xi = Xi - t_step*Kappa_P*L_Xi_hat;

            ac_avg_err      = abs(Xi - Xi_initial_average);
            ac_disagreement = abs(Xi - mean(Xi,2));
            conv_ac_avg_error_curve(iter)    = mean(max(ac_avg_err, [], 2), 'all');
            conv_ac_disagreement_curve(iter) = mean(max(ac_disagreement, [], 2), 'all');

            ac_state_hist(iter,:) = squeeze(Xi(conv_r,:,conv_m));
            ac_ref_hist(iter)     = mean(squeeze(Pi(conv_r,:,conv_m)));
            conv_curve_ac(iter)   = conv_ac_avg_error_curve(iter);

            for agent_i = 1:AgentQuantity
                neighbor_i = (L(agent_i,:) < 0);
                N_i = N_degree(agent_i);
                coeff_i = sigma_i_ac * a_param * (1 - a_param*N_i) / N_i;

                E_i = Xi_last_trigger(:,agent_i,:) - Xi(:,agent_i,:);
                e_norm_sq = squeeze(sum(E_i.^2, 1));

                Z_i = N_i*Xi(:,agent_i,:) - sum(Xi(:,neighbor_i,:), 2);
                z_norm_sq = squeeze(sum(Z_i.^2, 1));

                trigger_idx = (e_norm_sq(:).' > coeff_i * z_norm_sq(:).');
                if iter == 1
                    trigger_idx(:) = true;
                end

                if any(trigger_idx)
                    Xi_last_trigger(:,agent_i,trigger_idx) = Xi(:,agent_i,trigger_idx);
                    trigger_count_per_agent(agent_i,trigger_idx) = ...
                        trigger_count_per_agent(agent_i,trigger_idx) + 1;
                    broadcast_count_per_agent(agent_i) = broadcast_count_per_agent(agent_i) + 1;
                end
            end

            if max(abs(Xi(:)-Xi_prev(:))) < 1e-5
                break;
            end
        end

        t_ip_consensus = toc(tic_consensus);
        iter_converge = iter;

        conv_ac_avg_error_curve    = conv_ac_avg_error_curve(1:iter_converge);
        conv_ac_disagreement_curve = conv_ac_disagreement_curve(1:iter_converge);
        ac_state_hist = ac_state_hist(1:iter_converge,:);
        ac_ref_hist   = ac_ref_hist(1:iter_converge);
        conv_curve_ac = conv_curve_ac(1:iter_converge);
        conv_curve_dac = [];

        trigger_per_agent_point = mean(trigger_count_per_agent, 2);
        point_trigger_count_mean = mean(trigger_count_per_agent(:));
        broadcast_count_mean = mean(broadcast_count_per_agent);
        comm_train = point_trigger_count_mean;
        comm_test  = 0;
        trigger_ratio_train = comm_train / ConsensusBaselineIter;
        fprintf('  [IP-AC] 收敛步数:%d  平均更新次数/agent/point:%.1f\n', ...
            iter_converge, comm_train);
        fprintf('  [IP-AC broadcast stats] broadcast/agent:%.1f  point-trigger/agent/point:%.1f\n', ...
            broadcast_count_mean, point_trigger_count_mean);

        Xi_final = Xi;
    end

    %% 提取phi
    phi = zeros(y_dim,NumInducingPoints);
    for d = 1:y_dim
        xi1 = squeeze(Xi_final(2*d-1,1,:))';
        xi2 = squeeze(Xi_final(2*d,  1,:))';
        if ismember(base_method, {'gpoe','poe','bcm','rbcm'})
            den = xi2;
            small_den = abs(den) < 1e-4;
            den(small_den & den >= 0) = 1e-4;
            den(small_den & den <  0) = -1e-4;
            phi(d,:) = xi1 ./ den;
        else
            phi(d,:) = xi1 / AgentQuantity;
        end
    end

    %% 训练MaskedGP
    tic_maskedgp = tic;
    MaskedGP = LocalGP_MultiOutput(x_dim,y_dim,NumInducingPoints,1e-6,SigmaF,SigmaL);
    MaskedGP.add_Alldata(InducingPoints_Coordinates, phi);
    t_ip_maskedgp_train = toc(tic_maskedgp);

    t_ip_total_train = t_ip_local_gp_train + t_ip_inducing_prediction + t_ip_consensus + t_ip_maskedgp_train;
    t_train = t_ip_total_train;

    %% MaskedGP prediction on test points
    tic_test_total = tic;

    Num_Inducing = MaskedGP.DataQuantity;
    Alpha_Vec    = MaskedGP.alpha(1:Num_Inducing,:);
    Cholesky_L   = MaskedGP.L(1:Num_Inducing,1:Num_Inducing);

    tic_kernel = tic;
    K_star = MaskedGP.kernel(MaskedGP.X(:,1:Num_Inducing), X_eval');
    t_ip_test_kernel = toc(tic_kernel);

    tic_mean = tic;
    mu_normalized = (Alpha_Vec' * K_star)';
    %mu_pred = mu_normalized .* repmat(Y_std, N_eval, 1) + repmat(Y_mean, N_eval, 1);
    mu_pred = mu_normalized .* Y_std + Y_mean;
    t_ip_test_mean = toc(tic_mean);

    if compute_variance
    tic_variance = tic;

    V_matrix = Cholesky_L \ K_star;
    var_normalized = max(SigmaF^2 - sum(V_matrix.^2, 1)', SigmaN^2);

    % MATLAB implicit expansion:
    % var_normalized: [N_eval x 1]
    % Y_std.^2      : [1 x y_dim]
    var_pred = var_normalized .* (Y_std.^2);

    t_ip_test_variance = toc(tic_variance);
    else
        var_pred = NaN(N_eval, y_dim);
        t_ip_test_variance = 0;
    end

    t_test = toc(tic_test_total);

    %% Metrics
    err = Y_eval - mu_pred;
    Y_var_base_safe = max(Y_var_base, 1e-10);
    smse = mean(mean(err.^2) ./ Y_var_base_safe);
    rmse = mean(sqrt(mean(err.^2)));

    if compute_variance
        nlpd = mean(mean(0.5*(log(2*pi*var_pred) + err.^2 ./ var_pred)));

        % MSLL (Mean Standardized Log Loss), per Rasmussen & Williams,
        % GPML book, Eq. 2.34: "This loss can be standardized by
        % subtracting the loss that would be obtained under the trivial
        % model which predicts using a Gaussian with the MEAN AND
        % VARIANCE OF THE TRAINING DATA."
        % => baseline must use Y_mean / Y_std (training-set stats, already
        %    computed in original scale above), NOT the test set.
        % (Verified against the textbook screenshot; the earlier switch to
        % Y_eval-based baseline was incorrect and has been reverted.)
        Y_train_var_safe = max(Y_std.^2, 1e-10);   % [1 x y_dim], numerical safety
        err_trivial = Y_eval - Y_mean;             % implicit expansion: [N_eval x y_dim]
        nlpd_trivial = mean(mean(0.5*(log(2*pi*Y_train_var_safe) + ...
            err_trivial.^2 ./ Y_train_var_safe)));

        msll = nlpd - nlpd_trivial;
    else
        nlpd = NaN;
        msll = NaN;
    end

    t_train_per_point = (t_train/N_train) * 1000;
    t_test_per_point  = (t_test/N_eval) * 1000;

    if compute_variance
        fprintf('  SMSE=%.4f RMSE=%.4f NLPD=%.4f MSLL=%.4f Train:%.2fms/pt Test:%.2fms/pt\n', ...
            smse, rmse, nlpd, msll, t_train_per_point, t_test_per_point);
    else
        fprintf('  SMSE=%.4f RMSE=%.4f NLPD=NaN MSLL=NaN(mean-only) Train:%.2fms/pt Test:%.2fms/pt\n', ...
            smse, rmse, t_train_per_point, t_test_per_point);
    end

    fprintf('  Timing IP: LocalGP=%.3fs  IndPred=%.3fs  Consensus=%.3fs  MaskedGP=%.3fs  Total=%.3fs\n', ...
        t_ip_local_gp_train, t_ip_inducing_prediction, t_ip_consensus, t_ip_maskedgp_train, t_ip_total_train);
    fprintf('  IP test timing: Kernel=%.3fs  Mean=%.3fs  Variance=%.3fs  Total=%.3fs\n', ...
        t_ip_test_kernel, t_ip_test_mean, t_ip_test_variance, t_test);

    %% Compatibility variables expected by run_mc_dataset / plotting scripts
    t_train_total = t_train;
    t_test_total  = t_test;
    current_method = cur;
    t_total_train = t_ip_total_train;
    t_consensus = t_ip_consensus;

    err_sq_mean = mean(err.^2,2);
    smse_curve = cumsum(err_sq_mean) ./ (1:N_eval)' / mean(Y_var_base_safe);
    rmse_curve = sqrt(cumsum(err_sq_mean) ./ (1:N_eval)');
    event_count_mean = comm_train;

    save(fullfile(SaveFolder, sprintf('%s_M%d_tr%d_mc%d.mat', current_method, NumInducingPoints, tr_tag, seed)), ...
        'smse', 'rmse', 'nlpd', 'msll', ...
        't_train_total', 't_test_total', ...
        't_train_per_point', 't_test_per_point', ...
        't_ip_local_gp_train', 't_ip_inducing_prediction', ...
        't_ip_consensus', 't_ip_maskedgp_train', 't_ip_total_train', ...
        't_total_train', 't_consensus', ...
        'compute_variance', ...
        't_ip_test_kernel', 't_ip_test_mean', 't_ip_test_variance', ...
        'comm_train', 'comm_test', 'iter_converge', ...
        'event_count_mean', 'trigger_ratio_train', 'ConsensusBaselineIter', ...
        'current_method', 'seed', 'train_ratio', 'NumInducingPoints', 'N_train', 'N_eval', ...
        'smse_curve', 'rmse_curve', ...
        'trigger_count_per_agent', 'trigger_per_agent_point', ...
        'point_trigger_count_mean', 'broadcast_count_per_agent', 'broadcast_count_mean', ...
        'conv_curve_dac', 'conv_curve_ac', ...
        'conv_dac_tracking_curve', 'conv_dac_disagreement_curve', ...
        'conv_ac_avg_error_curve', 'conv_ac_disagreement_curve', ...
        'dac_state_hist', 'dac_ref_hist', ...
        'ac_state_hist', 'ac_ref_hist');
end
end

function L_X = laplacian_multiply_agent_dim(X, L)
%LAPLACIAN_MULTIPLY_AGENT_DIM Apply graph Laplacian along agent dimension.
%   X has size [p_dim x AgentQuantity x NumPoints].
%   L has size [AgentQuantity x AgentQuantity].
%   The result has the same size as X and satisfies
%       L_X(:,i,m) = sum_j L(i,j) * X(:,j,m).

    [p_dim, agent_quantity, num_points] = size(X);
    X_agent_first = permute(X, [2, 1, 3]);
    X_flat = reshape(X_agent_first, agent_quantity, []);
    L_X_flat = L * X_flat;
    L_X_agent_first = reshape(L_X_flat, agent_quantity, p_dim, num_points);
    L_X = permute(L_X_agent_first, [2, 1, 3]);
end
