function run_dac_dataset(DatasetName, CurrentMode)
%  DatasetName: 'KIN40K'
%  CurrentMode: 'moe' | 'gpoe' | 'moe_ac' | 'gpoe_ac' | 'all'
rng(0);

%% 1. 加载数据集
fprintf('加载数据集: %s\n', DatasetName);
switch upper(DatasetName)
    case 'KIN40K'
        train = load('KIN40K_train.mat');
        test  = load('KIN40K_test.mat');
        hp    = load('KIN40K_Hyperparameter.mat');
        X_train = train.x;      Y_train = train.y;
        X_test  = test.xtest;   Y_test  = test.ytest;
    case 'POL'
        train = load(fullfile('POL', 'POL_train.mat'));
        test  = load(fullfile('POL', 'POL_test.mat'));
        hp    = load(fullfile('POL', 'POL_Hyperparameter.mat'));
        X_train = train.x;      Y_train = train.y;
        X_test  = test.xtest;   Y_test  = test.ytest;
    case 'PUMADYN32NM'
        train = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_train.mat'));
        test  = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_test.mat'));
        hp    = load(fullfile('PUMADYN32NM', 'PUMADYN32NM_Hyperparameter.mat'));
        X_train = train.x;      Y_train = train.y;
        X_test  = test.xtest;   Y_test  = test.ytest;
    case 'SARCOS'
        train = load(fullfile('SARCOS', 'SARCOS_train.mat'));
        test  = load(fullfile('SARCOS', 'SARCOS_test.mat'));
        hp    = load(fullfile('SARCOS', 'SARCOS_GP_Hyperparameter.mat'));
        % 先用第1组
        hp.SigmaF = hp.SigmaF_set{1};
        hp.SigmaL = hp.SigmaL_set{1};
        hp.SigmaN = hp.SigmaN_set{1};
        X_train = train.sarcos_inv(:, 1:21);
        Y_train = train.sarcos_inv(:, 22:28);
        X_test  = test.sarcos_inv_test(:, 1:21);
        Y_test  = test.sarcos_inv_test(:, 22:28);
    otherwise
        error('未知数据集: %s', DatasetName);
end

[N_train, x_dim] = size(X_train);
y_dim      = size(Y_train, 2);
N_eval     = min(30000, size(X_test,1));
X_eval     = X_test(1:N_eval, :);
Y_eval     = Y_test(1:N_eval, :);
Y_var_base = var(Y_eval);
fprintf('Train: %d x %d, y_dim: %d, Eval: %d\n', N_train, x_dim, y_dim, N_eval);

%% 2. 参数
AgentQuantity     = 6;
NumInducingPoints = 400;
Kappa_P           = 10;
t_step            = 0.01;
SigmaF            = hp.SigmaF;
SigmaL            = hp.SigmaL; 
SigmaN            = hp.SigmaN;
prior_var         = SigmaF^2;
MaxDataPerAgent = min(floor(N_train / AgentQuantity), 1500);

%% 3. Topology
AgentQuantity = 6; LeaderQuantity = 1;
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity, LeaderQuantity);
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

%% 4. train local GP
fprintf('\n训练 %d 个局部GP...\n', AgentQuantity);
tic;
LocalGP_set = cell(AgentQuantity, 1);
for n = 1:AgentQuantity
    idx = (n-1)*MaxDataPerAgent+1 : min(n*MaxDataPerAgent, N_train);
    X_n = X_train(idx,:)';   % x_dim x N
    Y_n = Y_train(idx,:)';   % y_dim x N
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim, y_dim, MaxDataPerAgent, SigmaN, SigmaF, SigmaL);
    LocalGP_set{n}.add_Alldata(X_n, Y_n);
    LocalGP_set{n}.tau   = 1e-8;
    LocalGP_set{n}.delta = 0.01;
end
fprintf('训练耗时: %.2f 秒\n', toc);

%% 5. inducing_points
X_min = min(X_train)';
X_max = max(X_train)';
InducingPoints_Coordinates = X_min + (X_max - X_min) .* rand(x_dim, NumInducingPoints);

%% 6. 
dac_methods = {'moe', 'gpoe'};
ac_methods  = {'moe_ac', 'gpoe_ac'};
if strcmpi(CurrentMode, 'all')
    AllModes = [dac_methods, ac_methods];
else
    AllModes = {lower(CurrentMode)};
end

%% 7. p_dim（MoE/gPoE: 2*y_dim，每个维度[num, den]）
p_dim = 2 * y_dim;

%% 8. 计算P
fprintf('\n计算诱导点P...\n');

P_gpoe = zeros(p_dim, AgentQuantity, NumInducingPoints);
P_moe  = zeros(p_dim, AgentQuantity, NumInducingPoints);
for InducingPointIdx = 1:NumInducingPoints
    x_m = InducingPoints_Coordinates(:, InducingPointIdx);
    for n = 1:AgentQuantity
        [mu_n, var_n] = LocalGP_set{n}.predict(x_m);
        for d = 1:y_dim

            beta_nd = max(0.5 * (log(prior_var) - log(var_n(d))), eps);
            P_gpoe(2*d-1, n, InducingPointIdx) = AgentQuantity * beta_nd * mu_n(d) / var_n(d);
            P_gpoe(2*d,   n, InducingPointIdx) = AgentQuantity * beta_nd / var_n(d);
            
            P_moe(2*d-1, n, InducingPointIdx) = AgentQuantity * mu_n(d);
            P_moe(2*d,   n, InducingPointIdx) = AgentQuantity * (var_n(d) + mu_n(d)^2);
        end
    end
end

%% 9. 保存路径
SaveFolder = fullfile('Result', 'Dataset', DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end

results = zeros(numel(AllModes), 5);  % [SMSE, RMSE, NLPD, Train_T, Test_T]

for mi = 1:numel(AllModes)
    cur_method = AllModes{mi};
    fprintf('\n[%d/%d] %s\n', mi, numel(AllModes), cur_method);

    tic;
    switch lower(cur_method)
        case dac_methods
            
            base_method = cur_method;
            if strcmpi(cur_method, 'gpoe')
                P_inducing = P_gpoe;
            else  % moe
                P_inducing = P_moe;
            end
            Zeta = zeros(p_dim, AgentQuantity, NumInducingPoints);
            for dac_step = 1:400
                for InducingPointIdx = 1:NumInducingPoints    
                    P_InducingPoint    = P_inducing(:, :, InducingPointIdx);
                    Zeta_InducingPoint = Zeta(:, :, InducingPointIdx);
                    New_Consensus_Zeta_function = @(~, Zeta_ODE_Intern) Compute_New_Consensus_Derivative( ...
                        Zeta_ODE_Intern, P_InducingPoint, L, Kappa_P, AgentQuantity, p_dim);
                    [~, Zeta_ODE_Output] = ode45(New_Consensus_Zeta_function, ...
                        [0, t_step], Zeta_InducingPoint(:));
                    Zeta(:, :, InducingPointIdx) = reshape(Zeta_ODE_Output(end,:)', p_dim, AgentQuantity);
                end
            end
            Xi = P_inducing - Zeta;
            
            phi = zeros(y_dim, NumInducingPoints);
            for d = 1:y_dim
                xi1 = squeeze(Xi(2*d-1, 1, :))';  % 1 x M
                xi2 = squeeze(Xi(2*d,   1, :))';
                if strcmpi(cur_method, 'gpoe')
                    % gPoE: phi = xi1 / xi2 = sum(beta*mu/var) / sum(beta/var)
                    phi(d,:) = xi1 ./ xi2;
                else
                    % MoE: phi = xi1 / K  (Eq 3.36)
                    phi(d,:) = xi1 / AgentQuantity;
                end
            end

        case ac_methods            
            base_method = strrep(cur_method, '_ac', '');
            phi = zeros(y_dim, NumInducingPoints);
            for InducingPointIdx = 1:NumInducingPoints
                x_m = InducingPoints_Coordinates(:, InducingPointIdx);
                mu_all  = zeros(AgentQuantity, y_dim);
                var_all = zeros(AgentQuantity, y_dim);
                for n = 1:AgentQuantity
                    [mu_n, var_n] = LocalGP_set{n}.predict(x_m);
                    mu_all(n,:) = mu_n';
                    var_all(n,:) = var_n';
                end
                for d = 1:y_dim
                    switch base_method
                        case 'moe'
                            % MoE AC: 直接取均值 (Eq 3.36)
                            phi(d,InducingPointIdx) = mean(mu_all(:,d));
                        case 'gpoe'
                            % gPoE AC: 精度加权均值 (Eq 3.22)
                            beta_all = max(0.5*(log(prior_var) - log(var_all(:,d))), eps);
                            prec = sum(beta_all ./ var_all(:,d));
                            phi(d,InducingPointIdx) = sum(beta_all .* mu_all(:,d) ./ var_all(:,d)) / prec;
                    end
                end
            end
    end

    % 重建 MaskedGP
    MaskedGP = LocalGP_MultiOutput(x_dim, y_dim, NumInducingPoints, 1e-6, SigmaF, SigmaL);
    MaskedGP.add_Alldata(InducingPoints_Coordinates, phi);
    t_train = toc;

    % 预测
    mu_pred  = zeros(N_eval, y_dim);
    var_pred = ones(N_eval,  y_dim);
    tic;
    for n = 1:N_eval
        [mu_n, var_n] = MaskedGP.predict(X_eval(n,:)');
        mu_pred(n,:)  = mu_n';
        var_pred(n,:) = max(var_n', 1e-10);
        if mod(n, 5000) == 0
            fprintf('  %d/%d\n', n, N_eval);
        end
    end
    t_test = toc;

    % 计算指标
    error    = Y_eval - mu_pred;
    smse     = mean(mean(error.^2) ./ Y_var_base);
    rmse     = mean(sqrt(mean(error.^2)));
    nlpd     = mean(mean(0.5*(log(2*pi*var_pred) + error.^2./var_pred)));
    %coverage = mean(mean(abs(error) <= 1.96*sqrt(var_pred)));

    results(mi,:) = [smse, rmse, nlpd, t_train, t_test];

    save(fullfile(SaveFolder,[cur_method,'.mat']), ...
        'smse','rmse','nlpd','t_train','t_test','cur_method');
end

fprintf('完成。\n');
end

%% compute zeta_dot
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_Ref, L, Kappa, AgentQuantity, p_dim)

Zeta = reshape(Zeta_vec, p_dim, AgentQuantity);
dZeta_dt = Kappa * (P_Ref - Zeta) * L';
dZeta_dt = dZeta_dt(:);
end