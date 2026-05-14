function run_inducingpoint_dataset(DatasetName, CurrentMode, train_ratio, seed)
%  诱导点聚合方法
%  DatasetName  : 'KIN40K' | 'POL' | 'PUMADYN32NM' | 'SARCOS'
%  CurrentMode  : 'moe'|'gpoe'|'poe'|'bcm'|'rbcm'|'all'
%  train_ratio  : 训练集比例，如 0.4
%  seed         : 随机种子

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
N_eval = min(30000,size(X_test,1));
X_eval = X_test(1:N_eval,:);  Y_eval = Y_test(1:N_eval,:);
Y_var_base = var(Y_eval,0,1);
fprintf('Train=%d  Test=%d  x=%d  y=%d\n', N_train,N_eval,x_dim,y_dim);

%% 3. 分布式参数
AgentQuantity = 6;  Kappa_P = 10;  t_step = 0.01;
MaxDataPerAgent = min(floor(N_train/AgentQuantity), 3000);
switch upper(DatasetName)
    case {'SARCOS','POL'},  NumInducingPoints = 1500;
    otherwise,              NumInducingPoints = 800;
end
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity,1);
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;

%% 4. 训练局部 GP
tic;
LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    idx = (n-1)*MaxDataPerAgent+1 : min(n*MaxDataPerAgent,N_train);
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim,y_dim,MaxDataPerAgent,SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)', Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8;  LocalGP_set{n}.delta=0.01;
end
fprintf('局部GP: %.2fs\n', toc);

idx_ind = randperm(N_train, NumInducingPoints);
InducingPoints_Coordinates = X_train(idx_ind,:)';

%% 5. 方法
dac_methods = {'poe','gpoe','moe','bcm','rbcm'};
ac_methods  = {'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all'), AllModes=[dac_methods,ac_methods];
else, AllModes={lower(CurrentMode)}; end

%% 6. 预计算 P 矩阵
p_dim = 2*y_dim;
P_poe=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_gpoe=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_moe=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_bcm=zeros(p_dim,AgentQuantity,NumInducingPoints);
P_rbcm=zeros(p_dim,AgentQuantity,NumInducingPoints);
mu_ind=zeros(AgentQuantity,y_dim,NumInducingPoints);
var_ind=zeros(AgentQuantity,y_dim,NumInducingPoints);

for m = 1:NumInducingPoints
    x_m = InducingPoints_Coordinates(:,m);
    for n = 1:AgentQuantity
        [mu_n,var_n] = LocalGP_set{n}.predict(x_m);
        mu_ind(n,:,m)=mu_n'; var_ind(n,:,m)=var_n';
        for d = 1:y_dim
            b = 0.5*(log(prior_var)-log(var_n(d)));
            P_poe(2*d-1,n,m)  = AgentQuantity*mu_n(d)/var_n(d);
            P_poe(2*d,  n,m)  = AgentQuantity/var_n(d);

            P_gpoe(2*d-1,n,m) = AgentQuantity*b*mu_n(d)/var_n(d);
            P_gpoe(2*d,  n,m) = AgentQuantity*b/var_n(d);

            P_moe(2*d-1,n,m)  = AgentQuantity*mu_n(d);
            P_moe(2*d,  n,m)  = AgentQuantity*(var_n(d)+mu_n(d)^2);

            P_bcm(2*d-1,n,m)  = AgentQuantity*mu_n(d)/var_n(d);
            P_bcm(2*d,  n,m)  = AgentQuantity/var_n(d)-(AgentQuantity-1)/prior_var;

            P_rbcm(2*d-1,n,m) = AgentQuantity*b*mu_n(d)/var_n(d);
            P_rbcm(2*d,  n,m) = AgentQuantity*b/var_n(d)+(1-AgentQuantity*b)/prior_var;
        end
    end
end

SaveFolder = fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
tr_tag = round(train_ratio*100);

%% 7. 各方法主循环
for mi = 1:numel(AllModes)
    cur = AllModes{mi};
    fprintf('\n[%d/%d] %s\n', mi, numel(AllModes), cur);
    tic;

    switch lower(cur)
        case dac_methods
            if     strcmpi(cur,'gpoe'), Pi=P_gpoe;
            elseif strcmpi(cur,'poe'),  Pi=P_poe;
            elseif strcmpi(cur,'bcm'),  Pi=P_bcm;
            elseif strcmpi(cur,'rbcm'), Pi=P_rbcm;
            else,                       Pi=P_moe;
            end
            Zeta=zeros(p_dim,AgentQuantity,NumInducingPoints);
            for iter=1:3000
                prev=Zeta; diff=Pi-Zeta;
                for n=1:AgentQuantity
                    Zeta(:,n,:)=Zeta(:,n,:)+t_step*Kappa_P*sum(diff.*reshape(L(n,:),1,AgentQuantity,1),2);
                end
                if max(abs(Zeta(:)-prev(:)))<1e-5, fprintf('  收敛: %d步\n',iter); break; end
            end
            Xi=Pi-Zeta;
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
            base=strrep(lower(cur),'_ac','');
            phi=zeros(y_dim,NumInducingPoints);
            for m=1:NumInducingPoints
                for d=1:y_dim
                    mu_a=squeeze(mu_ind(:,d,m)); va=squeeze(var_ind(:,d,m));
                    b=max(0.5*(log(prior_var)-log(va)),eps);
                    switch base
                        case 'moe',  phi(d,m)=mean(mu_a);
                        case 'gpoe', phi(d,m)=sum(b.*mu_a./va)/max(sum(b./va),eps);
                        case 'poe',  phi(d,m)=sum(mu_a./va)/max(sum(1./va),eps);
                        case 'bcm'
                            prec=sum(1./va)-(AgentQuantity-1)/prior_var;
                            phi(d,m)=sum(mu_a./va)/max(prec,eps);
                        case 'rbcm'
                            prec=sum(b./va)+(1-sum(b))/prior_var;
                            phi(d,m)=sum(b.*mu_a./va)/max(prec,eps);
                    end
                end
            end
    end

    MaskedGP=LocalGP_MultiOutput(x_dim,y_dim,NumInducingPoints,1e-6,SigmaF,SigmaL);
    MaskedGP.add_Alldata(InducingPoints_Coordinates,phi);
    t_train=toc;

    tic;
    % 1. 提取 MaskedGP 已经训练好的内部参数
    Num_Inducing = MaskedGP.DataQuantity;
    Alpha_Vec    = MaskedGP.alpha(1:Num_Inducing, :);
    Cholesky_L   = MaskedGP.L(1:Num_Inducing, 1:Num_Inducing);

    % 2. 计算诱导点与所有测试点之间的核矩阵 K_star (Ks)
    K_star = MaskedGP.kernel(MaskedGP.X(:, 1:Num_Inducing), X_eval');

    % 3. 极速计算预测均值 (归一化空间)
    mu_normalized = (Alpha_Vec' * K_star)';

    % 4. 极速计算预测方差 (归一化空间)
    % 使用 Cholesky 因子做线性求解 V = L^(-1) * K_star，避免直接求逆
    V_matrix = Cholesky_L \ K_star;
    % prior_var (SigmaF^2) 减去信息增益，并用噪声方差 SigmaN^2 兜底防崩溃
    var_normalized = max(SigmaF^2 - sum(V_matrix.^2, 1)', SigmaN^2);

    % 5. 反归一化：将预测结果还原到真实的物理量纲
    mu_pred  = mu_normalized .* repmat(Y_std, N_eval, 1) + repmat(Y_mean, N_eval, 1);
    var_pred = repmat(var_normalized, 1, y_dim) .* repmat(Y_std.^2, N_eval, 1);

    t_test = toc;

    err=Y_eval-mu_pred;
    smse=mean(mean(err.^2)./Y_var_base);
    rmse=mean(sqrt(mean(err.^2)));
    nlpd=mean(mean(0.5*(log(2*pi*var_pred)+err.^2./var_pred)));
    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train=%.1fs  Test=%.1fs\n', ...
        smse,rmse,nlpd,t_train,t_test);

    % 累积曲线（向量化，无需循环）
    err_sq_mean  = mean(err.^2, 2);           % N_eval × 1
    smse_curve   = cumsum(err_sq_mean) ./ (1:N_eval)' / mean(Y_var_base);
    rmse_curve   = sqrt(cumsum(err_sq_mean) ./ (1:N_eval)');

    save(fullfile(SaveFolder,sprintf('%s_tr%d_mc%d.mat',cur,tr_tag,seed)), ...
        'smse','rmse','nlpd','t_train','t_test','cur','seed','train_ratio', ...
        'smse_curve','rmse_curve');
end
fprintf('\n[%s] 诱导点 done. seed=%d tr=%d%%\n\n', DatasetName, seed, tr_tag);
end