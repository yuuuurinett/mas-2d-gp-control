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

%% 2. 归一化
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

%% 3. 分布式参数与拓扑
AgentQuantity=6; Kappa_P=10; t_step=0.01;
MaxDataPerAgent=min(floor(N_train/AgentQuantity),3000);
MultiAgentSystem=Manipulator_2D_2DoF_SetMASTopology(AgentQuantity,1);
L=MultiAgentSystem.Agent_Topology.LaplacianMatrix;

N_degree=sum(L<0,2); N_max=max(N_degree);
sigma_i=0.1; a_param=0.9/N_max;
fprintf('ET 参数: sigma=%.2f  a=%.4f  (N_max=%d)\n',sigma_i,a_param,N_max);

%% 4. 训练局部 GP
t_start=tic;
LocalGP_set=cell(AgentQuantity,1);
for n=1:AgentQuantity
    idx=(n-1)*MaxDataPerAgent+1:min(n*MaxDataPerAgent,N_train);
    LocalGP_set{n}=LocalGP_MultiOutput(x_dim,y_dim,MaxDataPerAgent,SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)',Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8; LocalGP_set{n}.delta=0.01;
end
t_train=toc(t_start);
fprintf('局部GP训练完成: %.4fs\n',t_train);

%% 5. 方法列表
dac_methods={'moe','gpoe','poe','bcm','rbcm'};
ac_methods={'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all'), AllModes=[dac_methods,ac_methods];
else, AllModes={lower(CurrentMode)}; end

SaveFolder=fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
p_dim=2*y_dim; tr_tag=round(train_ratio*100);

%% 6. 批量预计算
fprintf('\n[批量预计算] 正在预计算 %d 个测试点...\n', N_eval);
tic;
Mu_Local_All  = zeros(AgentQuantity, y_dim, N_eval);
Var_Local_All = zeros(AgentQuantity, y_dim, N_eval);
X_eval_matrix = X_eval';

for n = 1:AgentQuantity
    [mn_batch, vn_batch] = batch_predict_external( ...
        LocalGP_set{n}, X_eval_matrix, SigmaN, SigmaF);
    for t = 1:N_eval
        Mu_Local_All(n,:,t)  = mn_batch(t,:);
        Var_Local_All(n,:,t) = repmat(vn_batch(t), 1, y_dim);
    end
end
Precompute_Time=toc;
fprintf('批量预计算完成！耗时: %.2fs (%.4f ms/pt)\n',Precompute_Time,Precompute_Time/N_eval*1000);

%% 7. 聚合方法主循环
for mi=1:numel(AllModes)
    Current_Method=AllModes{mi};
    Method_Base_Name=strrep(lower(Current_Method),'_ac','');
    fprintf('\n[%d/%d] 正在执行聚合: %s\n',mi,numel(AllModes),Current_Method);

    Final_Mean_Pred=zeros(N_eval,y_dim);
    Final_Var_Pred =zeros(N_eval,y_dim);
    tic;

    if ismember(lower(Current_Method),dac_methods)
        %% --- A. DAC ---
        P_Info_Matrix=zeros(p_dim,AgentQuantity,N_eval);
        for t=1:N_eval
            for n=1:AgentQuantity
                Local_Mu =Mu_Local_All(n,:,t)';
                Local_Var=Var_Local_All(n,:,t)';
                for d=1:y_dim
                    Local_Var_safe = max(Local_Var(d), SigmaN^2);
                    Beta = max(min(0.5*(log(prior_var)-log(Local_Var_safe)),10),eps);
                    switch lower(Current_Method)
                        case 'moe'
                            P_Info_Matrix(2*d-1,n,t)=AgentQuantity*Local_Mu(d);
                            P_Info_Matrix(2*d,  n,t)=AgentQuantity*(Local_Var_safe+Local_Mu(d)^2);
                        case 'gpoe'
                            P_Info_Matrix(2*d-1,n,t)=AgentQuantity*Beta*Local_Mu(d)/Local_Var_safe;
                            P_Info_Matrix(2*d,  n,t)=AgentQuantity*Beta/Local_Var_safe;
                        case 'poe'
                            P_Info_Matrix(2*d-1,n,t)=AgentQuantity*Local_Mu(d)/Local_Var_safe;
                            P_Info_Matrix(2*d,  n,t)=AgentQuantity/Local_Var_safe;
                        case 'bcm'
                            P_Info_Matrix(2*d-1,n,t)=AgentQuantity*Local_Mu(d)/Local_Var_safe;
                            P_Info_Matrix(2*d,  n,t)=AgentQuantity/Local_Var_safe-(AgentQuantity-1)/prior_var;
                        case 'rbcm'
                            P_Info_Matrix(2*d-1,n,t)=AgentQuantity*Beta*Local_Mu(d)/Local_Var_safe;
                            P_Info_Matrix(2*d,  n,t)=AgentQuantity*Beta/Local_Var_safe+(1-AgentQuantity*Beta)/prior_var;
                    end
                end
            end
        end

        Zeta  =zeros(p_dim,AgentQuantity,N_eval);
        Zeta_k=zeros(p_dim,AgentQuantity,N_eval);
        trigger_count_set=zeros(AgentQuantity,1);
        iter = 0;
        max_iter = 3000;  % 安全上限，防止死循环
        comm_train=0;

        while iter < max_iter
            iter = iter + 1;
            Prev_Zeta=Zeta;

            % 动力学更新（Pi - Zeta），保证收敛性
            diff=P_Info_Matrix-Zeta;
            for n=1:AgentQuantity
                Zeta(:,n,:)=Zeta(:,n,:)+t_step*Kappa_P*...
                    sum(diff.*reshape(L(n,:),1,AgentQuantity,1),2);
            end

            % ET 触发判定：统计实际通信次数
            for n=1:AgentQuantity
                e_i=Zeta_k(:,n,:)-Zeta(:,n,:);
                z_i=sum((Zeta(:,n,:)-Zeta).*reshape(L(n,:)<0,1,AgentQuantity,1),2);
                rho_i    =max(sum(e_i.^2,1));
                norm_z_sq=max(sum(z_i.^2,1));
                N_i=N_degree(n);
                rho_bar_i=(sigma_i*a_param*(1-a_param*N_i)/N_i)*norm_z_sq;
                if rho_i>rho_bar_i||iter==1
                    Zeta_k(:,n,:)=Zeta(:,n,:);
                    trigger_count_set(n)=trigger_count_set(n)+1;
                end
            end

            % 收敛判定
            if max(abs(Zeta(:)-Prev_Zeta(:)))<1e-5
                break;
            end
        end

        iter_converge=iter;  % while 结束后直接就是收敛步数
        comm_test=mean(trigger_count_set);
        fprintf('    -> 收敛步数: %d  平均通信: %.1f次\n',iter_converge,comm_test);

        Final_Xi=P_Info_Matrix-Zeta;
        for t=1:N_eval
            for d=1:y_dim
                Xi_1=Final_Xi(2*d-1,1,t); Xi_2=Final_Xi(2*d,1,t);
                if ismember(lower(Current_Method),{'gpoe','poe','bcm','rbcm'})
                    Final_Mean_Pred(t,d)=Xi_1/max(Xi_2,eps);
                else
                    Final_Mean_Pred(t,d)=Xi_1/AgentQuantity;
                end
            end
        end

    else
        %% --- B. AC ---
        comm_train=0; comm_test=1; iter_converge=1;
        for t=1:N_eval
            for d=1:y_dim
                All_Mu =Mu_Local_All(:,d,t);
                All_Var=Var_Local_All(:,d,t);
                Beta_Set=max(0.5*(log(prior_var)-log(All_Var)),eps);
                switch Method_Base_Name
                    case 'moe',  Final_Mean_Pred(t,d)=mean(All_Mu);
                    case 'gpoe', Final_Mean_Pred(t,d)=sum(Beta_Set.*All_Mu./All_Var)/max(sum(Beta_Set./All_Var),eps);
                    case 'poe',  Final_Mean_Pred(t,d)=sum(All_Mu./All_Var)/max(sum(1./All_Var),eps);
                    case 'bcm'
                        Prec=sum(1./All_Var)-(AgentQuantity-1)/prior_var;
                        Final_Mean_Pred(t,d)=sum(All_Mu./All_Var)/max(Prec,eps);
                    case 'rbcm'
                        Prec=sum(Beta_Set./All_Var)+(1-sum(Beta_Set))/prior_var;
                        Final_Mean_Pred(t,d)=sum(Beta_Set.*All_Mu./All_Var)/max(Prec,eps);
                end
            end
        end
    end

    %% --- C. 方差聚合 ---
    for t=1:N_eval
        for d=1:y_dim
            All_Var=Var_Local_All(:,d,t);
            Beta_Set=max(0.5*(log(prior_var)-log(All_Var)),eps);
            switch Method_Base_Name
                case 'moe',  Final_Var_Pred(t,d)=mean(All_Var);
                case 'gpoe', Final_Var_Pred(t,d)=1/max(sum(Beta_Set./All_Var),eps);
                case 'poe',  Final_Var_Pred(t,d)=1/max(sum(1./All_Var),eps);
                case 'bcm'
                    Prec=sum(1./All_Var)-(AgentQuantity-1)/prior_var;
                    Final_Var_Pred(t,d)=1/max(Prec,eps);
                case 'rbcm'
                    Prec=sum(Beta_Set./All_Var)+(1-sum(Beta_Set))/prior_var;
                    Final_Var_Pred(t,d)=1/max(Prec,eps);
            end
            Final_Var_Pred(t,d)=max(Final_Var_Pred(t,d),SigmaN^2);
        end
    end

    t_test=Precompute_Time/numel(AllModes)+toc;

    %% --- D. 反归一化与保存 ---
    mu_pred =Final_Mean_Pred.*repmat(Y_std,N_eval,1)+repmat(Y_mean,N_eval,1);
    var_pred=Final_Var_Pred .*repmat(Y_std.^2,N_eval,1);
    err=Y_eval-mu_pred;
    smse=mean(mean(err.^2)./Y_var_base);
    rmse=mean(sqrt(mean(err.^2)));
    nlpd=mean(mean(0.5*(log(2*pi*var_pred)+err.^2./var_pred)));
    t_train_per_point=(t_train/N_train)*1000;
    t_test_per_point =(t_test/N_eval)*1000;
    fprintf('  SMSE=%.4f  RMSE=%.4f  NLPD=%.4f  Train:%.2fms/pt  Test:%.2fms/pt\n',...
        smse,rmse,nlpd,t_train_per_point,t_test_per_point);

    err_sq_mean=mean(err.^2,2);
    smse_curve=cumsum(err_sq_mean)./(1:N_eval)'/mean(Y_var_base);
    rmse_curve=sqrt(cumsum(err_sq_mean)./(1:N_eval)');

    save(fullfile(SaveFolder,sprintf('%s_tp_tr%d_mc%d.mat',Current_Method,tr_tag,seed)),...
        'smse','rmse','nlpd','t_train','t_test',...
        't_train_per_point','t_test_per_point',...
        'comm_train','comm_test','iter_converge',...
        'Current_Method','seed','train_ratio','smse_curve','rmse_curve');
end
end