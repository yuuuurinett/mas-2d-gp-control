function run_testpoint_dataset(DatasetName, CurrentMode, train_ratio, seed)
if nargin < 3, train_ratio = 0.4; end
if nargin < 4, seed = 1; end
rng(seed);
fprintf('\n[测试点] %s  seed=%d  tr=%.0f%%\n', DatasetName, seed, train_ratio*100);

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
num_train_samples=round(size(train_x,1)*train_ratio);
train_indices=randperm(size(train_x,1),num_train_samples);
X_train=train_x(train_indices,:); Y_train=train_y(train_indices,:);
test_indices=randperm(size(test_x,1));
X_test=test_x(test_indices,:); Y_test=test_y(test_indices,:);
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
[num_train,input_dim]=size(X_train); output_dim=size(Y_train,2);
num_eval=min(3000,size(X_test,1));
X_eval=X_test(1:num_eval,:); Y_eval=Y_test(1:num_eval,:);
Y_var_base=var(Y_eval,0,1);
fprintf('Train=%d Test=%d x=%d y=%d\n',num_train,num_eval,input_dim,output_dim);

%% 3. 分布式参数
num_agents=6; dac_step_size=0.01; dac_gain=10;
max_data_per_agent=min(floor(num_train/num_agents),3000);
MultiAgentSystem=Manipulator_2D_2DoF_SetMASTopology(num_agents,1);
Laplacian=MultiAgentSystem.Agent_Topology.LaplacianMatrix;
neighbor_count_per_agent=sum(Laplacian<0,2);
num_directed_neighbor_links=sum(neighbor_count_per_agent);
max_neighbor_count=max(neighbor_count_per_agent);
neighbor_mask=(Laplacian<0)|(Laplacian'<0);
neighbor_mask(1:size(neighbor_mask,1)+1:end)=false;
% ET参数（用于TP-DAC）
et_sigma=0.5; et_a=0.5/max_neighbor_count;
fprintf('ET参数: sigma=%.2f a=%.4f\n',et_sigma,et_a);

%% 4. 训练局部GP
t_start=tic;
local_gp_set=cell(num_agents,1);
for agent_idx=1:num_agents
    data_idx=(agent_idx-1)*max_data_per_agent+1:min(agent_idx*max_data_per_agent,num_train);
    local_gp_set{agent_idx}=LocalGP_MultiOutput(input_dim,output_dim,...
        max_data_per_agent,SigmaN,SigmaF,SigmaL);
    local_gp_set{agent_idx}.add_Alldata(X_train(data_idx,:)',Y_train(data_idx,:)');
    local_gp_set{agent_idx}.tau=1e-8; local_gp_set{agent_idx}.delta=0.01;
end
t_train_gp=toc(t_start);
fprintf('局部GP: %.4fs\n',t_train_gp);

%% 5. 方法列表
dac_methods={'poe','gpoe','moe','bcm','rbcm'};
ac_methods={'moe_ac','gpoe_ac','poe_ac','bcm_ac','rbcm_ac'};
if strcmpi(CurrentMode,'all'), AllModes=[dac_methods,ac_methods];
else, AllModes={lower(CurrentMode)}; end
SaveFolder=fullfile('Result','Dataset',DatasetName);
if ~exist(SaveFolder,'dir'), mkdir(SaveFolder); end
p_dim=2*output_dim; tr_tag=round(train_ratio*100);

%% 6. 批量预计算测试点局部预测
fprintf('[预计算] %d个测试点...\n',num_eval); tic;
local_mu_at_testpoints=zeros(num_agents,output_dim,num_eval);
local_var_at_testpoints=zeros(num_agents,output_dim,num_eval);
for agent_idx=1:num_agents
    [mu_batch,var_batch]=batch_predict_external(local_gp_set{agent_idx},X_eval',SigmaN,SigmaF);
    for test_idx=1:num_eval
        local_mu_at_testpoints(agent_idx,:,test_idx)=mu_batch(test_idx,:);
        local_var_at_testpoints(agent_idx,:,test_idx)=repmat(var_batch(test_idx),1,output_dim);
    end
end
precompute_time=toc;
fprintf('预计算完成: %.2fs\n',precompute_time);

%% 7. 主循环
for method_idx=1:numel(AllModes)
    current_method=AllModes{method_idx};
    base_method_name=strrep(lower(current_method),'_ac','');
    fprintf('\n[%d/%d] %s\n',method_idx,numel(AllModes),current_method);
    final_mu_pred=zeros(num_eval,output_dim);
    final_var_pred=zeros(num_eval,output_dim);
    tic;

    % 构建P矩阵
    P_info_matrix=zeros(p_dim,num_agents,num_eval);
    for test_idx=1:num_eval
        for agent_idx=1:num_agents
            local_mu_i=local_mu_at_testpoints(agent_idx,:,test_idx)';
            local_var_i=local_var_at_testpoints(agent_idx,:,test_idx)';
            for dim_idx=1:output_dim
                sv=max(local_var_i(dim_idx),SigmaN^2);
                bg=max(min(0.5*(log(prior_var)-log(sv)),10),eps);
                switch base_method_name
                    case 'moe'
                        P_info_matrix(2*dim_idx-1,agent_idx,test_idx)=num_agents*local_mu_i(dim_idx);
                        P_info_matrix(2*dim_idx,agent_idx,test_idx)=num_agents*(sv+local_mu_i(dim_idx)^2);
                    case 'gpoe'
                        P_info_matrix(2*dim_idx-1,agent_idx,test_idx)=num_agents*bg*local_mu_i(dim_idx)/sv;
                        P_info_matrix(2*dim_idx,agent_idx,test_idx)=num_agents*bg/sv;
                    case 'poe'
                        P_info_matrix(2*dim_idx-1,agent_idx,test_idx)=num_agents*local_mu_i(dim_idx)/sv;
                        P_info_matrix(2*dim_idx,agent_idx,test_idx)=num_agents/sv;
                    case 'bcm'
                        P_info_matrix(2*dim_idx-1,agent_idx,test_idx)=num_agents*local_mu_i(dim_idx)/sv;
                        P_info_matrix(2*dim_idx,agent_idx,test_idx)=num_agents/sv-(num_agents-1)/prior_var;
                    case 'rbcm'
                        P_info_matrix(2*dim_idx-1,agent_idx,test_idx)=num_agents*bg*local_mu_i(dim_idx)/sv;
                        P_info_matrix(2*dim_idx,agent_idx,test_idx)=num_agents*bg/sv+(1-num_agents*bg)/prior_var;
                end
            end
        end
    end

    Xi=P_info_matrix; dac_iter=0; max_iters=3000;
    Pi_scale=max(1,mean(P_info_matrix(:).^2));

    if ismember(lower(current_method),dac_methods)
        %% TP-DAC：有 reference signal
        % 动力学：dZeta/dt = Kappa * (Pi - Zeta) * L'
        % Zeta 从零开始，追踪 Pi 的全局平均
        Zeta=zeros(p_dim,num_agents,num_eval);
        while dac_iter<max_iters
            dac_iter=dac_iter+1; Zeta_prev=Zeta;
            diff=P_info_matrix-Zeta;
            for agent_idx=1:num_agents
                Zeta(:,agent_idx,:)=Zeta(:,agent_idx,:)+dac_step_size*dac_gain*...
                    sum(diff.*reshape(Laplacian(agent_idx,:),1,num_agents,1),2);
            end
            if max(abs(Zeta(:)-Zeta_prev(:)))<1e-5, break; end
        end
        iter_converge=dac_iter;
        comm_test=iter_converge;
        comm_train=0;
        % 结果从 Pi - Zeta 提取
        Xi=P_info_matrix-Zeta;
        fprintf('  [TP-DAC] 收敛步数:%d 通信轮数:%d\n',iter_converge,comm_test);

    else
        %% TP-AC：无reference signal，纯共识 dXi/dt = -L*Xi
        % Xi初始值=Pi，收敛到所有agent的Pi的平均
        while dac_iter<max_iters
            dac_iter=dac_iter+1; Xi_prev=Xi;
            L_Xi=zeros(size(Xi));
            for agent_idx=1:num_agents
                L_Xi(:,agent_idx,:)=sum(Xi.*reshape(Laplacian(agent_idx,:),1,num_agents,1),2);
            end
            for agent_idx=1:num_agents
                Xi(:,agent_idx,:)=Xi(:,agent_idx,:)-dac_step_size*dac_gain*L_Xi(:,agent_idx,:);
            end
            % 用agent间disagreement判断收敛（不用前后差，防止第一步就收敛）
            Xi_mean = mean(Xi, 2);
            disagreement = max(abs(Xi - Xi_mean), [], 'all');
            if disagreement < 1e-5, break; end
        end
        iter_converge=dac_iter;
        comm_test=iter_converge;
        comm_train=0;
        fprintf('  [TP-AC] 收敛步数:%d 通信轮数:%d\n',iter_converge,comm_test);
    end

    % 提取预测
    % TP-DAC: Xi = P - Zeta，取第一个agent（收敛后所有agent值相同）
    % TP-AC:  Xi收敛到Pi的全局平均，取第一个agent
    Xi_consensus=Xi;
    for test_idx=1:num_eval
        for dim_idx=1:output_dim
            xi1=Xi_consensus(2*dim_idx-1,1,test_idx);  % 取agent 1
            xi2=Xi_consensus(2*dim_idx,  1,test_idx);
            if ismember(base_method_name,{'gpoe','poe','bcm','rbcm'})
                final_mu_pred(test_idx,dim_idx)=xi1/max(abs(xi2),1e-4);
            else
                final_mu_pred(test_idx,dim_idx)=xi1/num_agents;
            end
        end
    end

    % 方差聚合
    for test_idx=1:num_eval
        for dim_idx=1:output_dim
            va=local_var_at_testpoints(:,dim_idx,test_idx);
            ba=max(0.5*(log(prior_var)-log(va)),eps);
            switch base_method_name
                case 'moe', final_var_pred(test_idx,dim_idx)=mean(va);
                case 'gpoe', prec=sum(ba./va); final_var_pred(test_idx,dim_idx)=1/max(prec,eps);
                case 'poe',  prec=sum(1./va);  final_var_pred(test_idx,dim_idx)=1/max(prec,eps);
                case 'bcm',  prec=sum(1./va)-(num_agents-1)/prior_var; final_var_pred(test_idx,dim_idx)=1/max(prec,eps);
                case 'rbcm', prec=sum(ba./va)+(1-sum(ba))/prior_var; final_var_pred(test_idx,dim_idx)=1/max(prec,eps);
            end
            final_var_pred(test_idx,dim_idx)=max(final_var_pred(test_idx,dim_idx),SigmaN^2);
        end
    end

    t_test_total=precompute_time/numel(AllModes)+toc;
    mu_pred=final_mu_pred.*repmat(Y_std,num_eval,1)+repmat(Y_mean,num_eval,1);
    var_pred=final_var_pred.*repmat(Y_std.^2,num_eval,1);
    prediction_error=Y_eval-mu_pred;
    smse=mean(mean(prediction_error.^2)./Y_var_base);
    rmse=mean(sqrt(mean(prediction_error.^2)));
    nlpd=mean(mean(0.5*(log(2*pi*var_pred)+prediction_error.^2./var_pred)));
    t_train_per_point=(t_train_gp/num_train)*1000;
    t_test_per_point=(t_test_total/num_eval)*1000;
    fprintf('  SMSE=%.4f RMSE=%.4f NLPD=%.4f Train:%.2fms/pt Test:%.2fms/pt\n',...
        smse,rmse,nlpd,t_train_per_point,t_test_per_point);

    err_sq_mean=mean(prediction_error.^2,2);
    smse_curve=cumsum(err_sq_mean)./(1:num_eval)'/mean(Y_var_base);
    rmse_curve=sqrt(cumsum(err_sq_mean)./(1:num_eval)');
    event_count_test=comm_test;
    trigger_count_per_agent=zeros(num_agents,num_eval);

    save(fullfile(SaveFolder,sprintf('%s_tp_tr%d_mc%d.mat',current_method,tr_tag,seed)),...
        'smse','rmse','nlpd','t_train_gp','t_test_total',...
        't_train_per_point','t_test_per_point',...
        'comm_train','comm_test','iter_converge',...
        'neighbor_count_per_agent','num_directed_neighbor_links','neighbor_mask',...
        'trigger_count_per_agent','event_count_test',...
        'current_method','seed','train_ratio','smse_curve','rmse_curve');
end
end