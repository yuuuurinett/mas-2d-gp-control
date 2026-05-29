function profile_sarcos_timing()
% 专门针对 SARCOS 的逐步计时分析脚本
% 运行后会生成一张分解图，精确定位每一步的耗时占比
% 用法: profile_sarcos_timing()

rng(1);

%% ---- 加载 & 归一化 (复用你的标准流程) ----
fprintf('正在加载 SARCOS...\n');
tr     = load(fullfile('SARCOS','SARCOS_train.mat'));
te     = load(fullfile('SARCOS','SARCOS_test.mat'));
hp_raw = load(fullfile('SARCOS','SARCOS_GP_Hyperparameter.mat'));

hp.SigmaF = mean(cell2mat(hp_raw.SigmaF_set));
hp.SigmaN = mean(cell2mat(hp_raw.SigmaN_set));
hp.SigmaL = mean(cell2mat(hp_raw.SigmaL_set'), 2);

train_x = tr.sarcos_inv(:,1:21);  train_y = tr.sarcos_inv(:,22:28);
test_x  = te.sarcos_inv_test(:,1:21);

% 按 40% 抽取
n_train = round(size(train_x,1)*0.4);
idx_tr  = randperm(size(train_x,1), n_train);
X_train = train_x(idx_tr,:);  Y_train = train_y(idx_tr,:);
X_test  = test_x(randperm(size(test_x,1)),:);

SigmaL = hp.SigmaL(:);
SigmaF = hp.SigmaF;
SigmaN = hp.SigmaN;

X_mean = mean(X_train,1);  X_std = std(X_train,0,1);  X_std(X_std==0)=1;
X_train = (X_train-X_mean)./X_std;
X_test  = (X_test -X_mean)./X_std;
SigmaL  = SigmaL ./ X_std(:);

Y_mean = mean(Y_train,1);  Y_std = std(Y_train,0,1);  Y_std(Y_std==0)=1;
Y_train = (Y_train-Y_mean)./Y_std;
SigmaF  = SigmaF / mean(Y_std);
SigmaN  = SigmaN / mean(Y_std);
prior_var = SigmaF^2;

AgentQuantity    = 6;
N_train          = size(X_train,1);
y_dim            = size(Y_train,2);   % = 7
x_dim            = size(X_train,2);   % = 21
MaxDataPerAgent  = min(floor(N_train/AgentQuantity), 3000);
NumInducingPoints = 2500;
N_eval_list      = [50, 100, 200, 500, 1000, 2000, 3000];
N_eval_max       = max(N_eval_list);

fprintf('N_train=%d  y_dim=%d  x_dim=%d  MaxPerAgent=%d\n', ...
    N_train, y_dim, x_dim, MaxDataPerAgent);

%% ---- Step 1: 局部 GP 训练 ----
fprintf('\n[Step 1] 训练 %d 个局部 GP...\n', AgentQuantity);
t1 = tic;
LocalGP_set = cell(AgentQuantity,1);
for n = 1:AgentQuantity
    idx = (n-1)*MaxDataPerAgent+1 : min(n*MaxDataPerAgent, N_train);
    LocalGP_set{n} = LocalGP_MultiOutput(x_dim,y_dim,MaxDataPerAgent,SigmaN,SigmaF,SigmaL);
    LocalGP_set{n}.add_Alldata(X_train(idx,:)', Y_train(idx,:)');
    LocalGP_set{n}.tau=1e-8;  LocalGP_set{n}.delta=0.01;
end
T_gp_train = toc(t1);
fprintf('  LocalGP 训练总时间: %.3fs  (%.3f ms/pt)\n', T_gp_train, T_gp_train/N_train*1000);

%% ---- Step 2: 单点预测速度 vs 批量预测速度 ----
fprintf('\n[Step 2] 测量单点 vs 批量预测速度...\n');
x_single = X_test(1,:)';
X_batch  = X_test(1:100,:)';

% 单点预测 (100次取均值)
times_single = zeros(1,20);
for k = 1:20
    t_tmp = tic;
    LocalGP_set{1}.predict(x_single);
    times_single(k) = toc(t_tmp)*1000;
end
T_single_ms = mean(times_single(5:end));  % 去掉前几次 JIT 预热

% 批量预测 100 点
t_tmp = tic;
LocalGP_set{1}.predict(X_batch);
T_batch_100_ms = toc(t_tmp)*1000;

fprintf('  单点预测 (1 agent):      %.3f ms\n', T_single_ms);
fprintf('  批量预测 100点 (1 agent): %.3f ms  (%.3f ms/pt)\n', ...
    T_batch_100_ms, T_batch_100_ms/100);
fprintf('  批量加速比: %.1fx\n', T_single_ms*100 / T_batch_100_ms);

%% ---- Step 3: 预计算时间随 N_eval 的扩展 ----
fprintf('\n[Step 3] 预计算时间 vs N_eval...\n');
T_precompute = zeros(size(N_eval_list));

for ki = 1:numel(N_eval_list)
    N_e = N_eval_list(ki);
    X_e = X_test(1:N_e,:)';
    t_tmp = tic;
    for n = 1:AgentQuantity
        try
            [mn, vn] = LocalGP_set{n}.predict(X_e);
            % 消耗结果防止被优化掉
            if isempty(mn) || isempty(vn), disp('empty'); end
        catch
            for tt = 1:N_e
                LocalGP_set{n}.predict(X_e(:,tt));
            end
        end
    end
    T_precompute(ki) = toc(t_tmp);
    fprintf('  N_eval=%4d: %.3fs  (%.3f ms/pt)\n', N_e, T_precompute(ki), T_precompute(ki)/N_e*1000);
end

%% ---- Step 4: ET 迭代时间 (固定 N_eval=3000，测不同 y_dim) ----
fprintf('\n[Step 4] ET 迭代耗时 vs 状态维度 p_dim...\n');

% 用随机 P 矩阵模拟，测量纯 ET 迭代成本
MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology(AgentQuantity,1);
L = MultiAgentSystem.Agent_Topology.LaplacianMatrix;
N_degree = sum(L < 0, 2);
N_max    = max(N_degree);
sigma_i  = 0.5;
a_param  = 0.9 / N_max;
Kappa_P  = 10;  t_step = 0.01;

ydim_list  = [1, 2, 3, 5, 7];  % p_dim = 2*y_dim
T_et_iter  = zeros(size(ydim_list));
N_et_eval  = 500;  % 固定测试点数，只测维度影响

for ki = 1:numel(ydim_list)
    yd  = ydim_list(ki);
    pd  = 2*yd;
    Pi  = randn(pd, AgentQuantity, N_et_eval);

    Zeta   = zeros(pd, AgentQuantity, N_et_eval);
    Zeta_k = zeros(pd, AgentQuantity, N_et_eval);
    trigger_count_set = zeros(AgentQuantity,1);

    t_tmp = tic;
    for iter = 1:3000
        Zeta_prev = Zeta;
        diff_k = Pi - Zeta_k;
        for n = 1:AgentQuantity
            Zeta(:,n,:) = Zeta(:,n,:) + t_step * Kappa_P * ...
                sum(diff_k .* reshape(L(n,:),1,AgentQuantity,1), 2);
        end
        for n = 1:AgentQuantity
            e_i = Zeta_k(:,n,:) - Zeta(:,n,:);
            z_i = sum(reshape(L(n,:),1,AgentQuantity,1) .* Zeta_k, 2);
            rho_i     = max(sum(e_i.^2,1));
            norm_z_sq = max(sum(z_i.^2,1));
            N_i       = N_degree(n);
            rho_bar_i = (sigma_i * a_param * (1-a_param*N_i)/N_i) * norm_z_sq;
            if rho_i > rho_bar_i || iter==1
                Zeta_k(:,n,:) = Zeta(:,n,:);
                trigger_count_set(n) = trigger_count_set(n)+1;
            end
        end
        if max(abs(Zeta(:)-Zeta_prev(:))) < 1e-4
            break;
        end
    end
    T_et_iter(ki) = toc(t_tmp);
    fprintf('  y_dim=%d  p_dim=%2d: ET %.3fs  (%.3f ms/pt, %d步)\n', ...
        yd, pd, T_et_iter(ki), T_et_iter(ki)/N_et_eval*1000, iter);
end

%% ---- Step 5: IP 架构 MaskedGP 预测时间 ----
fprintf('\n[Step 5] IP 架构 MaskedGP 测试时间...\n');
idx_ind  = randperm(N_train, NumInducingPoints);
IP_coords = X_train(idx_ind,:)';
phi_dummy = randn(y_dim, NumInducingPoints);

MaskedGP = LocalGP_MultiOutput(x_dim,y_dim,NumInducingPoints,1e-6,SigmaF,SigmaL);
MaskedGP.add_Alldata(IP_coords, phi_dummy);

X_eval_3k = X_test(1:3000,:)';
t_tmp = tic;
Num_Ind   = MaskedGP.DataQuantity;
Alpha_Vec = MaskedGP.alpha(1:Num_Ind,:);
Chol_L    = MaskedGP.L(1:Num_Ind,1:Num_Ind);
K_star    = MaskedGP.kernel(MaskedGP.X(:,1:Num_Ind), X_eval_3k);
mu_norm   = (Alpha_Vec' * K_star)';
V_mat     = Chol_L \ K_star;
var_norm  = max(SigmaF^2 - sum(V_mat.^2,1)', SigmaN^2);
T_ip_test = toc(t_tmp);
fprintf('  IP MaskedGP 预测 3000点: %.3fs  (%.3f ms/pt)\n', T_ip_test, T_ip_test/3000*1000);

%% ---- 汇总打印 ----
fprintf('\n========== SARCOS 耗时汇总 ==========\n');
fprintf('LocalGP 训练 (6 agents):     %.3f ms/pt\n', T_gp_train/N_train*1000);
fprintf('预计算 6×LocalGP (3000pts):  %.3f ms/pt\n', T_precompute(end)/3000*1000);
fprintf('ET 迭代 (y_dim=7, 500pts):   %.3f ms/pt\n', T_et_iter(end)/N_et_eval*1000);
fprintf('IP MaskedGP 测试 (3000pts):  %.3f ms/pt\n', T_ip_test/3000*1000);
fprintf('ET 迭代占预计算比例:          %.1f%%\n', ...
    T_et_iter(end)/N_et_eval / (T_precompute(end)/3000) * 100);
fprintf('=====================================\n');

%% ---- 绘图 ----
fig = figure('Name','SARCOS Timing Profiling','Position',[100 100 1100 800]);

% --- 图1: 各步骤绝对耗时 (ms/pt) ---
subplot(2,3,1);
steps = {'LocalGP\n训练', '预计算\n(6×GP)', 'ET迭代\n(y=7)', 'IP测试\n(MaskedGP)'};
vals  = [T_gp_train/N_train*1000, ...
         T_precompute(end)/3000*1000, ...
         T_et_iter(end)/N_et_eval*1000, ...
         T_ip_test/3000*1000];
colors_bar = [0.32 0.45 0.72; 0.60 0.37 0.18; 0.80 0.30 0.20; 0.20 0.55 0.40];
b = bar(vals, 'FaceColor','flat');
b.CData = colors_bar;
set(gca,'XTickLabel',{'LocalGP训练','预计算(6GP)','ET迭代(y=7)','IP测试'}, ...
    'XTickLabelRotation',20,'FontSize',10);
ylabel('ms / point'); title('各步骤单点耗时');
grid on; box off;
for i=1:4
    text(i, vals(i)+max(vals)*0.02, sprintf('%.3f',vals(i)), ...
        'HorizontalAlignment','center','FontSize',9);
end

% --- 图2: 预计算时间 vs N_eval (线性 & 对数) ---
subplot(2,3,2);
plot(N_eval_list, T_precompute*1000, 'o-', 'LineWidth',1.5, 'Color',[0.32 0.45 0.72], ...
    'MarkerFaceColor',[0.32 0.45 0.72]);
xlabel('N_{eval}'); ylabel('总时间 (ms)'); title('预计算时间 vs 测试点数');
grid on; box off;
% 拟合斜率，判断是否线性
p = polyfit(N_eval_list, T_precompute*1000, 1);
hold on;
plot(N_eval_list, polyval(p,N_eval_list), '--', 'Color',[0.6 0.6 0.6]);
legend('实测','线性拟合','Location','northwest','FontSize',9);

% --- 图3: ET 迭代时间 vs y_dim ---
subplot(2,3,3);
plot(ydim_list, T_et_iter/N_et_eval*1000, 's-', 'LineWidth',1.5, ...
    'Color',[0.80 0.30 0.20], 'MarkerFaceColor',[0.80 0.30 0.20]);
xlabel('y\_dim (输出维度)'); ylabel('ms / point'); title('ET迭代耗时 vs 输出维度');
xline(7, '--k', 'SARCOS y=7', 'LabelVerticalAlignment','bottom','FontSize',9);
grid on; box off;

% --- 图4: 单点 vs 批量预测对比 ---
subplot(2,3,4);
labels = {'单点预测\n×100次', '批量预测\n100点'};
vals4  = [T_single_ms*100, T_batch_100_ms];
b4 = bar(vals4, 'FaceColor','flat');
b4.CData = [0.60 0.37 0.18; 0.20 0.55 0.40];
set(gca,'XTickLabel',{'单点×100','批量100点'},'FontSize',10);
ylabel('总时间 (ms)'); title(sprintf('批量加速比: %.1fx', T_single_ms*100/T_batch_100_ms));
grid on; box off;
for i=1:2
    text(i, vals4(i)+max(vals4)*0.02, sprintf('%.1fms',vals4(i)), ...
        'HorizontalAlignment','center','FontSize',9);
end

% --- 图5: 耗时饼图 (TP-DAC 测试阶段) ---
subplot(2,3,5);
% TP-DAC 测试 = 预计算 + ET迭代
t_precomp_3k = T_precompute(end)*1000;
t_et_3k      = T_et_iter(end) / N_et_eval * 3000 * 1000;  % 外推到3000点
pie_vals = [t_precomp_3k, t_et_3k];
pie_labels = {sprintf('预计算\n%.0fms',t_precomp_3k), sprintf('ET迭代\n%.0fms',t_et_3k)};
p5 = pie(pie_vals);
for i=1:2:length(p5)
    p5(i).FaceColor = colors_bar(i+1,:);
end
for i=2:2:length(p5)
    p5(i).FontSize = 9;
    p5(i).String   = pie_labels{(i+1)/2};
end
title('TP-DAC 测试耗时构成 (3000pts)');

% --- 图6: 与其他数据集对比 (预计算时间) ---
subplot(2,3,6);
% 用你结果表里的 CEN Test 数据作为预计算代理
datasets = {'KIN40K','POL','PUMADYN','SARCOS'};
cen_test = [1.80, 1.88, 0.99, 28.99];   % ms/pt，来自你的结果表
bar(cen_test, 'FaceColor','flat', 'CData', ...
    [0.32 0.45 0.72; 0.32 0.45 0.72; 0.32 0.45 0.72; 0.80 0.30 0.20]);
set(gca,'XTickLabel',datasets,'FontSize',10);
ylabel('ms / point'); title('CEN Test 耗时对比（基准参考）');
grid on; box off;
for i=1:4
    text(i, cen_test(i)+0.5, sprintf('%.2f',cen_test(i)), ...
        'HorizontalAlignment','center','FontSize',9);
end

sgtitle('SARCOS 耗时 Profiling 分析', 'FontSize',13, 'FontWeight','bold');

% 保存图片
saveas(fig, fullfile('Result','sarcos_timing_profile.png'));
fprintf('\n图已保存至 Result/sarcos_timing_profile.png\n');
end