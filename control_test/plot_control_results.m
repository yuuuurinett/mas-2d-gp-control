%% plot_control_results.m
clc; close all;

%% ===== 配置 =====
method        = 'poe';
SaveFolder    = fullfile('Result', 'Inducing_Point');
FigFolder     = fullfile('Result', 'Figures');
if ~exist(FigFolder, 'dir'), mkdir(FigFolder); end

fname_form   = fullfile(SaveFolder, sprintf('%s_formation.mat',   method));
fname_noform = fullfile(SaveFolder, sprintf('%s_noformation.mat', method));

%% ===== 加载数据 =====
if ~exist(fname_form,   'file'), error('找不到文件: %s', fname_form);   end
if ~exist(fname_noform, 'file'), error('找不到文件: %s', fname_noform); end

d_form   = load(fname_form);
d_noform = load(fname_noform);

t_set         = d_form.t_set;
AgentQuantity = 6;
x_dim         = 4;

%% ===== 图1：每个 agent — Tracking Error，有/无 formation 对比 =====
% Tracking error: ||ϑ_i(t)|| = ||x_i(t) - s_i(t) - x_l(t)||
for n = 1:AgentQuantity
    idx = (n-1)*x_dim + (1:x_dim);

    te_form   = sqrt(sum(d_form.vartheta_all_set(idx,:).^2,   1));
    te_noform = sqrt(sum(d_noform.vartheta_all_set(idx,:).^2, 1));

    fig = figure('Name', sprintf('Agent%d_Tracking', n), ...
        'Color', 'w', 'Position', [100, 100, 800, 350]);

    hold on; grid on; box on;
    plot(t_set, te_form,   'b-',  'LineWidth', 1.5, 'DisplayName', 'With Formation');
    plot(t_set, te_noform, 'r--', 'LineWidth', 1.5, 'DisplayName', 'No Formation');

    ylabel({'Tracking Error', '$\|\vartheta_i(t)\| = \|x_i - s_i - x_l\|$'}, ...
        'Interpreter', 'latex', 'FontSize', 11);
    xlabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 11);
    title(sprintf('Agent %d — %s: Tracking Error (With vs Without Formation)', ...
        n, upper(method)), 'FontSize', 11, 'FontName', 'Times New Roman');
    legend('Location', 'northeast', 'FontSize', 10);
    set(gca, 'FontSize', 10, 'FontName', 'Times New Roman');
    xlim([t_set(1), t_set(end)]);

    saveas(fig, fullfile(FigFolder, sprintf('%s_agent%d_tracking.png',  method, n)));
    savefig(fig, fullfile(FigFolder, sprintf('%s_agent%d_tracking.fig', method, n)));
end

%% ===== 图2：每个 agent — GP Prediction vs True Unknown Dynamics（有 formation）=====
% 两个维度分开画，直观展示 GP 预测和真实值的吻合程度
for n = 1:AgentQuantity
    f_true = squeeze(d_form.f_true_all_set(:, n, :));  % [2 × T]
    f_hat  = squeeze(d_form.f_hat_all_set(:,  n, :));  % [2 × T]

    fig = figure('Name', sprintf('Agent%d_Prediction', n), ...
        'Color', 'w', 'Position', [100, 100, 800, 550]);

    % --- 上图：第1维 ---
    subplot(2,1,1);
    hold on; grid on; box on;
    plot(t_set, f_true(1,:), 'k-',  'LineWidth', 1.5, ...
        'DisplayName', 'True $f_{i,1}$');
    plot(t_set, f_hat(1,:),  'b--', 'LineWidth', 1.5, ...
        'DisplayName', 'GP Prediction $\hat{f}_{i,1}$');
    ylabel('$f_{i,1}(t)$ (Joint 1)', 'Interpreter', 'latex', 'FontSize', 11);
    title(sprintf('Agent %d — %s: GP Prediction vs True Unknown Dynamics (With Formation)', ...
        n, upper(method)), 'FontSize', 11, 'FontName', 'Times New Roman');
    legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);
    set(gca, 'FontSize', 10, 'FontName', 'Times New Roman');
    xlim([t_set(1), t_set(end)]);

    % --- 下图：第2维 ---
    subplot(2,1,2);
    hold on; grid on; box on;
    plot(t_set, f_true(2,:), 'k-',  'LineWidth', 1.5, ...
        'DisplayName', 'True $f_{i,2}$');
    plot(t_set, f_hat(2,:),  'b--', 'LineWidth', 1.5, ...
        'DisplayName', 'GP Prediction $\hat{f}_{i,2}$');
    ylabel('$f_{i,2}(t)$ (Joint 2)', 'Interpreter', 'latex', 'FontSize', 11);
    xlabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 11);
    legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);
    set(gca, 'FontSize', 10, 'FontName', 'Times New Roman');
    xlim([t_set(1), t_set(end)]);

    %saveas(fig, fullfile(FigFolder, sprintf('%s_agent%d_prediction.png',  method, n)));
    savefig(fig, fullfile(FigFolder, sprintf('%s_agent%d_prediction.fig', method, n)));
end

fprintf('完成，共生成 %d 张图，已保存至 %s\n', AgentQuantity*2, FigFolder);