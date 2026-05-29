%% main_BatchTest
clc; close all;

%% 配置
CurrentMode = 'all';   % 'poe'|'gpoe'|'moe'|'bcm'|'rbcm' 或 'all'
TestType    = 'all';   % 'test'|'inducing'|'cen'|'nbr'|'all'

%% 路径
SaveFolder_Test     = fullfile('Result', 'Test_Point');
SaveFolder_Inducing = fullfile('Result', 'Inducing_Point');
SaveFolder_CEN      = fullfile('Result', 'CEN');
SaveFolder_NBR      = fullfile('Result', 'NBR');
for f = {SaveFolder_Test, SaveFolder_Inducing, SaveFolder_CEN, SaveFolder_NBR}
    if ~exist(f{1},'dir'), mkdir(f{1}); end
end

%% 模式列表
Modes_dac      = {'poe','gpoe','moe','bcm','rbcm'};
Modes_ac       = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};
Modes_baseline = {'local','exact'};

if strcmpi(CurrentMode,'all')
    AllModes_dac = Modes_dac;
else
    AllModes_dac = {CurrentMode};
end

%% 1. 运行仿真
run_inducing = ismember(lower(TestType), {'inducing','all'});
run_test     = ismember(lower(TestType), {'test','all'});
run_cen      = ismember(lower(TestType), {'cen','all'});
run_nbr      = ismember(lower(TestType), {'nbr','all'});

if run_inducing
    AllModes_ind = [Modes_dac, Modes_ac, Modes_baseline];
    fprintf('\n======== 诱导点聚合 (IP-DAC + IP-AC) ========\n');
    for m = 1:numel(AllModes_ind)
        fprintf('[%d/%d] %s\n', m, numel(AllModes_ind), AllModes_ind{m});
        run_simulation_inducing_point(AllModes_ind{m}, SaveFolder_Inducing, AllModes_ind{m});
    end
end

if run_test
    AllModes_test = [Modes_dac, Modes_ac, Modes_baseline];
    fprintf('\n======== 测试点聚合 (TP-DAC + TP-AC) ========\n');
    for m = 1:numel(AllModes_test)
        fprintf('[%d/%d] %s\n', m, numel(AllModes_test), AllModes_test{m});
        run_simulation_test_point(AllModes_test{m}, SaveFolder_Test, AllModes_test{m});
    end
end

if run_cen
    fprintf('\n======== 集中式聚合 (CEN) ========\n');
    for m = 1:numel(AllModes_dac)
        fname = sprintf('cen_%s', AllModes_dac{m});
        fprintf('[%d/%d] %s\n', m, numel(AllModes_dac), fname);
        run_simulation_cen(AllModes_dac{m}, SaveFolder_CEN, fname);
    end
end

if run_nbr
    fprintf('\n======== 邻域聚合 (NBR) ========\n');
    for m = 1:numel(AllModes_dac)
        fname = sprintf('nbr_%s', AllModes_dac{m});
        fprintf('[%d/%d] %s\n', m, numel(AllModes_dac), fname);
        run_simulation_nbr(AllModes_dac{m}, SaveFolder_NBR, fname);
    end
end

%% 2. 加载结果
% 读取时间轴
ref_candidates = {
    fullfile(SaveFolder_Inducing,'poe.mat');
    fullfile(SaveFolder_Test,    'poe.mat');
    fullfile(SaveFolder_CEN,     'cen_poe.mat');
    fullfile(SaveFolder_NBR,     'nbr_poe.mat')};
t_set = [];
for i = 1:numel(ref_candidates)
    if exist(ref_candidates{i},'file')
        tmp = load(ref_candidates{i},'t_set');
        t_set = tmp.t_set; break;
    end
end
if isempty(t_set), fprintf('无结果文件，退出。\n'); return; end
N = numel(t_set);
load_err = @(folder,fname) load_tracking_error(folder,fname,N);

% IP-DAC
Err_IP_DAC = nan(numel(Modes_dac),N);
for m=1:numel(Modes_dac)
    Err_IP_DAC(m,:) = load_err(SaveFolder_Inducing, Modes_dac{m});
end

% IP-AC
Err_IP_AC = nan(numel(Modes_ac),N);
for m=1:numel(Modes_ac)
    Err_IP_AC(m,:) = load_err(SaveFolder_Inducing, Modes_ac{m});
end

% TP-DAC
Err_TP_DAC = nan(numel(Modes_dac),N);
for m=1:numel(Modes_dac)
    Err_TP_DAC(m,:) = load_err(SaveFolder_Test, Modes_dac{m});
end

% TP-AC
Err_TP_AC = nan(numel(Modes_ac),N);
for m=1:numel(Modes_ac)
    Err_TP_AC(m,:) = load_err(SaveFolder_Test, Modes_ac{m});
end

% CEN
Err_CEN = nan(numel(Modes_dac),N);
for m=1:numel(Modes_dac)
    Err_CEN(m,:) = load_err(SaveFolder_CEN, sprintf('cen_%s',Modes_dac{m}));
end

% NBR
Err_NBR = nan(numel(Modes_dac),N);
for m=1:numel(Modes_dac)
    Err_NBR(m,:) = load_err(SaveFolder_NBR, sprintf('nbr_%s',Modes_dac{m}));
end

% Baseline
Err_Local = load_err(SaveFolder_Test, 'local');
Err_Exact = load_err(SaveFolder_Test, 'exact');

%% 3. 打印汇总表格
fprintf('\n%s\n', repmat('=',1,80));
fprintf('  Final Tracking Error ||e(T)||\n');
fprintf('  %-10s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n', ...
    'Method','IP-DAC','IP-AC','TP-DAC','TP-AC','CEN','NBR');
fprintf('  %s\n', repmat('-',1,76));
for m=1:numel(Modes_dac)
    fprintf('  %-10s  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-8.4f\n', ...
        Modes_dac{m}, ...
        Err_IP_DAC(m,end), Err_IP_AC(m,end), ...
        Err_TP_DAC(m,end), Err_TP_AC(m,end), ...
        Err_CEN(m,end),    Err_NBR(m,end));
end
fprintf('  %s\n', repmat('-',1,76));
fprintf('  %-10s  %-8s  %-8s  %-8.4f  %-8s  %-8s  %-8s\n', ...
    'local','-','-',Err_Local(end),'-','-','-');
fprintf('  %-10s  %-8s  %-8s  %-8.4f  %-8s  %-8s  %-8s\n', ...
    'exact','-','-',Err_Exact(end),'-','-','-');
fprintf('%s\n\n', repmat('=',1,80));

%% 4. 绘图：每种聚合方法一张图，对比所有框架
lw = 1.5;
for m=1:numel(Modes_dac)
    ac_name = [Modes_dac{m},'_ac'];
    ac_idx  = find(strcmp(Modes_ac, ac_name));

    fig = figure('Name',upper(Modes_dac{m}),'Color','w');
    hold on; grid on; box on;
    set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');

    % IP 框架
    if ~all(isnan(Err_IP_DAC(m,:)))
        plot(t_set, Err_IP_DAC(m,:), 'b-',  'LineWidth',lw, 'DisplayName','IP-DAC');
    end
    if ~isempty(ac_idx) && ~all(isnan(Err_IP_AC(ac_idx,:)))
        plot(t_set, Err_IP_AC(ac_idx,:), 'b--', 'LineWidth',lw, 'DisplayName','IP-AC');
    end

    % TP 框架
    if ~all(isnan(Err_TP_DAC(m,:)))
        plot(t_set, Err_TP_DAC(m,:), 'r-',  'LineWidth',lw, 'DisplayName','TP-DAC');
    end
    if ~isempty(ac_idx) && ~all(isnan(Err_TP_AC(ac_idx,:)))
        plot(t_set, Err_TP_AC(ac_idx,:), 'r--', 'LineWidth',lw, 'DisplayName','TP-AC');
    end

    % CEN / NBR
    if ~all(isnan(Err_CEN(m,:)))
        plot(t_set, Err_CEN(m,:), 'g-',  'LineWidth',lw, 'DisplayName','CEN');
    end
    if ~all(isnan(Err_NBR(m,:)))
        plot(t_set, Err_NBR(m,:), 'm-',  'LineWidth',lw, 'DisplayName','NBR');
    end

    % Baseline
    if ~all(isnan(Err_Local))
        plot(t_set, Err_Local, 'k--', 'LineWidth',1, 'DisplayName','Local');
    end
    if ~all(isnan(Err_Exact))
        plot(t_set, Err_Exact, 'k-',  'LineWidth',1, 'DisplayName','Exact');
    end

    ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
    xlabel('$t$','Interpreter','latex','FontSize',13);
    title(sprintf('Tracking Error - %s', upper(Modes_dac{m})), ...
        'FontSize',12,'FontName','Times New Roman');
    legend('Location','northeast','FontSize',9,'NumColumns',2);
    xlim([t_set(1), t_set(end)]);
    saveas(fig, fullfile('Result', sprintf('Comparison_%s.png', upper(Modes_dac{m}))));
end

fprintf('完成。\n');


%% 辅助函数
function err = load_tracking_error(folder, fname, N)
    fpath = fullfile(folder, [fname, '.mat']);
    if exist(fpath,'file')
        d = load(fpath,'TrackingError_vector');
        len = min(numel(d.TrackingError_vector), N);
        err = nan(1,N);
        err(1:len) = d.TrackingError_vector(1:len);
    else
        err = nan(1,N);
    end
end