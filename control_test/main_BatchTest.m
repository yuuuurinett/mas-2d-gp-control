%% main_BatchTest
clc; close all;

%% 配置
CurrentMode  = 'all';   % 'poe'|'gpoe'|'moe'|'bcm'|'rbcm' 或 'all'
TestType     = 'inducing';   % 'test'|'inducing'|'cen'|'nbr'|'all'
use_formation = true;   % true = 有formation，false = 无formation

% 文件名后缀，方便两次结果共存
form_tag = 'formation';
if ~use_formation, form_tag = 'noformation'; end

%% 路径
SaveFolder_Test     = fullfile('Result', 'Test_Point');
SaveFolder_Inducing = fullfile('Result', 'Inducing_Point');
SaveFolder_CEN      = fullfile('Result', 'CEN');
SaveFolder_NBR      = fullfile('Result', 'NBR');
SaveFolder_Figures  = fullfile('Result', 'Figures');
for f = {SaveFolder_Test, SaveFolder_Inducing, SaveFolder_CEN, SaveFolder_NBR, SaveFolder_Figures}
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
%run_inducing = ismember(lower(TestType), {'inducing','all'});
%run_test     = ismember(lower(TestType), {'test','all'});
%run_cen      = ismember(lower(TestType), {'cen','all'});
%run_nbr      = ismember(lower(TestType), {'nbr','all'});

run_inducing = false;
run_test     = false;
run_cen      = false;
run_nbr      = false;

if run_inducing
    AllModes_ind = [Modes_dac, Modes_ac, Modes_baseline];
    fprintf('\n======== 诱导点聚合 (IP-DAC + IP-AC) [%s] ========\n', form_tag);
    for m = 1:numel(AllModes_ind)
        fname = sprintf('%s_%s', AllModes_ind{m}, form_tag);
        fprintf('[%d/%d] %s\n', m, numel(AllModes_ind), fname);
        run_simulation_inducing_point(AllModes_ind{m}, SaveFolder_Inducing, fname, use_formation);
    end
end

if run_test
    AllModes_test = [Modes_dac, Modes_ac, Modes_baseline];
    fprintf('\n======== 测试点聚合 (TP-DAC + TP-AC) [%s] ========\n', form_tag);
    for m = 1:numel(AllModes_test)
        fname = sprintf('%s_%s', AllModes_test{m}, form_tag);
        fprintf('[%d/%d] %s\n', m, numel(AllModes_test), fname);
        run_simulation_test_point(AllModes_test{m}, SaveFolder_Test, fname, use_formation);
    end
end

if run_cen
    fprintf('\n======== 集中式聚合 (CEN) [%s] ========\n', form_tag);
    for m = 1:numel(AllModes_dac)
        fname = sprintf('cen_%s_%s', AllModes_dac{m}, form_tag);
        fprintf('[%d/%d] %s\n', m, numel(AllModes_dac), fname);
        run_simulation_cen(AllModes_dac{m}, SaveFolder_CEN, fname, use_formation);
    end
end

if run_nbr
    fprintf('\n======== 邻域聚合 (NBR) [%s] ========\n', form_tag);
    for m = 1:numel(AllModes_dac)
        fname = sprintf('nbr_%s_%s', AllModes_dac{m}, form_tag);
        fprintf('[%d/%d] %s\n', m, numel(AllModes_dac), fname);
        run_simulation_nbr(AllModes_dac{m}, SaveFolder_NBR, fname, use_formation);
    end
end

%% 2. 加载结果（用 formation 版本的时间轴）
ref_candidates = {
    fullfile(SaveFolder_Inducing, sprintf('poe_%s.mat', form_tag));
    fullfile(SaveFolder_Test,     sprintf('poe_%s.mat', form_tag));
    fullfile(SaveFolder_CEN,      sprintf('cen_poe_%s.mat', form_tag));
    fullfile(SaveFolder_NBR,      sprintf('nbr_poe_%s.mat', form_tag))};
t_set = [];
for i = 1:numel(ref_candidates)
    if exist(ref_candidates{i},'file')
        tmp = load(ref_candidates{i},'t_set');
        t_set = tmp.t_set; break;
    end
end
if isempty(t_set), fprintf('无结果文件，退出。\n'); return; end
N = numel(t_set);
load_err = @(folder,fname) load_tracking_error(folder, fname, N);

% IP-DAC / IP-AC
Err_IP_DAC = nan(numel(Modes_dac),N);
Err_IP_AC  = nan(numel(Modes_ac),N);
for m=1:numel(Modes_dac)
    Err_IP_DAC(m,:) = load_err(SaveFolder_Inducing, sprintf('%s_%s', Modes_dac{m}, form_tag));
end
for m=1:numel(Modes_ac)
    Err_IP_AC(m,:)  = load_err(SaveFolder_Inducing, sprintf('%s_%s', Modes_ac{m}, form_tag));
end

% TP-DAC / TP-AC
Err_TP_DAC = nan(numel(Modes_dac),N);
Err_TP_AC  = nan(numel(Modes_ac),N);
for m=1:numel(Modes_dac)
    Err_TP_DAC(m,:) = load_err(SaveFolder_Test, sprintf('%s_%s', Modes_dac{m}, form_tag));
end
for m=1:numel(Modes_ac)
    Err_TP_AC(m,:)  = load_err(SaveFolder_Test, sprintf('%s_%s', Modes_ac{m}, form_tag));
end

% CEN / NBR
Err_CEN = nan(numel(Modes_dac),N);
Err_NBR = nan(numel(Modes_dac),N);
for m=1:numel(Modes_dac)
    Err_CEN(m,:) = load_err(SaveFolder_CEN, sprintf('cen_%s_%s', Modes_dac{m}, form_tag));
    Err_NBR(m,:) = load_err(SaveFolder_NBR, sprintf('nbr_%s_%s', Modes_dac{m}, form_tag));
end

% Baseline
Err_Local = load_err(SaveFolder_Test, sprintf('local_%s', form_tag));
Err_Exact = load_err(SaveFolder_Test, sprintf('exact_%s', form_tag));

% IP-DAC ET 通信量
% Dataset-consistent definition:
% average number of triggered events per agent per inducing point
% Comm = mean(total_trigger_count) / NumInducingPoints
Comm_IP_DAC = nan(numel(Modes_dac), 1);
for m = 1:numel(Modes_dac)
    fpath = fullfile(SaveFolder_Inducing, sprintf('%s_%s.mat', Modes_dac{m}, form_tag));
    if exist(fpath, 'file')
        d = load(fpath);

        if isfield(d, 'trigger_count_per_agent_point')
            % Preferred: already saved by the updated simulation script
            Comm_IP_DAC(m) = d.trigger_count_per_agent_point;

        elseif isfield(d, 'total_trigger_count') && isfield(d, 'NumInducingPoints')
            % Backward-compatible fallback for older result files
            Comm_IP_DAC(m) = mean(d.total_trigger_count) / d.NumInducingPoints;

        elseif isfield(d, 'total_trigger_count') && isfield(d, 'P_inducing')
            % Fallback if NumInducingPoints was not saved, but P_inducing exists
            Comm_IP_DAC(m) = mean(d.total_trigger_count) / size(d.P_inducing, 3);

        else
            Comm_IP_DAC(m) = NaN;
        end
    end
end

%% 3. 打印汇总表格
fprintf('\n%s\n', repmat('=',1,90));
fprintf('  Final Tracking Error ||e(T)||  [%s]\n', form_tag);
fprintf('  %-10s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-12s\n', ...
    'Method','IP-DAC','IP-AC','TP-DAC','TP-AC','CEN','NBR','IP-DAC Comm/pt');
fprintf('  %s\n', repmat('-',1,86));
for m=1:numel(Modes_dac)
    if ~isnan(Comm_IP_DAC(m))
        comm_str = sprintf('%.2f', Comm_IP_DAC(m));
    else
        comm_str = '-';
    end
    fprintf('  %-10s  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-12s\n', ...
        Modes_dac{m}, ...
        Err_IP_DAC(m,end), Err_IP_AC(m,end), ...
        Err_TP_DAC(m,end), Err_TP_AC(m,end), ...
        Err_CEN(m,end),    Err_NBR(m,end), ...
        comm_str);
end
fprintf('%s\n\n', repmat('=',1,90));

% Local 和 Exact 单独打印
fprintf('  %-10s  %-8.4f\n', 'Local', Err_Local(end));
fprintf('  %-10s  %-8.4f\n', 'Exact', Err_Exact(end));
fprintf('%s\n\n', repmat('=',1,90));

%% 4. 绘图
lw = 1.5;
for m=1:numel(Modes_dac)
    ac_name = [Modes_dac{m},'_ac'];
    ac_idx  = find(strcmp(Modes_ac, ac_name));

    fig = figure('Name',upper(Modes_dac{m}),'Color','w');
    hold on; grid on; box on;
    set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');

    if ~all(isnan(Err_IP_DAC(m,:)))
        plot(t_set, Err_IP_DAC(m,:), 'b-',  'LineWidth',lw, 'DisplayName','IP-DAC');
    end
    if ~isempty(ac_idx) && ~all(isnan(Err_IP_AC(ac_idx,:)))
        plot(t_set, Err_IP_AC(ac_idx,:), 'b--', 'LineWidth',lw, 'DisplayName','IP-AC');
    end
    if ~all(isnan(Err_TP_DAC(m,:)))
        plot(t_set, Err_TP_DAC(m,:), 'r-',  'LineWidth',lw, 'DisplayName','TP-DAC');
    end
    if ~isempty(ac_idx) && ~all(isnan(Err_TP_AC(ac_idx,:)))
        plot(t_set, Err_TP_AC(ac_idx,:), 'r--', 'LineWidth',lw, 'DisplayName','TP-AC');
    end
    if ~all(isnan(Err_CEN(m,:)))
        plot(t_set, Err_CEN(m,:), 'g-',  'LineWidth',lw, 'DisplayName','CEN');
    end
    if ~all(isnan(Err_NBR(m,:)))
        plot(t_set, Err_NBR(m,:), 'm-',  'LineWidth',lw, 'DisplayName','NBR');
    end
    if ~all(isnan(Err_Local))
        plot(t_set, Err_Local, 'k--', 'LineWidth',1, 'DisplayName','Local');
    end
    if ~all(isnan(Err_Exact))
        plot(t_set, Err_Exact, 'k-',  'LineWidth',1, 'DisplayName','Exact');
    end

    ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
    xlabel('$t$ (s)','Interpreter','latex','FontSize',13);
    title(sprintf('Tracking Error - %s [%s]', upper(Modes_dac{m}), form_tag), ...
        'FontSize',12,'FontName','Times New Roman');
    legend('Location','northeast','FontSize',9,'NumColumns',2);
    xlim([t_set(1), t_set(end)]);

    fname_fig = sprintf('Comparison_%s_%s', upper(Modes_dac{m}), form_tag);
    saveas(fig, fullfile(SaveFolder_Figures, [fname_fig, '.png']));
    savefig(fig, fullfile(SaveFolder_Figures, [fname_fig, '.fig']));
end

fprintf('完成。图已保存至 %s\n', SaveFolder_Figures);

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