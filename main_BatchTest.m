%% main_BatchTest
clc; close all;

%% 配置 - 在这里修改参数
CurrentMode = 'all';        % 'poe'|'gpoe'|'moe'|'bcm'|'rbcm'|'local'|'exact'
                             % 'poe_ac'|'gpoe_ac'|'moe_ac'|'bcm_ac'|'rbcm_ac'
                             % 'all' 跑所有模式
TestType    = 'both';        % 'test'|'inducing'|'both'

%% 路径设置
SaveFolder_Test     = fullfile('Result', 'Test_Point');
SaveFolder_Inducing = fullfile('Result', 'Inducing_Point');
if ~exist(SaveFolder_Test,    'dir'), mkdir(SaveFolder_Test);    end
if ~exist(SaveFolder_Inducing,'dir'), mkdir(SaveFolder_Inducing); end

%% 模式列表
Modes_dac      = {'poe','gpoe','moe','bcm','rbcm'};
Modes_ac       = {'poe_ac','gpoe_ac','moe_ac','bcm_ac','rbcm_ac'};
Modes_baseline = {'local','exact'};

if strcmpi(CurrentMode, 'all')
    AllModes_test     = [Modes_dac, Modes_baseline];
    AllModes_inducing = [Modes_dac, Modes_ac, Modes_baseline];
else
    AllModes_test     = {CurrentMode};
    AllModes_inducing = {CurrentMode};
end

%% 1. 运行仿真
if ismember(lower(TestType), {'inducing','both'})
    fprintf('\n======== 诱导点聚合 ========\n');
    for m = 1:numel(AllModes_inducing)
        fprintf('[%d/%d] %s\n', m, numel(AllModes_inducing), AllModes_inducing{m});
        run_simulation_inducing_point(AllModes_inducing{m}, SaveFolder_Inducing, AllModes_inducing{m});
    end
end

if ismember(lower(TestType), {'test','both'})
    fprintf('\n======== 测试点聚合 ========\n');
    for m = 1:numel(AllModes_test)
        fprintf('[%d/%d] %s\n', m, numel(AllModes_test), AllModes_test{m});
        run_simulation_test_point(AllModes_test{m}, SaveFolder_Test, AllModes_test{m});
    end
end

%% 2. 加载结果
has_inducing = ~isempty(AllModes_inducing) && ...
    exist(fullfile(SaveFolder_Inducing,[AllModes_inducing{1},'.mat']),'file');
has_test = ~isempty(AllModes_test) && ...
    exist(fullfile(SaveFolder_Test,[AllModes_test{1},'.mat']),'file');

if has_inducing
    temp = load(fullfile(SaveFolder_Inducing,[AllModes_inducing{1},'.mat']),'t_set');
    t_set = temp.t_set;
    N_ind = numel(AllModes_inducing);
    Err_Inducing = zeros(N_ind, numel(t_set));
    for m = 1:N_ind
        d = load(fullfile(SaveFolder_Inducing,[AllModes_inducing{m},'.mat']),'TrackingError_vector');
        len = min(numel(d.TrackingError_vector), numel(t_set));
        Err_Inducing(m,1:len) = d.TrackingError_vector(1:len);
    end
end

if has_test
    temp = load(fullfile(SaveFolder_Test,[AllModes_test{1},'.mat']),'t_set');
    t_set = temp.t_set;
    N_test = numel(AllModes_test);
    Err_Test = zeros(N_test, numel(t_set));
    for m = 1:N_test
        d = load(fullfile(SaveFolder_Test,[AllModes_test{m},'.mat']),'TrackingError_vector');
        len = min(numel(d.TrackingError_vector), numel(t_set));
        Err_Test(m,1:len) = d.TrackingError_vector(1:len);
    end
end

%% 3. 打印表格
fprintf('\n%s\n', repmat('=',1,44));
if has_inducing
    fprintf('Inducing-Point:\n');
    fprintf('  %-12s  %16s\n','Method','Final ||e||');
    fprintf('  %-12s  %16s\n','------','-----------');
    for m = 1:N_ind
        fprintf('  %-12s  %16.4f\n', AllModes_inducing{m}, Err_Inducing(m,end));
    end
    fprintf('\n');
end
if has_test
    fprintf('Test-Point:\n');
    fprintf('  %-12s  %16s\n','Method','Final ||e||');
    fprintf('  %-12s  %16s\n','------','-----------');
    for m = 1:N_test
        fprintf('  %-12s  %16.4f\n', AllModes_test{m}, Err_Test(m,end));
    end
end
fprintf('%s\n\n', repmat('=',1,44));

%% 4. 绘图
if ~exist('t_set','var'), return; end
t_start = t_set(1); t_end = t_set(end);
styles = {'b-','r-','g-','m-','y-','b--','r--','g--','m--','y--','k-','k:','k-.','k--'};

if has_inducing && numel(AllModes_inducing) > 1
    % 找出DAC、AC、baseline的索引
    idx_dac  = find(ismember(AllModes_inducing, Modes_dac));
    idx_ac   = find(ismember(AllModes_inducing, Modes_ac));
    idx_base = find(ismember(AllModes_inducing, Modes_baseline));

    % Figure 1: 诱导点 DAC
    if ~isempty(idx_dac)
        figure('Name','Inducing-Point DAC','Color','w');
        hold on; grid on; box on;
        set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
        for m = idx_dac
            si = mod(m-1,numel(styles))+1;
            plot(t_set, Err_Inducing(m,:), styles{si}, 'LineWidth',1.5, ...
                'DisplayName', AllModes_inducing{m});
        end
        for m = idx_base
            plot(t_set, Err_Inducing(m,:), 'k-', 'LineWidth',1.5, ...
                'DisplayName', AllModes_inducing{m});
            if m ~= idx_base(end)
                % 用不同黑色线型区分local和exact
            end
        end
        ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
        xlabel('$t$','Interpreter','latex','FontSize',13);
        title('Inducing-Point Aggregation (DAC)','FontSize',12,'FontName','Times New Roman');
        legend('Location','northeast','FontSize',9,'NumColumns',2);
        xlim([t_start, t_end]);
        saveas(gcf, fullfile('Result','Inducing_DAC_Comparison.png'));
    end

    % Figure 2: 诱导点 AC
    if ~isempty(idx_ac)
        figure('Name','Inducing-Point AC','Color','w');
        hold on; grid on; box on;
        set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
        styles_ac = {'b--','r--','g--','m--','y--'};
        for i = 1:numel(idx_ac)
            m = idx_ac(i);
            si = mod(i-1,numel(styles_ac))+1;
            plot(t_set, Err_Inducing(m,:), styles_ac{si}, 'LineWidth',1.5, ...
                'DisplayName', AllModes_inducing{m});
        end
        styles_base = {'k-','k:'};
        for i = 1:numel(idx_base)
            m = idx_base(i);
            plot(t_set, Err_Inducing(m,:), styles_base{i}, 'LineWidth',1.5, ...
                'DisplayName', AllModes_inducing{m});
        end
        ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
        xlabel('$t$','Interpreter','latex','FontSize',13);
        title('Inducing-Point Aggregation (AC)','FontSize',12,'FontName','Times New Roman');
        legend('Location','northeast','FontSize',9,'NumColumns',2);
        xlim([t_start, t_end]);
        saveas(gcf, fullfile('Result','Inducing_AC_Comparison.png'));
    end
end

if has_test && numel(AllModes_test) > 1
    figure('Name','Test-Point','Color','w');
    hold on; grid on; box on;
    set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
    for m = 1:N_test
        si = mod(m-1,numel(styles))+1;
        plot(t_set, Err_Test(m,:), styles{si}, 'LineWidth',1.5, ...
            'DisplayName', AllModes_test{m});
    end
    ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
    xlabel('$t$','Interpreter','latex','FontSize',13);
    title('Test-Point Aggregation (DAC)','FontSize',12,'FontName','Times New Roman');
    legend('Location','northeast','FontSize',9,'NumColumns',2);
    xlim([t_start, t_end]);
    saveas(gcf, fullfile('Result','Test_Point_Comparison.png'));
end

fprintf('完成。\n');