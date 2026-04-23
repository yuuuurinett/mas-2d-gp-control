function main_BatchTest(mode, framework)
%% main_BatchTest
%  用法:
%   main_BatchTest()                     % 跑所有模式，两个框架
%   main_BatchTest('rbcm', 'inducing')   % 只跑rbcm，诱导点
%   main_BatchTest('bcm',  'test')       % 只跑bcm，测试点
%   main_BatchTest('poe',  'both')       % 只跑poe，两个框架

%% 配置
Modes_5        = {'poe','gpoe','moe','bcm','rbcm'};
Modes_baseline = {'local','exact'};

SaveFolder_Test     = fullfile('Result','Test_Point');
SaveFolder_Inducing = fullfile('Result','Inducing_Point');

line_styles = {'b-','r-','g-','m-','y-','k--','k:'};

if nargin == 0
    AllModes  = [Modes_5, Modes_baseline];
    framework = 'both';
else
    AllModes  = {mode};
end
framework = lower(framework);

%% 1. 运行仿真
do_simulation = true;
if do_simulation
    if ismember(framework, {'inducing','both'})
        fprintf('\n======== 诱导点聚合 ========\n');
        for m = 1:numel(AllModes)
            fprintf('[%d/%d] %s\n', m, numel(AllModes), AllModes{m});
            run_simulation_inducing_point(AllModes{m}, SaveFolder_Inducing, AllModes{m});
        end
    end

    if ismember(framework, {'test','both'})
        fprintf('\n======== 测试点聚合 ========\n');
        for m = 1:numel(AllModes)
            fprintf('[%d/%d] %s\n', m, numel(AllModes), AllModes{m});
            run_simulation_test_point(AllModes{m}, SaveFolder_Test, AllModes{m});
        end
    end
end

%% 2. 加载结果
N = numel(AllModes);

has_inducing = exist(fullfile(SaveFolder_Inducing, [AllModes{1},'.mat']), 'file');
has_test     = exist(fullfile(SaveFolder_Test,     [AllModes{1},'.mat']), 'file');

if has_inducing
    temp = load(fullfile(SaveFolder_Inducing,[AllModes{1},'.mat']), 't_set');
    t_set = temp.t_set;
    Err_Inducing = zeros(N, numel(t_set));
    for m = 1:N
        d = load(fullfile(SaveFolder_Inducing,[AllModes{m},'.mat']),'TrackingError_vector');
        Err_Inducing(m,:) = d.TrackingError_vector;
    end
elseif has_test
    temp = load(fullfile(SaveFolder_Test,[AllModes{1},'.mat']), 't_set');
    t_set = temp.t_set;
end

if has_test
    Err_Test = zeros(N, numel(t_set));
    for m = 1:N
        d = load(fullfile(SaveFolder_Test,[AllModes{m},'.mat']),'TrackingError_vector');
        Err_Test(m,:) = d.TrackingError_vector;
    end
end

%% 3. 打印表格
fprintf('\n%s\n', repmat('=',1,44));
fprintf('  %-10s', 'Method');
if has_inducing, fprintf('  %16s', 'Inducing ||e||'); end
if has_test,     fprintf('  %16s', 'Test ||e||');     end
fprintf('\n  %-10s', '------');
if has_inducing, fprintf('  %16s', '---------------'); end
if has_test,     fprintf('  %16s', '---------------'); end
fprintf('\n');
for m = 1:N
    fprintf('  %-10s', AllModes{m});
    if has_inducing, fprintf('  %16.4f', Err_Inducing(m,end)); end
    if has_test,     fprintf('  %16.4f', Err_Test(m,end));     end
    fprintf('\n');
end
fprintf('%s\n\n', repmat('=',1,44));

%% 4. 绘图
t_start = t_set(1); t_end = t_set(end);

if has_inducing && has_test
    % 两张图并排对比
    figure('Color','w','Position',[50 100 1400 450]);

    subplot(1,2,1);
    hold on; grid on; box on;
    set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
    for m = 1:N
        plot(t_set, Err_Inducing(m,:), line_styles{m}, 'LineWidth', 1.5);
    end
    ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
    xlabel('$t$','Interpreter','latex','FontSize',13);
    title('Inducing-Point Aggregation (DAC)','FontSize',12,'FontName','Times New Roman');
    legend(AllModes,'Location','northeast','FontSize',10);
    xlim([t_start, t_end]);

    subplot(1,2,2);
    hold on; grid on; box on;
    set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
    for m = 1:N
        plot(t_set, Err_Test(m,:), line_styles{m}, 'LineWidth', 1.5);
    end
    ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
    xlabel('$t$','Interpreter','latex','FontSize',13);
    title('Test-Point Aggregation (DAC)','FontSize',12,'FontName','Times New Roman');
    legend(AllModes,'Location','northeast','FontSize',10);
    xlim([t_start, t_end]);

    saveas(gcf, fullfile(SaveFolder_Inducing,'Comparison_Both.fig'));
    saveas(gcf, fullfile(SaveFolder_Inducing,'Comparison_Both.png'));

elseif has_inducing
    figure('Color','w');
    hold on; grid on; box on;
    set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
    for m = 1:N
        plot(t_set, Err_Inducing(m,:), line_styles{m}, 'LineWidth', 1.5);
    end
    ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
    xlabel('$t$','Interpreter','latex','FontSize',13);
    title('Inducing-Point Aggregation (DAC)','FontSize',12,'FontName','Times New Roman');
    legend(AllModes,'Location','northeast','FontSize',10);
    xlim([t_start, t_end]);
    saveas(gcf, fullfile(SaveFolder_Inducing,'Inducing_Point_Comparison.fig'));
    saveas(gcf, fullfile(SaveFolder_Inducing,'Inducing_Point_Comparison.png'));

elseif has_test
    figure('Color','w');
    hold on; grid on; box on;
    set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
    for m = 1:N
        plot(t_set, Err_Test(m,:), line_styles{m}, 'LineWidth', 1.5);
    end
    ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
    xlabel('$t$','Interpreter','latex','FontSize',13);
    title('Test-Point Aggregation (DAC)','FontSize',12,'FontName','Times New Roman');
    legend(AllModes,'Location','northeast','FontSize',10);
    xlim([t_start, t_end]);
    saveas(gcf, fullfile('Result','Test_Point_Comparison.fig'));
    saveas(gcf, fullfile('Result','Test_Point_Comparison.png'));
end

fprintf('完成。\n');
end