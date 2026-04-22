clc; clear; close all;

%% 1. 配置参数

SimulationModes = {'distributed_POE', 'masked_POE', 'distributed_RBCM', 'local', 'none', 'exact'};


SaveFolderName = fullfile('Result', 'Main_Simulation_Results');

%% 2. 
do_simulation = true; 
if do_simulation
    for ModeNr = 1:numel(SimulationModes)
        CurrentMode = SimulationModes{ModeNr};
        SaveFileName = CurrentMode;

        fprintf('开始运行模式 [%d/%d]: %s\n', ModeNr, numel(SimulationModes), CurrentMode);

        run_main_simulation_mode(CurrentMode, SaveFolderName, SaveFileName);
    end
end

%% 3. 
NumModes = numel(SimulationModes);

FirstFile = fullfile(SaveFolderName, [SimulationModes{1}, '.mat']);
%if ~exist(FirstFile, 'file')
   % error('未找到数据文件，请确保 do_simulation = true 运行了至少一次。');
%end

temp_data = load(FirstFile, 't_set', 'bound_local', 'bound_distributed', 'bound_exact');
t_set = temp_data.t_set;
bound_local = temp_data.bound_local;
bound_distributed = temp_data.bound_distributed;
bound_exact = temp_data.bound_exact;


TrackingError_matrix = zeros(NumModes, numel(t_set));

for ModeNr = 1:NumModes
    CurrentMode = SimulationModes{ModeNr};
    DataFile = fullfile(SaveFolderName, [CurrentMode, '.mat']);
    LoadedResult = load(DataFile, 'TrackingError_vector');
    
  
    TrackingError_matrix(ModeNr, :) = LoadedResult.TrackingError_vector;
end

%% 4. Plot
figure('Color','w','Position',[100 100 700 550]);
t_start = t_set(1);
t_end = t_set(end);

% 上图: ||e||

subplot(2,1,1);
styles = {'r-', 'g-', 'm-', 'k-', 'b-', 'b--'};
hold on; grid on; box on;
set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');

for ModeNr = 1:NumModes
    style_idx = mod(ModeNr - 1, numel(styles)) + 1; 
    plot(t_set, TrackingError_matrix(ModeNr,:), styles{style_idx}, 'LineWidth',1.5);
end

ylabel('$\|e\|$','Interpreter','latex','FontSize',12);
xlim([t_start, t_end]);
LegendNames = strrep(SimulationModes, '_', ' '); 
legend(LegendNames, 'Location','northeast','FontSize',10);
set(gca,'XTickLabel',[]);

% 下图: v(t) 误差上界

subplot(2,1,2);
hold on; grid on; box on;
set(gca,'FontSize',11,'FontName','Times New Roman');

plot(t_set, bound_local, 'k-',  'LineWidth',1.5, 'DisplayName','local');
plot(t_set, bound_distributed,  'r-',  'LineWidth',1.5, 'DisplayName','distributed');
plot(t_set, bound_exact, 'b--', 'LineWidth',1.5, 'DisplayName','exact');

ylabel('$v(t)$','Interpreter','latex','FontSize',12);
xlabel('$t$','Interpreter','latex','FontSize',12);
xlim([t_start, t_end]);
legend('Location','northeast','FontSize',10);


%saveas(gcf, fullfile(SaveFolderName, 'Simulation_Plot.png'));

saveas(gcf, fullfile(SaveFolderName, 'Simulation_Plot.fig'));