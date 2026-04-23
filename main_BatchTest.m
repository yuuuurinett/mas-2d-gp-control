clc; clear; close all;

%% setting
Modes_5        = {'poe','gpoe','moe','bcm','rbcm'};
Modes_baseline = {'local','exact'};
AllModes       = [Modes_5, Modes_baseline];

SaveFolder_Test     = fullfile('Result','Test_Point');
SaveFolder_Inducing = fullfile('Result','Inducing_Point');

%% 1. simulation
do_simulation = true;
if do_simulation
   
    fprintf('\n======== 测试点聚合 ========\n');
    for m = 1:numel(AllModes)
        fprintf('[%d/%d] %s\n', m, numel(AllModes), AllModes{m});
        run_simulation_test_point(AllModes{m}, SaveFolder_Test, AllModes{m});
    end

    

    fprintf('\n======== 诱导点聚合 ========\n');
    for m = 1:numel(AllModes)
        fprintf('[%d/%d] %s\n', m, numel(AllModes), AllModes{m});
        run_simulation_inducing_point(AllModes{m}, SaveFolder_Inducing, AllModes{m});
    end
    
end

%% 2. 
N = numel(AllModes);
temp = load(fullfile(SaveFolder_Test,[AllModes{1},'.mat']), ...
    't_set','bound_local','bound_distributed','bound_exact');
temp = load(fullfile(SaveFolder_Test,[AllModes{1},'.mat']), 't_set');
t_set             = temp.t_set;
%bound_local       = temp.bound_local;
%bound_distributed = temp.bound_distributed;
%bound_exact       = temp.bound_exact;

%Err_Test     = zeros(N, numel(t_set));
Err_Inducing = zeros(N, numel(t_set));
for m = 1:N
    d1 = load(fullfile(SaveFolder_Test,     [AllModes{m},'.mat']),'TrackingError_vector');
    d2 = load(fullfile(SaveFolder_Inducing, [AllModes{m},'.mat']),'TrackingError_vector');
    Err_Test(m,:)     = d1.TrackingError_vector;
    Err_Inducing(m,:) = d2.TrackingError_vector;
end

%% 3. plot
fprintf('\n%s\n', repmat('=',1,52));
fprintf('  %-10s  %16s  %16s\n','Method','Test Final||e||','Inducing Final||e||');
fprintf('  %-10s  %16s  %16s\n','------','---------------','-------------------');
for m = 1:numel(Modes_5)
    fprintf('  %-10s  %16.4f  %16.4f\n', AllModes{m}, ...      
    Err_Inducing(m,end));
    %Err_Test(m,end));
end
fprintf('  %s\n', repmat('-',1,48));
for m = numel(Modes_5)+1:N
     fprintf('  %-10s  %16.4f  %16.4f\n', AllModes{m}, ...      
      Err_Inducing(m,end));
        %Err_Test(m,end), Err_Inducing(m,end));
end
fprintf('%s\n\n', repmat('=',1,52));

%{
%% 4. 图1: 测试点模式对比
colors = {[0.0,0.45,0.74],[0.85,0.33,0.10],[0.47,0.67,0.19],...
          [0.63,0.08,0.18],[0.93,0.69,0.13],[0.5,0.5,0.5],[0,0,0]};
styles  = {'-','-','-','-','-','--',':'};
lw      = [1.8,1.8,1.8,1.8,1.8,1.5,1.5];
LegendNames = {'PoE','gPoE','MoE','BCM','rBCM','Local','Exact'};
t_start = t_set(1); t_end = t_set(end);

figure('Color','w','Position',[50 100 700 520]);
%subplot(2,1,1); 
hold on; grid on; box on;
set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
for m = 1:N
    plot(t_set,Err_Test(m,:),'Color',colors{m},'LineStyle',styles{m},'LineWidth',lw(m));
end
ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
title('Test-Point Aggregation (DAC)','FontSize',12,'FontName','Times New Roman');
legend(LegendNames,'Location','northeast','FontSize',9,'NumColumns',2);
xlim([t_start,t_end]); set(gca,'XTickLabel',[]);

%subplot(2,1,2); hold on; grid on; box on;
set(gca,'FontSize',11,'FontName','Times New Roman');
%plot(t_set,bound_local,'k-','LineWidth',1.5,'DisplayName','v(t) local');
%plot(t_set,bound_distributed,'r-','LineWidth',1.5,'DisplayName','v(t) distributed');
%plot(t_set,bound_exact,'b--','LineWidth',1.5,'DisplayName','v(t) exact');
%ylabel('$v(t)$','Interpreter','latex','FontSize',13);
%xlabel('$t$','Interpreter','latex','FontSize',13);
xlim([t_start,t_end]); legend('Location','northeast','FontSize',10);
saveas(gcf,fullfile(SaveFolder_Test,'Test_Point_Comparison.fig'));
saveas(gcf,fullfile(SaveFolder_Test,'Test_Point_Comparison.png'));
%}


%% 5. 图2: 诱导点模式对比
figure('Color','w','Position',[800 100 700 520]);
hold on; grid on; box on;
set(gca,'YScale','log','FontSize',11,'FontName','Times New Roman');
for m = 1:N
    plot(t_set,Err_Inducing(m,:),'Color',colors{m},'LineStyle',styles{m},'LineWidth',lw(m));
end
ylabel('$\|e\|$','Interpreter','latex','FontSize',13);
title('Inducing-Point Aggregation (DAC)','FontSize',12,'FontName','Times New Roman');
legend(LegendNames,'Location','northeast','FontSize',9,'NumColumns',2);
xlim([t_start,t_end]); set(gca,'XTickLabel',[]);
saveas(gcf,fullfile(SaveFolder_Inducing,'Inducing_Point_Comparison.fig'));
fprintf('所有图已保存。\n');