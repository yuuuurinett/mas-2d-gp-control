function [MeanTable,StdTable,LongTable] = ...
    run_mc50_all_aggregation_methods( ...
    do_run,seed_subset,do_summarize)
%RUN_MC50_ALL_AGGREGATION_METHODS Final 50-seed control comparison.
%
% Default configuration:
%   seeds                 = 1:50  (全新运行，不保留/复用之前的)
%   simulation horizon    = 30 s
%   final IP setting      = M=600
%   common TP/CEN/NBR/Offline/Exact baseline setting = M=400 folder

% 💡【新增】：启动计时器
tStart = tic;

% 1. 参数初始化
if nargin < 1 || isempty(do_run), do_run = true; end
if nargin < 2, seed_subset = []; end
if nargin < 3 || isempty(do_summarize), do_summarize = ~do_run; end

mc_count = 50;

% Running all remaining seeds inside one MATLAB process caused substantial
% memory accumulation overnight.  An empty subset now completes only the
% next missing M600 seed and returns.  For an unattended full resume use
% run_mc50_m600_safe.ps1, which launches a fresh MATLAB process per seed.
if do_run && isempty(seed_subset)
    completed_seed = run_mc50_next_ip_seed();
    MeanTable = table(); StdTable = table(); LongTable = table();
    if isempty(completed_seed), return; end
    return;
end
if do_run && numel(seed_subset) ~= 1
    error(['For memory safety, run one seed per MATLAB process. ' ...
        'Use run_mc50_m600_safe.ps1 for a seed range.']);
end

if do_run
    fprintf('\n🚀 开始仿真，当前运行的种子范围: %s\n', mat2str(seed_subset));
end

% 2. M400 common baselines are already complete.  Do not reopen them when
% resuming M600 IP seeds.

% 3. 运行 IP-DAC/IP-AC 并执行最终的 Summarize - M=600
[MeanTable,StdTable,LongTable] = ...
    run_mc10_all_aggregation_methods( ...
    do_run, seed_subset, do_summarize, 600, true, mc_count);

% 💡【新增】：停止计时并打印耗时
elapsedTime = toc(tStart); % 获取总秒数
h = floor(elapsedTime / 3600); % 计算小时
m = floor(mod(elapsedTime, 3600) / 60); % 计算分钟
s = mod(elapsedTime, 60); % 计算剩余秒数

fprintf('\n====================================================\n');
if do_run
    fprintf('✅ 仿真运行完毕！\n');
    fprintf('本次运行种子: %s\n', mat2str(seed_subset));
    fprintf('总共耗时: %d 小时 %d 分钟 %.2f 秒\n', h, m, s);
else
    fprintf('✅ 数据汇总完毕！\n');
    fprintf('汇总耗时: %.2f 秒\n', elapsedTime);
end
fprintf('====================================================\n\n');

end
