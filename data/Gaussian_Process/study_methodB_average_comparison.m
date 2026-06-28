%% study_methodB_average_comparison.m
% Method B simplified study:
%   No percentage threshold.
%   No score function.
%
% Idea:
%   1. Use the coarse region selected from the Low/Medium/High table.
%   2. For each M in that region, average SMSE, MSLL, and training time
%      over the five aggregation methods: POE, GPOE, MOE, BCM, RBCM.
%   3. Output one comparison table per dataset.
%   4. Mark:
%        - best averaged SMSE: lowest AvgSMSE
%        - best averaged MSLL: lowest AvgMSLL
%        - fastest: lowest AvgTrainTime
%
% This script does NOT automatically force a single selection rule.
% It gives a clean table for manual trade-off explanation.
%
% Put this file in:
%   data/Gaussian_Process/
%
% Run:
%   study_methodB_average_comparison
%
% Output:
%   Result/Figures/M_ablation_methodB_average/
%       MethodB_average_comparison.md

clear; clc; close all;

%% ===================== User settings =====================

Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
Method_label = {'POE','GPOE','MOE','BCM','RBCM'}; %#ok<NASGU>

Dataset_list = {'KIN40K','POL','PUMADYN32NM','SARCOS'};

M_list      = 100:100:2500;
seeds       = 1:3;
train_ratio = 0.4;
tr_tag      = round(train_ratio * 100);

ProjectRoot = fileparts(mfilename('fullpath'));

% If needed, force it manually:
% ProjectRoot = 'C:\Users\Yurou Du\Desktop\mas-2d-gp-control\data\Gaussian_Process';

ResultRoot = fullfile(ProjectRoot, 'Result', 'Dataset');

OutFolder = fullfile(ProjectRoot, 'Result', 'Figures', 'M_ablation_methodB_average');
if ~exist(OutFolder, 'dir')
    mkdir(OutFolder);
end

%% ===================== Coarse regions =====================
% These regions come from the Low/Medium/High table.

FineRegion.KIN40K      = [1000 2000];  % Medium region
FineRegion.POL         = [1000 2000];  % Medium region
FineRegion.PUMADYN32NM = [100 1000];   % Low region
FineRegion.SARCOS      = [2000 2500];  % High region

CoarseRegionName.KIN40K      = 'Medium';
CoarseRegionName.POL         = 'Medium';
CoarseRegionName.PUMADYN32NM = 'Low';
CoarseRegionName.SARCOS      = 'High';

%% ===================== Read data =====================

data = struct();
MissingCounter = 0;

for ci = 1:numel(Method_list)
    method = Method_list{ci};

    for di = 1:numel(Dataset_list)
        dataset = Dataset_list{di};

        ResultFolder = fullfile(ResultRoot, dataset);

        SMSE_all      = nan(numel(M_list), numel(seeds));
        MSLL_all      = nan(numel(M_list), numel(seeds));
        TrainTime_all = nan(numel(M_list), numel(seeds));   % ms/pt

        for mi = 1:numel(M_list)
            M = M_list(mi);

            for si = 1:numel(seeds)
                seed = seeds(si);

                fname = sprintf('%s_M%d_tr%d_mc%d.mat', method, M, tr_tag, seed);
                fpath = fullfile(ResultFolder, fname);

                if ~exist(fpath, 'file')
                    MissingCounter = MissingCounter + 1;
                    fprintf('[Missing] %s\n', fpath);
                    continue;
                end

                S = load(fpath, 'smse', 'msll', ...
                    't_train_per_point', 't_train_total', 'N_train');

                if isfield(S, 'smse')
                    SMSE_all(mi, si) = S.smse;
                end

                if isfield(S, 'msll')
                    MSLL_all(mi, si) = S.msll;
                end

                if isfield(S, 't_train_per_point')
                    TrainTime_all(mi, si) = S.t_train_per_point;
                elseif isfield(S, 't_train_total') && isfield(S, 'N_train') && S.N_train > 0
                    TrainTime_all(mi, si) = (S.t_train_total / S.N_train) * 1000;
                end
            end
        end

        data.(method).(dataset).SMSE_mean      = mean(SMSE_all, 2, 'omitnan');
        data.(method).(dataset).MSLL_mean      = mean(MSLL_all, 2, 'omitnan');
        data.(method).(dataset).TrainTime_mean = mean(TrainTime_all, 2, 'omitnan');
    end
end

fprintf('\nData loading finished. Missing files: %d\n', MissingCounter);

%% ===================== Average comparison within each coarse region =====================

AllTables = struct();

fprintf('\n============================================================\n');
fprintf('Method B simplified: average comparison without thresholds\n');
fprintf('============================================================\n');

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};

    region = FineRegion.(dataset);
    coarse_name = CoarseRegionName.(dataset);

    idx_region = M_list >= region(1) & M_list <= region(2);
    M_region = M_list(idx_region);

    smse_mat = nan(numel(M_region), numel(Method_list));
    msll_mat = nan(numel(M_region), numel(Method_list));
    time_mat = nan(numel(M_region), numel(Method_list));

    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);

        smse_mat(:,ci) = d.SMSE_mean(idx_region);
        msll_mat(:,ci) = d.MSLL_mean(idx_region);
        time_mat(:,ci) = d.TrainTime_mean(idx_region);
    end

    AvgSMSE = mean(smse_mat, 2, 'omitnan');
    AvgMSLL = mean(msll_mat, 2, 'omitnan');
    AvgTime = mean(time_mat, 2, 'omitnan');

    [~, idx_best_smse] = min(AvgSMSE);
    [~, idx_best_msll] = min(AvgMSLL);
    [~, idx_fastest]   = min(AvgTime);

    BestSMSEFlag = repmat({''}, numel(M_region), 1);
    BestMSLLFlag = repmat({''}, numel(M_region), 1);
    FastestFlag  = repmat({''}, numel(M_region), 1);

    BestSMSEFlag{idx_best_smse} = 'best SMSE';
    BestMSLLFlag{idx_best_msll} = 'best MSLL';
    FastestFlag{idx_fastest}    = 'fastest';

    Recommendation = make_recommendation(dataset, M_region, AvgSMSE, AvgMSLL, AvgTime, ...
        idx_best_smse, idx_best_msll, idx_fastest);

    T = table(M_region(:), AvgSMSE(:), AvgMSLL(:), AvgTime(:), ...
        BestSMSEFlag, BestMSLLFlag, FastestFlag, ...
        'VariableNames', {'M','AvgSMSE','AvgMSLL','AvgTrainTime_mspt', ...
                          'BestSMSE','BestMSLL','Fastest'});

    AllTables.(dataset).Table = T;
    AllTables.(dataset).CoarseRegion = coarse_name;
    AllTables.(dataset).RegionStart = region(1);
    AllTables.(dataset).RegionEnd = region(2);
    AllTables.(dataset).Recommendation = Recommendation;

    fprintf('\nDataset: %s\n', dataset);
    fprintf('Coarse region: %s, fine range: M = %d to %d\n', ...
        coarse_name, region(1), region(2));
    disp(T);
    fprintf('Recommendation: %s\n', Recommendation);
end

%% ===================== Save Markdown report =====================

md_path = fullfile(OutFolder, 'MethodB_average_comparison.md');
fid = fopen(md_path, 'w');

fprintf(fid, '# Method B: Average Comparison without Thresholds\n\n');

fprintf(fid, 'This report uses a simplified Method B. No percentage threshold and no score function are used.\n\n');

fprintf(fid, 'For each dataset, the coarse region is first selected from the Low/Medium/High table. Within that region, SMSE, MSLL, and training time are averaged over the five aggregation methods: POE, GPOE, MOE, BCM, and RBCM.\n\n');

fprintf(fid, 'Lower SMSE and lower MSLL indicate better predictive performance. Lower training time indicates lower computational cost.\n\n');

fprintf(fid, '## Summary\n\n');
fprintf(fid, '| Dataset | Coarse Region | Fine Range | Best Avg SMSE | Best Avg MSLL | Fastest | Recommendation |\n');
fprintf(fid, '|---|---|---:|---:|---:|---:|---|\n');

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};
    T = AllTables.(dataset).Table;

    idx_best_smse = find(strcmp(T.BestSMSE, 'best SMSE'), 1);
    idx_best_msll = find(strcmp(T.BestMSLL, 'best MSLL'), 1);
    idx_fastest   = find(strcmp(T.Fastest, 'fastest'), 1);

    range_text = sprintf('%d-%d', ...
        AllTables.(dataset).RegionStart, ...
        AllTables.(dataset).RegionEnd);

    fprintf(fid, '| %s | %s | %s | M=%d | M=%d | M=%d | %s |\n', ...
        dataset, ...
        AllTables.(dataset).CoarseRegion, ...
        range_text, ...
        T.M(idx_best_smse), ...
        T.M(idx_best_msll), ...
        T.M(idx_fastest), ...
        AllTables.(dataset).Recommendation);
end

fprintf(fid, '\n\n## Detailed Average Tables\n\n');

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};
    T = AllTables.(dataset).Table;

    fprintf(fid, '\n\n### %s\n\n', dataset);
    fprintf(fid, 'Coarse region: **%s**, fine range: **M = %d-%d**.\n\n', ...
        AllTables.(dataset).CoarseRegion, ...
        AllTables.(dataset).RegionStart, ...
        AllTables.(dataset).RegionEnd);

    fprintf(fid, '| M | Avg SMSE | Avg MSLL | Avg Train Time (ms/pt) | Best SMSE | Best MSLL | Fastest |\n');
    fprintf(fid, '|---:|---:|---:|---:|---|---|---|\n');

    for r = 1:height(T)
        fprintf(fid, '| %d | %.5g | %.5g | %.5g | %s | %s | %s |\n', ...
            T.M(r), ...
            T.AvgSMSE(r), ...
            T.AvgMSLL(r), ...
            T.AvgTrainTime_mspt(r), ...
            T.BestSMSE{r}, ...
            T.BestMSLL{r}, ...
            T.Fastest{r});
    end

    fprintf(fid, '\nRecommendation: %s\n', AllTables.(dataset).Recommendation);
end

fclose(fid);

fprintf('\nMarkdown report saved to:\n%s\n', md_path);

% Optional: automatically open the Markdown file after running.
% open(md_path);

%% ========================================================================
%% Local helper function
%% ========================================================================

function recommendation = make_recommendation(dataset, M_region, AvgSMSE, AvgMSLL, AvgTime, ...
    idx_best_smse, idx_best_msll, idx_fastest)

    M_best_smse = M_region(idx_best_smse);
    M_best_msll = M_region(idx_best_msll);
    M_fastest   = M_region(idx_fastest);

    switch dataset
        case 'KIN40K'
            recommendation = sprintf(['Choose a medium-region M near the best MSLL instead of the largest M. ', ...
                'Best Avg SMSE is at M=%d, best Avg MSLL is at M=%d, and the fastest point is M=%d. ', ...
                'The final choice should prioritize the MSLL-stable medium region while keeping SMSE low.'], ...
                M_best_smse, M_best_msll, M_fastest);

        case 'POL'
            recommendation = sprintf(['Choose a medium-region M near the best MSLL/SMSE compromise. ', ...
                'Best Avg SMSE is at M=%d, best Avg MSLL is at M=%d, and the fastest point is M=%d. ', ...
                'Avoid choosing a larger M if it mainly improves SMSE but worsens MSLL or costs more.'], ...
                M_best_smse, M_best_msll, M_fastest);

        case 'PUMADYN32NM'
            recommendation = sprintf(['Choose the low-region efficient point. ', ...
                'Best Avg SMSE is at M=%d, best Avg MSLL is at M=%d, and the fastest point is M=%d. ', ...
                'Since the metrics are nearly saturated, the smallest/fastest M is preferred.'], ...
                M_best_smse, M_best_msll, M_fastest);

        case 'SARCOS'
            recommendation = sprintf(['Choose the high-region point with the best predictive performance. ', ...
                'Best Avg SMSE is at M=%d, best Avg MSLL is at M=%d, and the fastest point is M=%d. ', ...
                'Since performance continues to improve at high M, a larger M is justified.'], ...
                M_best_smse, M_best_msll, M_fastest);

        otherwise
            recommendation = sprintf(['Best Avg SMSE is at M=%d, best Avg MSLL is at M=%d, ', ...
                'and the fastest point is M=%d.'], ...
                M_best_smse, M_best_msll, M_fastest);
    end

    % Keep variables referenced for clarity and future extension.
    %#ok<NASGU>
    AvgSMSE = AvgSMSE;
    AvgMSLL = AvgMSLL;
    AvgTime = AvgTime;
end
