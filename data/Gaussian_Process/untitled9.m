%% Cost-aware fine selection of M
% This block selects M after the coarse Low/Medium/High screening.
% Criterion:
%   score = normalized SMSE + normalized MSLL + lambda * normalized train time
% Lower score is better.

lambda_time = 0.5;

FineRegion.KIN40K      = [1000 2000];
FineRegion.POL         = [1000 2000];
FineRegion.PUMADYN32NM = [100 1000];
FineRegion.SARCOS      = [2000 2500];

fprintf('\n============================================================\n');
fprintf('Cost-aware fine selection of M\n');
fprintf('score = norm(SMSE) + norm(MSLL) + %.2f * norm(Time)\n', lambda_time);
fprintf('============================================================\n');

for di = 1:numel(Dataset_list)
    dataset = Dataset_list{di};
    range_M = FineRegion.(dataset);

    idx_region = M_list >= range_M(1) & M_list <= range_M(2);
    M_region = M_list(idx_region);

    % Average over methods to select dataset-level M.
    smse_mat = [];
    msll_mat = [];
    time_mat = [];

    for ci = 1:numel(Method_list)
        method = Method_list{ci};
        d = data.(method).(dataset);

        smse_mat(:,ci) = d.SMSE_mean(idx_region);
        msll_mat(:,ci) = d.MSLL_mean(idx_region);
        time_mat(:,ci) = d.TrainTime_mean(idx_region);
    end

    avg_smse = mean(smse_mat, 2, 'omitnan');
    avg_msll = mean(msll_mat, 2, 'omitnan');
    avg_time = mean(time_mat, 2, 'omitnan');

    norm_smse = normalize_minmax_cost(avg_smse);
    norm_msll = normalize_minmax_cost(avg_msll);
    norm_time = normalize_minmax_cost(avg_time);

    score = norm_smse + norm_msll + lambda_time * norm_time;

    [best_score, best_idx] = min(score);
    selected_M = M_region(best_idx);

    fprintf('\nDataset: %s\n', dataset);
    fprintf('Fine region: M = %d to %d\n', range_M(1), range_M(2));
    fprintf('Selected M = %d, score = %.4f\n', selected_M, best_score);

    T_select = table(M_region(:), avg_smse(:), avg_msll(:), avg_time(:), score(:), ...
        'VariableNames', {'M','AvgSMSE','AvgMSLL','AvgTrainTime_mspt','Score'});

    disp(T_select);
end

function y_norm = normalize_minmax_cost(y)
    y = y(:);
    y_min = min(y, [], 'omitnan');
    y_max = max(y, [], 'omitnan');

    if abs(y_max - y_min) < 1e-12
        y_norm = zeros(size(y));
    else
        y_norm = (y - y_min) ./ (y_max - y_min);
    end
end