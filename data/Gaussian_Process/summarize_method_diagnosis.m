function T = summarize_method_diagnosis(DatasetName, M, seed, train_ratio)

if nargin < 1, DatasetName = 'POL'; end
if nargin < 2, M = 2500; end
if nargin < 3, seed = 1; end
if nargin < 4, train_ratio = 0.4; end

Method_list = {'poe','gpoe','moe','bcm','rbcm'};
tr_tag = round(train_ratio * 100);

MainFuncPath = which('run_inducingpoint_dataset_trade_off');
ProjectRoot = fileparts(MainFuncPath);
ResultFolder = fullfile(ProjectRoot, 'Result', 'Dataset', DatasetName);

fprintf('[Summary read folder]\n%s\n', ResultFolder);

Rows = {};

for mm = 1:numel(Method_list)
    method = Method_list{mm};
    fname = sprintf('%s_M%d_tr%d_mc%d.mat', method, M, tr_tag, seed);
    fpath = fullfile(ResultFolder, fname);

    if ~isfile(fpath)
        warning('Missing file: %s', fpath);
        continue;
    end

    S = load(fpath);

    Rows(end+1,:) = { ...
        DatasetName, method, M, seed, ...
        getfield_safe(S,'smse'), ...
        getfield_safe(S,'rmse'), ...
        getfield_safe(S,'nlpd'), ...
        getfield_safe(S,'msll'), ...
        getfield_safe(S,'Mean_var_pred'), ...
        getfield_safe(S,'Median_var_pred'), ...
        getfield_safe(S,'Min_var_pred'), ...
        getfield_safe(S,'Max_var_pred'), ...
        getfield_safe(S,'Mean_err2'), ...
        getfield_safe(S,'Median_err2'), ...
        getfield_safe(S,'Mean_err2_over_var'), ...
        getfield_safe(S,'Median_err2_over_var') ...
    };
end

if isempty(Rows)
    warning('No diagnosis files found for Dataset=%s, M=%d, seed=%d.', DatasetName, M, seed);
    T = table();
    return;
end

T = cell2table(Rows, 'VariableNames', { ...
    'Dataset','Method','M','Seed', ...
    'SMSE','RMSE','NLPD','MSLL', ...
    'MeanVarPred','MedianVarPred','MinVarPred','MaxVarPred', ...
    'MeanErr2','MedianErr2','MeanErr2OverVar','MedianErr2OverVar' ...
});

disp(T);

outname = sprintf('diagnosis_%s_M%d_seed%d.csv', DatasetName, M, seed);
writetable(T, fullfile(ResultFolder, outname));

end

function val = getfield_safe(S, fieldname)
if isfield(S, fieldname)
    val = S.(fieldname);
else
    val = NaN;
end
end