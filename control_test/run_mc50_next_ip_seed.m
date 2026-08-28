function completed_seed = run_mc50_next_ip_seed()
%RUN_MC50_NEXT_IP_SEED Resume exactly one incomplete M600 MC50 seed.
% Run this from MATLAB when a short, restartable unit of work is desired.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
result_folder = fullfile(repo_root,'result', ...
    'mc50_all_aggregation_m600_T30_per_agent');
methods = {'poe','gpoe','moe','bcm','rbcm'};

completed_seed = [];
for seed = 1:50
    seed_folder = fullfile(result_folder,sprintf('seed_%03d',seed));
    missing = false;
    for method_i = 1:numel(methods)
        missing = missing || ~isfile(fullfile(seed_folder, ...
            [methods{method_i} '_ip_dac.mat']));
        missing = missing || ~isfile(fullfile(seed_folder, ...
            [methods{method_i} '_ip_ac.mat']));
    end
    if missing
        fprintf('Resuming MC50 M600 seed %d only.\n',seed);
        run_mc10_all_aggregation_methods(true,seed,false,600,true,50);
        completed_seed = seed;
        fprintf('Seed %d complete. Restart MATLAB before the next seed ',seed);
        fprintf('if memory usage has grown.\n');
        return;
    end
end

fprintf(['All 50 M600 IP seeds are complete. Summary was not generated ' ...
    'because the common M400 batch must first pass a configuration audit.\n']);
end
