%% Quick method diagnosis for MOE abnormal MSLL behavior
% Put this file in mas_2D_test/data/Gaussian_Process/ and run it there.
% It compares POE/GPOE/MOE/BCM/RBCM at selected M values.

clear; clc;

Dataset_list = {'KIN40K','POL'};
Method_list  = {'poe','gpoe','moe','bcm','rbcm'};
M_list       = [1000 1500 2000 2500];
seed         = 1;
train_ratio  = 0.4;

for dd = 1:numel(Dataset_list)
    DatasetName = Dataset_list{dd};

    for ii = 1:numel(M_list)
        M = M_list(ii);

        fprintf('\n============================================================\n');
        fprintf('Diagnosis dataset=%s | M=%d | seed=%d\n', DatasetName, M, seed);
        fprintf('============================================================\n');

        for mm = 1:numel(Method_list)
            CurrentMode = Method_list{mm};

            fprintf('\n------------------------------\n');
            fprintf('Run method=%s | M=%d\n', CurrentMode, M);
            fprintf('------------------------------\n');

            run_inducingpoint_dataset_trade_off(DatasetName, CurrentMode, train_ratio, seed, M);
        end
    end
end
