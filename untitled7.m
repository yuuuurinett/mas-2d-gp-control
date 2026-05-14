%% test_dac_quick.m
% 快速验证 DAC 向量化版本是否正确
% 预期结果（对比之前 MC=10 的均值）：
%   KIN40K  DAC-MOE  SMSE ≈ 0.11
%   PUMADYN DAC-MOE  SMSE ≈ 0.048
% 只跑 seed=1，约 30 秒

clc;
fprintf('=== 快速验证：DAC 向量化版本 ===\n\n');

%% 测试1：KIN40K（主要验证基本功能）
fprintf('--- KIN40K ---\n');
run_dac_dataset('KIN40K', 'moe', 0.4, 1);

%% 测试2：PUMADYN（验证高维数据）
fprintf('--- PUMADYN32NM ---\n');
run_dac_dataset('PUMADYN32NM', 'moe', 0.4, 1);

fprintf('\n=== 验证完成 ===\n');
fprintf('如果 SMSE 与历史结果相近，说明向量化版本正确。\n');