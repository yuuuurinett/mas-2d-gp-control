tr = load('POL_train.mat');
X = tr.x;
N = size(X, 1);
AgentQuantity = 6;
n_per_agent = floor(N / AgentQuantity);
k = 5;  % 找第k个最近邻

% 先标准化（消除量纲影响）
X_norm = (X - mean(X)) / std(X);

% 全部数据的平均最近邻距离（模拟LoG，数据统一管理）
rng(0);
idx_sample = randperm(N, 500);  % 抽500个点算，省时间
D_full = pdist2(X_norm(idx_sample,:), X_norm, 'euclidean');
D_full_sorted = sort(D_full, 2);
avg_dist_full = mean(D_full_sorted(:, k+1));  % 第1列是自己(=0)

% 每个agent的平均最近邻距离（模拟分布式方法）
avg_dist_agent = zeros(AgentQuantity, 1);
for a = 1:AgentQuantity
    idx_agent = (a-1)*n_per_agent+1 : a*n_per_agent;
    X_agent = X_norm(idx_agent, :);
    n_a = size(X_agent, 1);
    idx_s = randperm(n_a, min(100, n_a));
    D_a = pdist2(X_agent(idx_s,:), X_agent, 'euclidean');
    D_a_sorted = sort(D_a, 2);
    avg_dist_agent(a) = mean(D_a_sorted(:, k+1));
end

fprintf('全部数据的平均第%d近邻距离: %.4f\n', k, avg_dist_full);
fprintf('每个agent的平均第%d近邻距离:\n', k);
for a = 1:AgentQuantity
    fprintf('  Agent %d: %.4f (比全局稀疏 %.1fx)\n', ...
        a, avg_dist_agent(a), avg_dist_agent(a)/avg_dist_full);
end
fprintf('平均稀疏化倍数: %.1fx\n', mean(avg_dist_agent)/avg_dist_full);