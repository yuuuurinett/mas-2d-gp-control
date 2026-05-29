function [mu_batch, var_batch] = batch_predict_external(gp, X_query, SigmaN, SigmaF)
% 外部批量预测函数，绕过 LocalGP_MultiOutput.predict() 的单点限制
%
% 数学上完全等价于逐点调用 gp.predict(x)，但对所有查询点一次性完成矩阵运算
% 不修改导师代码，仅通过 gp 对象的公开属性和 kernel 方法实现
%
% 输入:
%   gp      : LocalGP_MultiOutput 对象（已调用 add_Alldata 完成训练）
%   X_query : [x_dim × N_query]  查询点矩阵（列为点）
%   SigmaN  : 噪声标准差（标量，用于方差下界截断）
%   SigmaF  : 信号标准差（标量，用于先验方差）
%
% 输出:
%   mu_batch  : [N_query × y_dim]  预测均值
%   var_batch : [N_query × 1]      预测方差（各输出维度相同，与单点版一致）
%
% 用法示例（替换预计算循环中对 predict 的调用）:
%   [mn, vn] = batch_predict_external(LocalGP_set{n}, X_eval', SigmaN, SigmaF);
%   % mn: [N_eval × y_dim],  vn: [N_eval × 1]

N_data  = gp.DataQuantity;
alpha   = gp.alpha(1:N_data, :);          % [N_data × y_dim]
Chol_L  = gp.L(1:N_data, 1:N_data);      % [N_data × N_data] 下三角
X_train = gp.X(:, 1:N_data);             % [x_dim × N_data]

% K_star: [N_data × N_query]
% gp.kernel 的调用方式与内部一致: kernel(X1, X2) → [size(X1,2) × size(X2,2)]
K_star = gp.kernel(X_train, X_query);    % [N_data × N_query]

% 预测均值: mu = K_star' * alpha
mu_batch = (alpha' * K_star)';           % [N_query × y_dim]

% 预测方差: var = sigma_f^2 - diag(V'V), V = L\K_star
% Chol_L 是下三角，所以 L\K_star 用 mldivide 即可
V        = Chol_L \ K_star;              % [N_data × N_query]
var_vec  = max(SigmaF^2 - sum(V.^2, 1)', SigmaN^2);  % [N_query × 1]

% 与单点版一致：方差对所有 y_dim 相同（各输出共享同一核函数）
var_batch = var_vec;                     % [N_query × 1]

end