function MaskedGP = gp_masked_gp_ac( ...
    LocalGP_set, AdjMatrix, NeighbourSet, ...
    AgentQuantity, ...
    NumInducingPoints, InducingPoints_Coordinates, SigmaF, SigmaL, ...
    x_dim, y_dim)

M = NumInducingPoints;
NumNumeratorDenominator = 2 * y_dim; 
ZetaDim = NumNumeratorDenominator * M; 

% 计算 P
P_ReferenceSignal_matrix = zeros(ZetaDim, AgentQuantity);
for AgentNr = 1:AgentQuantity
    for InducingPointIdx = 1:M
        CurrentInducingPoint = InducingPoints_Coordinates(:, InducingPointIdx);
        [mu_n, var_n] = LocalGP_set{AgentNr}.predict(CurrentInducingPoint);
        
        P_ReferenceSignal_matrix(InducingPointIdx, AgentNr) = ...
            AgentQuantity * mu_n(1) / var_n(1);
        P_ReferenceSignal_matrix(M + InducingPointIdx, AgentNr) = ...
            AgentQuantity * mu_n(2) / var_n(2);
        P_ReferenceSignal_matrix(2*M + InducingPointIdx, AgentNr) = ...
            AgentQuantity / var_n(1);
        P_ReferenceSignal_matrix(3*M + InducingPointIdx, AgentNr) = ...
            AgentQuantity / var_n(2);
    end
end

% 论文式(2.34)(2.35)：静态 AC
L_ac = diag(sum(AdjMatrix, 2)) - AdjMatrix;
lambda_K = max(eig(L_ac));
alpha = 1.9 / lambda_K;  % 论文式(2.35)：0 < alpha < 2/lambda_K
W = eye(AgentQuantity) - alpha * L_ac;  % 混合矩阵


X_consensus = P_ReferenceSignal_matrix;
tol = 1e-10;
for iter = 1:100000
    X_prev = X_consensus;
    X_consensus = X_consensus * W';  % 论文式(2.34)：x(k+1) = W*x(k)
    if norm(X_consensus - X_prev, 'fro') < tol
        fprintf('AC converged at iter %d\n', iter);
        break;
    end
end

% 验证收敛
for i = 1:AgentQuantity
    for j = i+1:AgentQuantity
        fprintf('Xi(%d) - Xi(%d) norm = %.4e\n', i, j, ...
            norm(X_consensus(:,i) - X_consensus(:,j)));
    end
end

% 解码 Xi
Xi_matrix = X_consensus;

idx_x_numerator   = 1:M;
idx_y_numerator   = M+1:2*M;
idx_x_denominator = 2*M+1:3*M;
idx_y_denominator = 3*M+1:4*M;

MaskedGP_TrainingOutput = zeros(y_dim, M, AgentQuantity);
for AgentNr = 1:AgentQuantity
    MaskedGP_TrainingOutput(1, :, AgentNr) = ...
        Xi_matrix(idx_x_numerator, AgentNr)' ./ Xi_matrix(idx_x_denominator, AgentNr)';
    MaskedGP_TrainingOutput(2, :, AgentNr) = ...
        Xi_matrix(idx_y_numerator, AgentNr)' ./ Xi_matrix(idx_y_denominator, AgentNr)';
end

% 建立 MaskedGP
MaskedGP = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    MaskedGP{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-6, SigmaF, SigmaL);
    MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, ...
        MaskedGP_TrainingOutput(:, :, AgentNr));
end
end