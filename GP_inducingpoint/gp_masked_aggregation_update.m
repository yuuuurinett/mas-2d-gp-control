function [MaskedGP, Zeta_vector] = gp_masked_aggregation_update( ...
    P, Zeta_vector, L, ...
    Kappa_P, AgentQuantity, NumInducingPoints, TimeStep, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, method, p_dim)

method  = lower(method);
M       = NumInducingPoints;
y_dim   = 2;
prior_var = SigmaF^2;   

%% Step 1: DAC integration over [0, TimeStep]
if TimeStep > 0      
    for InducingPointIdx = 1:M
        P_InducingPoint    = P(:, :, InducingPointIdx);
        Zeta_InducingPoint = Zeta_vector(:, :, InducingPointIdx);

        New_Consensus_Zeta_function = @(~, Zeta_ODE_Intern) Compute_New_Consensus_Derivative( ...
            Zeta_ODE_Intern, P_InducingPoint, L, Kappa_P, AgentQuantity, p_dim);

        [~, Zeta_ODE_Output] = ode45(New_Consensus_Zeta_function, ...
            [0, TimeStep], Zeta_InducingPoint(:));

        Zeta_vector(:, :, InducingPointIdx) = reshape(Zeta_ODE_Output(end,:)', ...
            p_dim, AgentQuantity);
    end
end

%% Step 2: compute xi

Xi_all = P - Zeta_vector;   % p_dim x AgentQuantity x M

switch method
    case {'poe', 'gpoe', 'moe'}
        num1 = squeeze(Xi_all(1, :, :));  % AgentQuantity x M
        den1 = squeeze(Xi_all(2, :, :));
        num2 = squeeze(Xi_all(3, :, :));
        den2 = squeeze(Xi_all(4, :, :));

        phi1 = num1 ./ den1;   % AgentQuantity x M
        phi2 = num2 ./ den2;

    case 'bcm'
        num1 = squeeze(Xi_all(1, :, :));
        num2 = squeeze(Xi_all(2, :, :));
        den1 = squeeze(Xi_all(3, :, :));
        den2 = squeeze(Xi_all(4, :, :));
 
        prior_correction = (1 - AgentQuantity) / prior_var;
        den1_fused = den1 + prior_correction;
        den2_fused = den2 + prior_correction;
        
        % 🌟 优雅回退：如果分母小于 0.01（极度不自信或过零），直接输出 0！
        phi1 = zeros(size(num1));
        phi2 = zeros(size(num2));
        mask1 = den1_fused > 1e-2;
        mask2 = den2_fused > 1e-2;
        
        phi1(mask1) = num1(mask1) ./ den1_fused(mask1);
        phi2(mask2) = num2(mask2) ./ den2_fused(mask2);

    case 'rbcm'
       
        num1  = squeeze(Xi_all(1, :, :));
        den1  = squeeze(Xi_all(2, :, :));
        beta1 = squeeze(Xi_all(3, :, :));
        num2  = squeeze(Xi_all(4, :, :));
        den2  = squeeze(Xi_all(5, :, :));
        beta2 = squeeze(Xi_all(6, :, :));

        den1_fused = den1 + (1 - beta1) / prior_var;
        den2_fused = den2 + (1 - beta2) / prior_var;
        
       
        phi1 = zeros(size(num1));
        phi2 = zeros(size(num2));
        mask1 = den1_fused > 1e-2;
        mask2 = den2_fused > 1e-2;
        
        phi1(mask1) = num1(mask1) ./ den1_fused(mask1);
        phi2(mask2) = num2(mask2) ./ den2_fused(mask2);
end

%% Step 3: masked GP
MaskedGP = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    Y_agent = [phi1(AgentNr, :); phi2(AgentNr, :)];  % 2 x M
    MaskedGP{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-4, SigmaF, SigmaL);
    MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, Y_agent);
end

end

%% compute zeta_dot
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_Ref, L, Kappa, AgentQuantity, p_dim)

Zeta = reshape(Zeta_vec, p_dim, AgentQuantity);
dZeta_dt = Kappa * (P_Ref - Zeta) * L';
dZeta_dt = dZeta_dt(:);
end