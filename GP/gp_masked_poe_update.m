function [MaskedGP, Zeta_vector] = gp_masked_poe_update( ...
    P, Zeta_vector,  L, ...
    Kappa_P, AgentQuantity, NumInducingPoints, TimeStep, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim)

M=NumInducingPoints;
y_dim=2;

%% Step1: update dac for each inducing point separately
if TimeStep > 0
    for InducingPointIdx = 1:M
        P_InducingPoint    = P(:, :, InducingPointIdx);
        Zeta_InducingPoint = Zeta_vector(:, :, InducingPointIdx);

        New_Consensus_Zeta_function = @(~, Zeta_ODE_Intern) Compute_New_Consensus_Derivative( ...
            Zeta_ODE_Intern, P_InducingPoint, L, Kappa_P, AgentQuantity);

        [~, Zeta_ODE_Output] = ode45(New_Consensus_Zeta_function, ...
            [0, TimeStep], Zeta_InducingPoint(:));

        Zeta_vector(:, :, InducingPointIdx) = reshape(Zeta_ODE_Output(end,:)', ...
            4, AgentQuantity);
    end
end
%% Step2: compute xi 
Xi_all = P - Zeta_vector;  % 4 x AgentQuantity x M

numerator_x   = squeeze(Xi_all(1,:,:));  % AgentQuantity x M
denominator_x = squeeze(Xi_all(3,:,:));
numerator_y   = squeeze(Xi_all(2,:,:));
denominator_y = squeeze(Xi_all(4,:,:));

phi_x = numerator_x ./ denominator_x;
phi_y = numerator_y ./ denominator_y;

MaskedGP = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    Y_agent = [phi_x(AgentNr, :); phi_y(AgentNr, :)];
    MaskedGP{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, ...
        M, 1e-6, SigmaF, SigmaL);
    MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, Y_agent);
end

end
%%
function dZeta_dt = Compute_New_Consensus_Derivative(...
    Zeta_vec, P_InducingPoint, L, Kappa, AgentQuantity)

Zeta_Matrix = reshape(Zeta_vec, 4, AgentQuantity);
dZeta_dt = Kappa * (P_InducingPoint - Zeta_Matrix) * L';
dZeta_dt = dZeta_dt(:);

end
