function [MaskedGP, Zeta_vector, Zeta_last_trigger, trigger_count] = gp_masked_aggregation_update( ...
    P, Zeta_vector, L, ...
    Kappa_P, AgentQuantity, NumInducingPoints, TimeStep, ...
    InducingPoints_Coordinates, SigmaF, SigmaL, x_dim, method, p_dim, ...
    Zeta_last_trigger, neighbor_count)
%GP_MASKED_AGGREGATION_UPDATE_FIXED
% Dataset-consistent IP-DAC update for the control task.
%
% IMPORTANT VARIABLE CONVENTION
%   P(:,i,m)                 : local GP information vector of agent i at inducing point m
%   Zeta_vector(:,i,m)       : DAC auxiliary state zeta_i^m
%   Xi_now(:,i,m)            : current DAC output, xi_i^m = P_i^m - zeta_i^m
%   Zeta_last_trigger(:,i,m) : last transmitted xi_hat_i^m
%
% The name Zeta_last_trigger is kept only for compatibility with the
% existing run_simulation_inducing_point.m interface. Mathematically it is
% Xi_last_trigger, not zeta_last_trigger.
%
% DAC dynamics used here:
%   dot{zeta} = Kappa_P * L * xi_hat
%   xi        = P - zeta
%
% Event trigger:
%   point-wise per agent and per inducing point.
%   The computation is vectorized over all inducing points, but the trigger
%   decision is still made separately for each inducing point.

method    = lower(method);
M         = NumInducingPoints;
y_dim     = 2;
prior_var = SigmaF^2;

trigger_count = zeros(AgentQuantity, 1);

%% Step 1: Dataset-consistent IP-DAC update
% Xi_last_trigger is stored in Zeta_last_trigger for backward compatibility.
% Consensus dynamics uses last-triggered xi_hat, not current xi.
if TimeStep > 0
    L_Xi_hat = laplacian_multiply_agent_dim(Zeta_last_trigger, L);
    Zeta_vector = Zeta_vector + TimeStep * Kappa_P * L_Xi_hat;
end

% Current DAC output xi = P - zeta
Xi_now = P - Zeta_vector;

%% Step 2: Point-wise event-triggered communication
% This is NOT one trigger for all inducing points. For each agent i, the
% trigger_idx below is a logical vector of length M.
%
% We use the undirected-graph Kia-style trigger structure based on the
% broadcast disagreement term:
%   ||xi_hat_i^m - xi_i^m||^2
%      > 1/(4 d_i) * sum_j ||xi_hat_i^m - xi_hat_j^m||^2
%        + 1/(4 d_i) * epsilon_i^2
%
% This is evaluated independently for every inducing point m.
epsilon_i = 0.01;

if TimeStep > 0
    for agent_i = 1:AgentQuantity
        d_i_out = neighbor_count(agent_i);
        neighbor_mask_i = (L(agent_i,:) < 0);

        if d_i_out <= 0 || ~any(neighbor_mask_i)
            continue;
        end

        % e_i^m = xi_hat_i^m - xi_i^m
        E_i = Zeta_last_trigger(:, agent_i, :) - Xi_now(:, agent_i, :);
        e_norm_sq = squeeze(sum(E_i.^2, 1));        % [M x 1] or [1 x M]
        e_norm_sq = e_norm_sq(:).';                 % [1 x M]

        % Kia-style broadcast neighbor disagreement:
        % sum_{j in N_i} ||xi_hat_i^m - xi_hat_j^m||^2
        Xi_hat_i = Zeta_last_trigger(:, agent_i, :);              % [p_dim x 1 x M]
        Xi_hat_neighbors = Zeta_last_trigger(:, neighbor_mask_i, :); % [p_dim x d_i x M]
        Diff_hat = Xi_hat_i - Xi_hat_neighbors;                   % implicit expansion
        neighbor_diff_sq_each = squeeze(sum(Diff_hat.^2, 1));     % [d_i x M]

        if isvector(neighbor_diff_sq_each)
            % If there is only one neighbor, squeeze may return [M x 1].
            z_hat_sq = neighbor_diff_sq_each(:).';
        else
            z_hat_sq = sum(neighbor_diff_sq_each, 1);             % [1 x M]
        end

        threshold = (1/(4*d_i_out)) * z_hat_sq + ...
                    (1/(4*d_i_out)) * epsilon_i^2;

        trigger_idx = e_norm_sq > threshold;

        if any(trigger_idx)
            % Only update the triggered inducing-point snapshots.
            Zeta_last_trigger(:, agent_i, trigger_idx) = Xi_now(:, agent_i, trigger_idx);
            trigger_count(agent_i) = nnz(trigger_idx);
        end
    end
end

%% Step 3: Reconstruct MaskedGP from aggregated xi information
Xi_all = Xi_now;

switch method
    case {'poe', 'gpoe', 'moe'}
        num1 = squeeze(Xi_all(1, :, :));
        den1 = squeeze(Xi_all(2, :, :));
        num2 = squeeze(Xi_all(3, :, :));
        den2 = squeeze(Xi_all(4, :, :));
        phi1 = safe_divide(num1, den1);
        phi2 = safe_divide(num2, den2);

    case 'bcm'
        num1 = squeeze(Xi_all(1, :, :));
        num2 = squeeze(Xi_all(2, :, :));
        den1 = squeeze(Xi_all(3, :, :));
        den2 = squeeze(Xi_all(4, :, :));
        prior_correction = (1 - AgentQuantity) / prior_var;
        den1_fused = den1 + prior_correction;
        den2_fused = den2 + prior_correction;
        phi1 = zeros(size(num1));
        phi2 = zeros(size(num2));
        mask1 = abs(den1_fused) > 1e-8;
        mask2 = abs(den2_fused) > 1e-8;
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
        mask1 = abs(den1_fused) > 1e-8;
        mask2 = abs(den2_fused) > 1e-8;
        phi1(mask1) = num1(mask1) ./ den1_fused(mask1);
        phi2(mask2) = num2(mask2) ./ den2_fused(mask2);

    otherwise
        error('Unknown aggregation method: %s', method);
end

phi1(~isfinite(phi1)) = 0;
phi2(~isfinite(phi2)) = 0;

% Debug safety clamp. You can loosen/remove this after stability is verified.
phi_clip = 30;
phi1 = max(-phi_clip, min(phi_clip, phi1));
phi2 = max(-phi_clip, min(phi_clip, phi2));

MaskedGP = cell(AgentQuantity, 1);
for AgentNr = 1:AgentQuantity
    Y_agent = [phi1(AgentNr, :); phi2(AgentNr, :)];
    MaskedGP{AgentNr} = LocalGP_MultiOutput(x_dim, y_dim, M, 1e-4, SigmaF, SigmaL);
    MaskedGP{AgentNr}.add_Alldata(InducingPoints_Coordinates, Y_agent);
end

end

%% Helper: L acts only on the agent dimension
function L_X = laplacian_multiply_agent_dim(X, L)
% X: [p_dim x AgentQuantity x NumPoints]
% L: [AgentQuantity x AgentQuantity]
[p_dim, AgentQuantity, NumPoints] = size(X);

X_agent_first = permute(X, [2, 1, 3]);              % [Agent x p_dim x M]
X_flat = reshape(X_agent_first, AgentQuantity, []); % [Agent x (p_dim*M)]

L_X_flat = L * X_flat;

L_X_agent_first = reshape(L_X_flat, AgentQuantity, p_dim, NumPoints);
L_X = permute(L_X_agent_first, [2, 1, 3]);          % [p_dim x Agent x M]
end

%% Helper: numerically safe division
function out = safe_divide(num, den)
out = zeros(size(num));
mask = abs(den) > 1e-8;
out(mask) = num(mask) ./ den(mask);
end
