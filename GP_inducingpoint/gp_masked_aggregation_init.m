function [P, p_dim, diag_info] = gp_masked_aggregation_init( ...
    LocalGP_set, AgentQuantity, ...
    NumInducingPoints, InducingPoints_Coordinates, method)
%    InducingPoints_Coordinates : x_dim x M
%    P      : p_dim x AgentQuantity x M
%    p_dim  : dimension of P's first axis (4 or 6)
%    diag_info.min_var(AgentNr)         : min posterior variance across
%                                          all M inducing points and both
%                                          output channels, for this agent
%    diag_info.near_floor_count(AgentNr): number of (channel, point)
%                                          entries within 10x of
%                                          posterior_var_floor

method = lower(method);
M = NumInducingPoints;
prior_var = LocalGP_set{1}.SigmaF^2;  % sigma_0^²
aggregation_cfg = control_aggregation_parameters();
posterior_var_floor = aggregation_cfg.posterior_var_floor;
rbcm_beta_max = aggregation_cfg.rbcm_beta_max;
bcm_prior_scale = aggregation_cfg.bcm_prior_scale;

switch method
    case {'poe', 'gpoe', 'moe', 'bcm', 'rbcm'}
        % Unified layout: [numerator_1; denominator_1; ...].  Scaling each
        % local contribution by AgentQuantity makes the DAC average equal
        % to the centralized sum required by the product-based methods.
        p_dim = 4;
    otherwise
        error('Unknown aggregation method: %s. Choose from poe/gpoe/moe/bcm/rbcm.', method);
end

P = zeros(p_dim, AgentQuantity, M);
diag_info.min_var = inf(AgentQuantity,1);
diag_info.near_floor_count = zeros(AgentQuantity,1);

for AgentNr = 1:AgentQuantity
    [mu_all, var_all] = predict_inducing_points_batch( ...
        LocalGP_set{AgentNr}, InducingPoints_Coordinates);

    % Diagnostic snapshot of this agent's raw posterior variance BEFORE
    % the posterior_var_floor clamp below is applied per-channel, so we
    % can tell whether variance is approaching the floor over time.
    diag_info.min_var(AgentNr) = min(var_all(:));
    diag_info.near_floor_count(AgentNr) = nnz(var_all(:) < posterior_var_floor*10);

    for InducingPointIdx = 1:M
        mu_n = mu_all(:,InducingPointIdx);
        var_n = var_all(:,InducingPointIdx);

        mu1  = max(-30, min(30, mu_n(1)));
        var1 = max(var_n(1), posterior_var_floor);
        mu2  = max(-30, min(30, mu_n(2)));
        var2 = max(var_n(2), posterior_var_floor);

        switch method

            case 'poe'                        
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * mu1 / var1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity / var1;
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * mu2 / var2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity / var2;

            case 'gpoe'        
                beta1 = max(min(0.5 * (log(prior_var) - log(var1)), ...
                    aggregation_cfg.gpoe_beta_max),eps);
                beta2 = max(min(0.5 * (log(prior_var) - log(var2)), ...
                    aggregation_cfg.gpoe_beta_max),eps);
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * beta1 * mu1 / var1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity * beta1 / var1;
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * beta2 * mu2 / var2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity * beta2 / var2;

            case 'moe'
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * mu1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity * (var1 + mu1^2);
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * mu2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity * (var2 + mu2^2);

            case 'bcm'
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * mu1 / var1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity / var1 - ...
                    bcm_prior_scale * (AgentQuantity - 1) / prior_var;
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * mu2 / var2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity / var2 - ...
                    bcm_prior_scale * (AgentQuantity - 1) / prior_var;

            case 'rbcm'
                raw_beta1 = 0.5 * (log(prior_var) - log(var1));
                raw_beta2 = 0.5 * (log(prior_var) - log(var2));
                beta1 = min(max(raw_beta1, eps), rbcm_beta_max);
                beta2 = min(max(raw_beta2, eps), rbcm_beta_max);
                P(1, AgentNr, InducingPointIdx) = AgentQuantity * beta1 * mu1 / var1;
                P(2, AgentNr, InducingPointIdx) = AgentQuantity * beta1 / var1 + ...
                    (1 - AgentQuantity * beta1) / prior_var;
                P(3, AgentNr, InducingPointIdx) = AgentQuantity * beta2 * mu2 / var2;
                P(4, AgentNr, InducingPointIdx) = AgentQuantity * beta2 / var2 + ...
                    (1 - AgentQuantity * beta2) / prior_var;
        end
    end
end
end

function [mu, var] = predict_inducing_points_batch(gp, X_query)
% Vectorized equivalent of repeated gp.predict calls for mean/variance.
M = size(X_query,2);
if gp.DataQuantity == 0
    mu = zeros(gp.y_dim,M);
    var = gp.SigmaF^2 * ones(gp.y_dim,M);
    return;
end

N = gp.DataQuantity;
X_train = gp.X(:,1:N);
K_star = gp.kernel(X_train,X_query);
L = gp.L(1:N,1:N);
V = L \ K_star;

mu = gp.alpha(1:N,:)' * K_star;
var_scalar = max(gp.SigmaF^2 - sum(V.^2,1),0);
var = repmat(var_scalar,gp.y_dim,1);
end
