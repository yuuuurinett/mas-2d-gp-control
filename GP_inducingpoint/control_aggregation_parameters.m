function cfg = control_aggregation_parameters()
%CONTROL_AGGREGATION_PARAMETERS Shared control-task aggregation constants.
% All IP, TP, centralized, and neighbor implementations must use these
% values so that method labels differ only in support/query architecture.

cfg.posterior_var_floor = positive_env( ...
    'CONTROL_AGGREGATION_VARIANCE_FLOOR',1e-3);
cfg.precision_floor = positive_env( ...
    'CONTROL_AGGREGATION_PRECISION_FLOOR',1e-4);
cfg.gpoe_beta_max = positive_env('CONTROL_GPOE_BETA_MAX',10);
cfg.rbcm_beta_max = positive_env('CONTROL_RBCM_BETA_MAX',0.25);
cfg.bcm_prior_scale = positive_env('CONTROL_BCM_PRIOR_SCALE',0.25);
end

function value = positive_env(name,default_value)
value = str2double(getenv(name));
if ~(isfinite(value) && value > 0)
    value = default_value;
end
end
