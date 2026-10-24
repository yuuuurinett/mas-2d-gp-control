function cfg = control_aggregation_parameters(SigmaN)
%CONTROL_AGGREGATION_PARAMETERS Shared control-task aggregation constants.
% All IP, TP, centralized, and neighbor implementations must use these
% values so that method labels differ only in support/query architecture.
% Following the dataset implementation, inverse-variance aggregation uses
% the observation-noise variance SigmaN^2 as its default variance floor.

if nargin < 1 || isempty(SigmaN)
    SigmaN = 0.05;
end
if ~(isscalar(SigmaN) && isfinite(SigmaN) && SigmaN > 0)
    error('SigmaN must be a finite positive scalar.');
end

cfg.posterior_var_floor = positive_env( ...
    'CONTROL_AGGREGATION_VARIANCE_FLOOR',SigmaN^2);
cfg.posterior_var_floor_default = SigmaN^2;
cfg.posterior_var_floor_uses_override = ...
    isfinite(str2double(getenv('CONTROL_AGGREGATION_VARIANCE_FLOOR'))) && ...
    str2double(getenv('CONTROL_AGGREGATION_VARIANCE_FLOOR')) > 0;
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
