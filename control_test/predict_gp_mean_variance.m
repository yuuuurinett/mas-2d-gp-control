function [mu, var] = predict_gp_mean_variance(gp, X_query)
% GP posterior mean/variance without computing unused error bounds.

query_count = size(X_query,2);
if gp.DataQuantity == 0
    mu = zeros(gp.y_dim,query_count);
    var = gp.SigmaF^2 * ones(gp.y_dim,query_count);
    return;
end

N = gp.DataQuantity;
K_star = gp.kernel(gp.X(:,1:N),X_query);
V = gp.L(1:N,1:N) \ K_star;

mu = gp.alpha(1:N,:)' * K_star;
var_scalar = max(gp.SigmaF^2 - sum(V.^2,1),0);
var = repmat(var_scalar,gp.y_dim,1);
end
