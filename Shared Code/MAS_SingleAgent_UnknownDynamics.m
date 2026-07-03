function f = MAS_SingleAgent_UnknownDynamics(x,q_dim,SystemOrder)
x1 = x(1,:);
x2 = x(2,:);
f = sin(x1) + 0.5 ./ (1 + exp(- x2 / 10));
end