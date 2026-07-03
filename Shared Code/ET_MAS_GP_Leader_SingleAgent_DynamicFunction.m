function dx = ET_MAS_GP_Leader_SingleAgent_DynamicFunction(t,x,u,q_dim,SystemOrder)
x_set = reshape(x,[q_dim,SystemOrder]);
dx_set = [x_set(:,2:SystemOrder),nan(q_dim,1)];
fx = ET_MAS_GP_Leader_UnknownDynamics(x,q_dim,SystemOrder);
dx_set(:,end) = fx + u;
dx = reshape(dx_set,[q_dim * SystemOrder,1]);
end