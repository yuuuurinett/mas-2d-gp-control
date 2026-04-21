function dx = single_agent_dynamics(t,x,u,q_dim,SystemOrder) 
x_set = reshape(x,[q_dim,SystemOrder]); 
dx_set = [x_set(:,2:SystemOrder),nan(q_dim,1)]; 
dx_set(:,end) = h_unknown(x_set) + u;
dx = reshape(dx_set,[q_dim * SystemOrder,1]); 
end