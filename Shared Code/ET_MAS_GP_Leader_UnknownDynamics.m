function f = ET_MAS_GP_Leader_UnknownDynamics(x,q_dim,SystemOrder)
x_i_cell = cell(SystemOrder,1);
for SystemOrderNr = 1:SystemOrder
	x_i_cell{SystemOrderNr} = x((SystemOrderNr-1)*q_dim+(1:q_dim),:);
end
% x1 = x(1,:);
% x2 = x(2,:);
f = 5 * sin(10 * x_i_cell{1}) + 0.5 ./ (1 + exp(- x_i_cell{2} / 10)) + 5;
end