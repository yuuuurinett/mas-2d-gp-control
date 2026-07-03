function [x_l_set, x_l_r_set] = ET_MAS_GP_Leader_LeaderDynamics(t_set)
t_set = reshape(t_set,[1,numel(t_set)]);
A = 1;
w = 0.4;
x_l_set = [
	A * sin(w * t_set); 
	A * w * cos(w * t_set)];
x_l_r_set = - A * w^2 * sin(w * t_set);
end
