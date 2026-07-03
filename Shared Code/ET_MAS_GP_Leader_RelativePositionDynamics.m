function [s_i_set, s_i_r_set] = ET_MAS_GP_Leader_RelativePositionDynamics(t_set, AgentNr)
t_set = reshape(t_set,[1,numel(t_set)]);
A = 0.01;
w = 6;
s_i_set = [
	A * sin(w * t_set + AgentNr * pi / 2); 
	A * w * cos(w * t_set + AgentNr * pi / 2)];
s_i_r_set = - A * w^2 * sin(w * t_set + AgentNr * pi / 2);
end