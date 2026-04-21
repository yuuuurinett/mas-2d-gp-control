function [si,sir] = Manipulator_2D_2DoF_RelativePositionDynamics(t,AgentNr,AgentQuantity)
w = 1.5;
L = 0.2;
phi = AgentNr / AgentQuantity * 2 * pi;
%% Leader State
rxl = L * cos(w * t + phi);
ryl = L * sin(w * t + phi);
drxl = - w * L * sin(w * t + phi);
dryl =   w * L * cos(w * t + phi);
si = [rxl;ryl;drxl;dryl];
%% Leader xlr
xlr_x = - w^2 * L * cos(w * t + phi);
xlr_y = - w^2 * L * sin(w * t + phi);
sir = [xlr_x;xlr_y];
end