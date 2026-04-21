function [xl,xlr,w] = Manipulator_2D_2DoF_LeaderDynamics(t,L1)
w = 0.5;
%% Leader State
rxl = L1 * cos(w * t);
ryl = L1 * sin(w * t);
drxl = - w * L1 * sin(w * t);
dryl =   w * L1 * cos(w * t);
xl = [rxl;ryl;drxl;dryl];
%% Leader xlr
xlr_x = - w^2 * L1 * cos(w * t);
xlr_y = - w^2 * L1 * sin(w * t);
xlr = [xlr_x;xlr_y];
end