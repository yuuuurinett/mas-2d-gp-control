function [h,g,f] = Manipulator_2D_2DoF_get_Dynamics_h_g_f(x,L1,L2,m1,m2)
r = x(1:2);
dr = x(3:4);
%% Convert to joint coordinate
q = Manipulator_2D_2DoF_get_q_from_r(r, L1, L2);
q1 = q(1);
q2 = q(2);

invT = Manipulator_2D_2DoF_get_invT(q, L1, L2);
dq = invT * dr;
dq1 = dq(1);
dq2 = dq(2);
%% Mass and Force in joint coordinate
[Mass_q,Force_q] = Manipulator_2D_2DoF_get_MassForce_q(q,dq,L1,L2,m1,m2);

Force_Coordinate_Convertion_1 = L1 * cos(q1) * dq1^2 + L2 * cos(q1 + q2) * (dq1 + dq2)^2;
Force_Coordinate_Convertion_2 = L1 * sin(q1) * dq1^2 + L2 * sin(q1 + q2) * (dq1 + dq2)^2;
Force_Coordinate_Convertion = [Force_Coordinate_Convertion_1; Force_Coordinate_Convertion_2];
%%
Mass_r = invT' * Mass_q * invT;
Mass_r_safe = Mass_r + 1e-4 * eye(size(Mass_r));

%%
h = Mass_r_safe \ (- invT' * Mass_q * invT * Force_Coordinate_Convertion - ...
	invT' * Force_q);
g = Mass_r_safe \ invT';
f = Manipulator_2D_2DoF_UnknownDynamics(x);

end