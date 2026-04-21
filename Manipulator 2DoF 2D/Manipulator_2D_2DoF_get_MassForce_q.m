function [Mass_q,Force_q] = Manipulator_2D_2DoF_get_MassForce_q(q,dq,L1,L2,m1,m2)
g = 9.8;
%%
q1 = q(1);
q2 = q(2);
dq1 = dq(1);
dq2 = dq(2);
%% Mass Matrix
Mass_q_11 = 1/3 * m1 * L1^2 + m2 * L1^2 + 1/3 * m2 * L2^2 + m2 * L1 * L2 * cos(q1);
Mass_q_12 = 1/3 * m2 * L2^2 + 1/2 * m2 * L1 * L2 * cos(q1);
Mass_q_21 = Mass_q_12;
Mass_q_22 = 1/3 * m2 * L2^2;
Mass_q = [Mass_q_11, Mass_q_12; Mass_q_21, Mass_q_22];
%% Force Vector
Force_Damping_1 = - (2 * dq1 * dq2 + dq2^2);
Force_Damping_2 = dq1^2;
Force_Damping = 1/2 * m2 * L1 * L2 * sin(q2) * [Force_Damping_1; Force_Damping_2];

Force_Gravity_1 = 1/2 * m1 * L1 * cos(q1) + m2 * L1 * cos(q1) + 1/2 * m2 * L2 * cos(q1 + q2);
Force_Gravity_2 = 1/2 * m2 * L2 * cos(q1 + q2);
Force_Gravity = [Force_Gravity_1; Force_Gravity_2] * g;

Force_q = Force_Damping + Force_Gravity;
end