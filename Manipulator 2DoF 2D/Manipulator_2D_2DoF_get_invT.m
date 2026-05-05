function invT = Manipulator_2D_2DoF_get_invT(q, L1, L2)

q1 = q(1);
q2 = q(2);

invT_11 = L2 * cos(q1 + q2);
invT_12 = L2 * sin(q1 + q2);
invT_21 = - L1 * cos(q1) - L2 * cos(q1 + q2);
invT_22 = - L1 * sin(q1) - L2 * sin(q1 + q2);

sin_q2 = sin(q2);
sin_q2_safe = sign(sin_q2 + eps) * max(abs(sin_q2), 0.05);
invT = 1 / (L1 * L2 * sin_q2_safe) * [invT_11, invT_12; invT_21, invT_22];

%invT = 1 / (L1 * L2 * sin(q2)) * [invT_11, invT_12; invT_21, invT_22];
%invT = 1 / (L1 * L2 * (sin(q2) ) * [invT_11, invT_12; invT_21, invT_22];
end