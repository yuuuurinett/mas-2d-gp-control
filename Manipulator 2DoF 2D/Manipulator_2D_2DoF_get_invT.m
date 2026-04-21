function invT = Manipulator_2D_2DoF_get_invT(q, L1, L2)

q1 = q(1);
q2 = q(2);

invT_11 = L2 * cos(q1 + q2);
invT_12 = L2 * sin(q1 + q2);
invT_21 = - L1 * cos(q1) - L2 * cos(q1 + q2);
invT_22 = - L1 * sin(q1) - L2 * sin(q1 + q2);

invT = 1 / (L1 * L2 * sin(q2)) * [invT_11, invT_12; invT_21, invT_22];

end