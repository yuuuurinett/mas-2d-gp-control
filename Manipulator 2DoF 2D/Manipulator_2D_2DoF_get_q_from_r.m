function q = Manipulator_2D_2DoF_get_q_from_r(r, L1, L2)

x = r(1);
y = r(2);
L = norm(r);
La = (L^2 + L1^2 - L2^2) / (2 * L);
Lb = (L^2 - L1^2 + L2^2) / (2 * L);
a = max(0,real(sqrt(L1^2 - La^2)));
e = r / L;
R = [0,1;-1,0];
ea = R * e;
rB = La * e + a * ea;
q1 = atan2(rB(2), rB(1));

R1 = get_R_z(-q1);
R1 = R1(1:2,1:2);
r_rotated = R1 * r;
dr = r_rotated - R1 * rB;
q2 = atan2(dr(2), dr(1));

q = [q1;q2];

% FigureObj_r_to_q = figure('Name','Test from r to q');
% AxesObj_r_to_q = axes(FigureObj_r_to_q);
% rB = L1 * [cos(q1); sin(q1)];
% rC = rB + L2 * [cos(q1 + q2); sin(q1 + q2)];
% hold(AxesObj_r_to_q,'on');
% plot(AxesObj_r_to_q,[0,x],[0,y],'k.--');
% plot(AxesObj_r_to_q,[0,rB(1),rC(1)],[0,rB(2),rC(2)],'r-o');
% axis(AxesObj_r_to_q,[-(L1 + L2), (L1 + L2),-(L1 + L2), (L1 + L2)]);
end