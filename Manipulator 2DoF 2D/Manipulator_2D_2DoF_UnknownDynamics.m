function f = Manipulator_2D_2DoF_UnknownDynamics(x)

r = x(1:2,:);
dr = x(3:4,:);

rx = r(1,:);
ry = r(2,:);
drx = dr(1,:);
dry = dr(2,:);

f_1 = 5 * sin(rx) + 3 * cos(rx) + drx.^2 + 6;
f_2 = 3 * cos(ry) + 5 * sin(ry) + dry.^2 + 10;
f = [f_1; f_2];

end 