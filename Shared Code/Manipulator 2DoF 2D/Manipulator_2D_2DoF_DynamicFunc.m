function dx = Manipulator_2D_2DoF_DynamicFunc(t,x,u,L1,L2,m1,m2,unknown_scale)

if nargin < 8 || isempty(unknown_scale)
	unknown_scale = 1;
end

dr = x(3:4);

[h,g,f] = Manipulator_2D_2DoF_get_Dynamics_h_g_f(x,L1,L2,m1,m2);

ddr = h + g * u + unknown_scale * f;

dx = [dr;ddr];
end
