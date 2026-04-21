function Nu_n = controller(...
    r_n, Epsilon_vector,  Kappa_C, lambda_set, SystemOrder)
Nu_n = Kappa_C * r_n;
for SystemOrderNr = 1:SystemOrder-1
    Nu_n=Nu_n+lambda_set(SystemOrderNr) * Epsilon_vector(SystemOrderNr+1);
end
end