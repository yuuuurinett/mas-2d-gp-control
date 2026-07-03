function [s_i_set, s_i_r_set] = ET_MAS_GP_Leader_RelativePositionDynamics_2D(t_set, AgentNr)
    t_set = reshape(t_set, [1, numel(t_set)]);

    A = 0.2;        
    w = 1.5;        
    phase = AgentNr * pi / 3; 

    s_i1 = [A * cos(w * t_set + phase); 
            A * sin(w * t_set + phase)];

    s_i2 = [-A * w * sin(w * t_set + phase); 
             A * w * cos(w * t_set + phase)];

    s_i_set = [s_i1; s_i2]; 

    s_i_r_set = [-A * w^2 * cos(w * t_set + phase); 
                 -A * w^2 * sin(w * t_set + phase)];
end