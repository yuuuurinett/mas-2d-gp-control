function [Epsilon_vector,r_n] = consensus_error( ...
    AgentNr, AgentState_matrix, LeaderState_vector, ...
    AdjMatrix, NeighbourSet, B,lambda_set, SystemOrder, q_dim)
    b_n=B(AgentNr);
    AgentNeighbourSet = NeighbourSet{AgentNr};
    AgentNeighbourQuantity = numel(AgentNeighbourSet);
    r_n = zeros(q_dim, 1); 
    Epsilon_vector = zeros(q_dim,SystemOrder);
   
    for SystemOrderNr = 1:SystemOrder
    e_ik = -b_n * (AgentState_matrix(SystemOrderNr*q_dim-1:SystemOrderNr*q_dim, AgentNr) ...
                 - LeaderState_vector(SystemOrderNr*q_dim-1:SystemOrderNr*q_dim));
    
    for NeighborNr = 1:AgentNeighbourQuantity
        NeighbourAgentNr = AgentNeighbourSet(NeighborNr);
        a_ij = AdjMatrix(AgentNr, NeighbourAgentNr);
        e_ik = e_ik - a_ij * (AgentState_matrix(SystemOrderNr*q_dim-1:SystemOrderNr*q_dim, AgentNr) ...
                             - AgentState_matrix(SystemOrderNr*q_dim-1:SystemOrderNr*q_dim, NeighbourAgentNr));
    end
    
    Epsilon_vector(:, SystemOrderNr) = e_ik;
    lambda_k = lambda_set(SystemOrderNr);
    r_n = r_n + lambda_k * e_ik;
    end
   
end