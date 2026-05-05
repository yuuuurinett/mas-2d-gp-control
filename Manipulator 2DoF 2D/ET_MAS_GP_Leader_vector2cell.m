function a_cell = ET_MAS_GP_Leader_vector2cell(a_all, AgentQuantity, SystemOrder)
a_cell = cell(AgentQuantity,SystemOrder);
a_matrix = reshape(a_all,[],AgentQuantity);
for AgentNr = 1:AgentQuantity
	a_i = a_matrix(:,AgentNr);
	a_i_matrix = reshape(a_i,[],SystemOrder);
	for SystemOrderNr = 1:SystemOrder
		a_cell{AgentNr,SystemOrderNr} = a_i_matrix(:,SystemOrderNr);
	end
end
end