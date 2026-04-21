classdef Agent_Leader_Topology_Class < handle

	properties
		AgentQuantity;
		LeaderQuantity;
		ConnectionMatrix;
		AgentNeighbour;
	end

	methods
		function obj = Agent_Leader_Topology_Class(AgentQuantity_in,LeaderQuantity_in)
			obj.AgentQuantity = AgentQuantity_in;
			obj.LeaderQuantity = LeaderQuantity_in;
			obj.ConnectionMatrix = zeros(obj.AgentQuantity,obj.LeaderQuantity);
			obj.AgentNeighbour = cell(obj.AgentQuantity,1);
		end
		%% Add Connection %每条边要调用一次
		function addConnection(obj,NodeType1,NodeNr1,NodeType2,NodeNr2, ...
				ConnectionWeight)
			if nargin < 6  %Number of Arguments IN
				ConnectionWeight = 1;
			end
			% Check whether the connection is valid
			if strcmpi(NodeType1,'Leader') && strcmpi(NodeType2,'Agent')
				LeaderNodeNr = NodeNr1;
				AgentNodeNr = NodeNr2;
			elseif strcmpi(NodeType1,'Agent') && strcmpi(NodeType2,'Leader')
				LeaderNodeNr = NodeNr2;
				AgentNodeNr = NodeNr1;
			else
				error('Not a valid connection between leader and agent!');
			end
			% Check wether the node number is valid(对比输入的编号是否超过了最初定义的数量)
			if LeaderNodeNr > obj.LeaderQuantity || AgentNodeNr > obj.AgentQuantity
				warning('The node number exceed the quantity, no connection is built!');
				return;
			end
			%
			obj.ConnectionMatrix(AgentNodeNr,LeaderNodeNr) = ConnectionWeight; %行：agent，列：leader(存储 agent 和 leader 之间的连接权重)
			% Configure the Agent neighbour information
			ExistedConnectionNr = find(obj.AgentNeighbour{AgentNodeNr} == LeaderNodeNr,1); %所有agent追踪的leader
			if isempty(ExistedConnectionNr)
				obj.AgentNeighbour{AgentNodeNr} = [obj.AgentNeighbour{AgentNodeNr};LeaderNodeNr];
			end
		end
		%% Add Connection Set
		function addConnectionSet(obj,ConnectionConfigurationSet)
			for ConnectionNr = 1:numel(ConnectionConfigurationSet)
				ConnectionConfiguration = ConnectionConfigurationSet{ConnectionNr};
				NodeType1 = ConnectionConfiguration{1};
				NodeNr1   = ConnectionConfiguration{2};
				NodeType2 = ConnectionConfiguration{3};
				NodeNr2   = ConnectionConfiguration{4};
				weight    = ConnectionConfiguration{5};
				obj.addConnection(NodeType1,NodeNr1,NodeType2,NodeNr2,weight);
			end
		end
	end
end
