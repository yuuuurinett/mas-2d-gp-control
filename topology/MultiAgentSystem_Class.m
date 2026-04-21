classdef MultiAgentSystem_Class < handle
	properties
		AgentQuantity;
		LeaderQuantity;
		q_dim;
		SystemOrder;
		x_dim;
		x;
		t;
		t_step;

		Agent_Topology;
		AgentSet;

		Agent_Leader_Topology;
		LeaderSet;
		Extended_Topology;

		ControlGainSet;
	end

	methods
		function obj = MultiAgentSystem_Class(AgentQuantity_in,LeaderQuantity_in)
			obj.AgentQuantity = AgentQuantity_in;
			obj.LeaderQuantity = LeaderQuantity_in;
			obj.Agent_Topology = Topology_Class(obj.AgentQuantity);
			obj.Agent_Leader_Topology = Agent_Leader_Topology_Class( ...
				obj.AgentQuantity,obj.LeaderQuantity);
			obj.Extended_Topology = Topology_Class(obj.AgentQuantity + obj.LeaderQuantity);
			obj.AgentSet = cell(obj.AgentQuantity,1);
			obj.LeaderSet = cell(obj.LeaderQuantity,1);
		end
		%% Intitialize Agent Set
		function initialize_AgentSet(obj,q_dim_in,SystemOrder_in,DynamicsFunction_in, ...
				t_start,t_step,x_init)
			obj.q_dim = q_dim_in;
			obj.SystemOrder = SystemOrder_in;

			obj.ControlGainSet.AverageConsensusGain = ones(obj.SystemOrder,1);
			obj.ControlGainSet.LeaderFollowerGain = ones(obj.SystemOrder,1);

			obj.x_dim = obj.q_dim * obj.SystemOrder;
			obj.x = x_init;
			x_init = reshape(x_init,[obj.q_dim * obj.AgentQuantity,obj.SystemOrder]);
			for AgentNr = 1:obj.AgentQuantity
				xi_init = x_init((AgentNr - 1) * obj.q_dim + (1:obj.q_dim),:);
				xi_init = reshape(xi_init,[obj.x_dim,1]);

				obj.AgentSet{AgentNr} = Agent_Class(obj.q_dim,obj.SystemOrder,DynamicsFunction_in);
				obj.AgentSet{AgentNr}.set_AgentInitialState(t_start,xi_init);
				obj.AgentSet{AgentNr}.TimeStep = t_step;
				obj.AgentSet{AgentNr}.ControlGainSet = obj.ControlGainSet;
				obj.AgentSet{AgentNr}.AgentNumber = AgentNr;
				obj.AgentSet{AgentNr}.prepare_StateSharing;
			end
			
			obj.t = t_start;
			obj.t_step = t_step;
		end
		%% Get Extended Topology
		function get_ExtendedTopology(obj) %(N+1)*(N+1)
			obj.Extended_Topology.AdjacencyMatrix(1:obj.AgentQuantity,1:obj.AgentQuantity) = ...
				obj.Agent_Topology.AdjacencyMatrix; %N*N
			obj.Extended_Topology.AdjacencyMatrix(1:obj.AgentQuantity,(obj.AgentQuantity+1):end) = ...
				obj.Agent_Leader_Topology.ConnectionMatrix; %1*(N+1)/最后一列是B
			obj.Extended_Topology.AdjacencyMatrix((obj.AgentQuantity+1):end,1:obj.AgentQuantity) = ...
				transpose(obj.Agent_Leader_Topology.ConnectionMatrix); %(N+1)*1/B转置

			obj.Extended_Topology.DegreeMatrix = diag(sum(obj.Extended_Topology.AdjacencyMatrix,2));

			obj.Extended_Topology.LaplacianMatrix = obj.Extended_Topology.DegreeMatrix - ...
				obj.Extended_Topology.AdjacencyMatrix; %L_tilde 
		end
		%% Set Control Gain
		function set_ControlGainSet(obj,ControlGainSet_in)
			obj.ControlGainSet = ControlGainSet_in;
			for AgentNr = 1:obj.AgentQuantity
				obj.AgentSet{AgentNr}.ControlGainSet = obj.ControlGainSet;
			end
		end
		%% Intitialize Leader Set
		function initialize_LeaderSet(obj,q_ref_set)
			for LeaderNr = 1:obj.LeaderQuantity
				obj.LeaderSet{LeaderNr} = Leader_Class( ...
					obj.q_dim,obj.SystemOrder,q_ref_set{LeaderNr});
			end
		end
		%% Local Computation
		function do_AgentsLocalComputation(obj)
			for AgentNr = 1:obj.AgentQuantity
				obj.AgentSet{AgentNr}.do_LocalComputation;
			end
		end
		%% Share Information
		function share_Information(obj)
			for AgentNr = 1:obj.AgentQuantity
				obj.share_AgentState_for_SingleAgent(AgentNr);
				obj.share_LeaderState_for_SingleAgent(AgentNr);
				obj.share_NeighbourPrediction_for_SingleAgent(AgentNr);
			end
			for AgentNr = 1:obj.AgentQuantity
				obj.AgentSet{AgentNr}.SharedInformation = [];
			end
		end
		%% Share Information
		function share_Information_2(obj)
			for AgentNr = 1:obj.AgentQuantity
				obj.share_NeighbourPrediction_for_SingleAgent(AgentNr)
			end
			for AgentNr = 1:obj.AgentQuantity
				obj.AgentSet{AgentNr}.SharedInformation = [];
			end
		end
		%% Share Information for Single Agent
% 		function share_Information_for_SingleAgent(obj,AgentNr)
% 			obj.share_AgentState_for_SingleAgent(AgentNr);
% 			obj.share_LeaderState_for_SingleAgent(AgentNr);
% 			obj.share_NeighbourPrediction_for_SingleAgent(AgentNr);
% 		end
		%% Share Agent States for Single Agent
		function share_AgentState_for_SingleAgent(obj,AgentNr)
			NeighbourNrSet = obj.Agent_Topology.NeighbourSet{AgentNr};
			NeighbourStateSet = cell(numel(NeighbourNrSet),1);
			%
			for NeighbourNr = numel(NeighbourNrSet):-1:1
				NeighbourAgentNr = NeighbourNrSet(NeighbourNr);
				NeighbourAgent = obj.AgentSet{NeighbourAgentNr}; 
				
				if isfield(NeighbourAgent.SharedInformation,'LocalState')
					NeighbourStateSet{NeighbourNr} = ...
						NeighbourAgent.SharedInformation.LocalState;
					NeighbourStateSet{NeighbourNr}.ConnectionWeight = ...
						obj.Agent_Topology.AdjacencyMatrix(AgentNr,NeighbourAgentNr);
				else
					
				end
			end
			obj.AgentSet{AgentNr}.ReceivedInformation.NeighbourStateSet = NeighbourStateSet;
		end
		%% Share Leader States for Single Agent
		function share_LeaderState_for_SingleAgent(obj,AgentNr)
			LeaderNeighbourNrSet = obj.Agent_Leader_Topology.AgentNeighbour{AgentNr};
			LeaderStateSet = cell(numel(LeaderNeighbourNrSet),1);
			for LeaderNeighbourNr = 1:numel(LeaderNeighbourNrSet)
				LeaderNr = LeaderNeighbourNrSet(LeaderNeighbourNr);
				Leader = obj.LeaderSet{LeaderNr};

				LeaderStateSet{LeaderNeighbourNr}.State = Leader.ReferenceTrajectoryFunction(obj.t);
				LeaderStateSet{LeaderNeighbourNr}.ConnectionWeight = ...
					obj.Agent_Leader_Topology.ConnectionMatrix(AgentNr,LeaderNr);
				LeaderStateSet{LeaderNeighbourNr}.StateDerivative = Leader.dqmdt_Function(obj.t);
			end
			obj.AgentSet{AgentNr}.ReceivedInformation.LeaderStateSet = LeaderStateSet;
		end
		%% Share Neighbour Prediction for Single Agent
		function share_NeighbourPrediction_for_SingleAgent(obj,AgentNr)
			NeighbourNrSet = obj.Agent_Topology.NeighbourSet{AgentNr};
			NeighbourPredictionSet = cell(numel(NeighbourNrSet),1);
			%
			for NeighbourNr = numel(NeighbourNrSet):-1:1
				NeighbourAgentNr = NeighbourNrSet(NeighbourNr);
				NeighbourAgent = obj.AgentSet{NeighbourAgentNr}; 
				
				if isfield(NeighbourAgent.SharedInformation,'NeighbourPrediction')
					NeighbourPrediction = NeighbourAgent.SharedInformation.NeighbourPrediction;
					NeighbourPredictionPosition = ...
						find(NeighbourPrediction.sequence == AgentNr,1);
					if ~isempty(NeighbourPredictionPosition)
						NeighbourPredictionSet{NeighbourNr} = ...
							NeighbourPrediction.prediction{NeighbourPredictionPosition};
					else
						NeighbourPredictionSet{NeighbourNr} = [];
					end

				else
					
				end
			end
			obj.AgentSet{AgentNr}.ReceivedInformation.NeighbourPredictionSet = NeighbourPredictionSet;

		end
		%% Cross Computation
		function do_AgentsCrossComputation(obj)
			for AgentNr = 1:obj.AgentQuantity
				obj.AgentSet{AgentNr}.do_CrossComputation;
			end
		end
		%% 
		function do_AgentsAggregationComputation(obj)
			for AgentNr = 1:obj.AgentQuantity
				obj.AgentSet{AgentNr}.do_AggregationComputation;
			end
		end
		%% get Control Input
		function set_ControlInput(obj)
			for AgentNr = 1:obj.AgentQuantity
				obj.AgentSet{AgentNr}.set_ControlInput;
			end
		end
		%% Envolove Multi Agent System
		function x_end = envolve_MultiAgentSystem(obj)
			x_set = nan(obj.q_dim * obj.AgentQuantity, obj.SystemOrder);
			for AgentNr = 1:obj.AgentQuantity
				obj.AgentSet{AgentNr}.envolve_AgentState;
				xi = obj.AgentSet{AgentNr}.State;
				x_set((AgentNr - 1) * obj.q_dim + (1:obj.q_dim), 1:obj.SystemOrder) = ...
					reshape(xi,[obj.q_dim,obj.SystemOrder]);
			end
			obj.x = reshape(x_set,[obj.x_dim * obj.AgentQuantity,1]);
			obj.t = obj.t + obj.t_step;

			x_end = obj.x;
		end
	end

end