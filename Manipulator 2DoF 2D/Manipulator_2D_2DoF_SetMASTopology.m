function MultiAgentSystem = Manipulator_2D_2DoF_SetMASTopology( ...
	AgentQuantity,LeaderQuantity)
%% Agent Topology Set
AgentConnectionConfigurationSet = {
	{1,2,'Directed',1};
	{2,3,'Undirected',1};
	{3,4,'Undirected',1};
	{4,5,'Directed',1};
	{5,6,'Undirected',1};
	{1,6,'Undirected',1};
    % {1,2,'Directed',0.5}; 
    % {1,6,'Directed',0.5};
    % {2,3,'Directed',0.5}; 
    % {2,4,'Directed',0.5};
    % {3,2,'Directed',0.5}; 
    % {3,1,'Directed',0.5};
    % {4,3,'Directed',0.5}; 
    % {4,5,'Directed',0.5};
    % {5,6,'Directed',0.5}; 
    % {5,4,'Directed',0.5};
    % {6,5,'Directed',0.5}; 
    % {6,1,'Directed',0.5}
};
%% Agent Leader Topology Set
AgentLeaderConnectionConfigurationSet = {
	{'Leader',1,'Agent',1,1};
% 	{'Leader',1,'Agent',2,1};
% 	{'Leader',1,'Agent',3,1};
	{'Leader',1,'Agent',4,1};
% 	{'Leader',1,'Agent',5,1};
% 	{'Leader',1,'Agent',6,1};
};
%% Create Multi-Agent System Object
MultiAgentSystem = MultiAgentSystem_Class(AgentQuantity,LeaderQuantity);
%% Agent Topology
MultiAgentSystem.Agent_Topology.addConnectionSet(AgentConnectionConfigurationSet);
% MultiAgentSystem.Agent_Topology.check_Connectivity(true);
%% Agent Leader Tolology
MultiAgentSystem.Agent_Leader_Topology.addConnectionSet( ...
	AgentLeaderConnectionConfigurationSet);
% 
MultiAgentSystem.get_ExtendedTopology;
%% Show: Topology
%MultiAgentSystem.Agent_Topology.show_topology;

end