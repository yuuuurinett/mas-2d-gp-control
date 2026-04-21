classdef Topology_Class < handle

	properties
		AgentQuantity;
		VertexSet;
		EdgeSet = [];
		NeighbourSet;
		AdjacencyMatrix;
		DegreeMatrix;
		LaplacianMatrix;
	end

	methods
		function obj = Topology_Class(AgentQuantity_in)
			obj.AgentQuantity = AgentQuantity_in;
			obj.VertexSet = 1:obj.AgentQuantity;
			obj.NeighbourSet = cell(obj.AgentQuantity,1);
			obj.AdjacencyMatrix = zeros(obj.AgentQuantity);
			obj.DegreeMatrix = zeros(obj.AgentQuantity);
			obj.LaplacianMatrix = obj.DegreeMatrix - obj.AdjacencyMatrix;
		end
		%% Add Connection
		function addConnection(obj,NodeNr1,NodeNr2,ConeectionType,ConnectionWeight)
			if nargin < 5 %如果调用时没有传第5个参数（即没有指定权重),就默认设为1
				ConnectionWeight = 1;
			end
			% Edge from Node 1 to Node 2
			if NodeNr1 > obj.AgentQuantity || NodeNr2 > obj.AgentQuantity
				fprintf('Add Connection Fail! The number of Node exceed the upper bound!\n');
				return;
			end
			switch ConeectionType
				case 'Undirected'
					addDirectedConnection(obj,NodeNr1,NodeNr2,ConnectionWeight)
					addDirectedConnection(obj,NodeNr2,NodeNr1,ConnectionWeight)
				case 'Directed'
					addDirectedConnection(obj,NodeNr1,NodeNr2,ConnectionWeight);
				otherwise
					fprintf('Wrong Connection Type!\n');
			end
			obj.LaplacianMatrix = obj.DegreeMatrix - obj.AdjacencyMatrix;
		end
		%% Add Connection Set
		function addConnectionSet(obj,ConnectionConfigurationSet)
			for ConnectionNr = 1:numel(ConnectionConfigurationSet)
				ConnectionConfiguration = ConnectionConfigurationSet{ConnectionNr};
				NodeNr1 = ConnectionConfiguration{1};
				NodeNr2 = ConnectionConfiguration{2};
				ConeectionType = ConnectionConfiguration{3};
				weight = ConnectionConfiguration{4};
				obj.addConnection(NodeNr1,NodeNr2,ConeectionType,weight);
			end
		end
		%% Add Directed Connection
		function addDirectedConnection(obj,NodeNr1,NodeNr2,ConnectionWeight)
			if isempty(obj.EdgeSet)
				obj.EdgeSet = [obj.EdgeSet;[NodeNr1,NodeNr2]];
				obj.NeighbourSet{NodeNr2} = [obj.NeighbourSet{NodeNr2},NodeNr1];
			else
				ExistedEdgeNr = find(obj.EdgeSet(:,1) == NodeNr1 & obj.EdgeSet(:,2) == NodeNr2, 1);
				if isempty(ExistedEdgeNr)
					obj.EdgeSet = [obj.EdgeSet;[NodeNr1,NodeNr2]];
					obj.NeighbourSet{NodeNr2} = [obj.NeighbourSet{NodeNr2},NodeNr1];
				end
			end
			obj.DegreeMatrix(NodeNr2,NodeNr2) = ...
				obj.DegreeMatrix(NodeNr2,NodeNr2) - obj.AdjacencyMatrix(NodeNr2,NodeNr1) + ConnectionWeight;
			obj.AdjacencyMatrix(NodeNr2,NodeNr1) = ConnectionWeight;
		end
		%% Show Topology
		function show_topology(obj)
			TopologyFigureObj = figure('Name','Topology');
			TopologyAxesObj = axes(TopologyFigureObj);
			hold(TopologyAxesObj,'on');
			% Node and Label
			r_set = nan(2,obj.AgentQuantity);
			MarkerRadius = 0.1;
			MarkerSmoothFactor = 20;
			CircleAngleSet = linspace(0,2*pi,MarkerSmoothFactor);
			CircleData = MarkerRadius * [cos(CircleAngleSet);sin(CircleAngleSet)];
			for AgentNr = 1:obj.AgentQuantity
				phi = -(AgentNr - 1) * 2 * pi / obj.AgentQuantity;
				r_set(:,AgentNr) = 1 * [sin(phi);cos(phi)];

				
				plot(TopologyAxesObj, ...
					CircleData(1,:) + r_set(1,AgentNr), ...
					CircleData(2,:) + r_set(2,AgentNr),'r-');

				r_label = [sin(phi);cos(phi)] + [-0.015;0];
				text(TopologyAxesObj,r_label(1),r_label(2),num2str(AgentNr));
			end
			% Edge
			for EdgeNr = 1:size(obj.EdgeSet,1)
				EdgeNodeNr1 = obj.EdgeSet(EdgeNr,1);
				EdgeNodeNr2 = obj.EdgeSet(EdgeNr,2);
				NodePosition1 = r_set(:,EdgeNodeNr1); %2*N:每列是一个节点
				NodePosition2 = r_set(:,EdgeNodeNr2);
				DeltaPosition = NodePosition2 - NodePosition1;
				DeltaDirection = DeltaPosition / norm(DeltaPosition);

				ArrowStartPosition = NodePosition1 + ...
					DeltaDirection * MarkerRadius;
% 				ArrowEndPosition = NodePosition2 - ...
% 					DeltaDirection * MarkerRadius;
% 				ArrowLine = [ArrowStartPosition,ArrowEndPosition];
% 				plot(TopologyAxesObj,ArrowLine(1,:),ArrowLine(2,:),'b-');

				ArrowDeltaPosition = DeltaDirection * (norm(DeltaPosition) - 2 * MarkerRadius);
				quiver(TopologyAxesObj, ...
					ArrowStartPosition(1),ArrowStartPosition(2),ArrowDeltaPosition(1),ArrowDeltaPosition(2),0, ...
					'color','b','AutoScale','off');
			end
			%
			axis(TopologyAxesObj,[-1.2,1.2,-1.2,1.2]);
		end
		function show_topology_3D(obj,TopologyHeight)
			if nargin < 2
				TopologyHeight = 0;
			end
			TopologyFigureObj = figure('Name','Topology');
			TopologyAxesObj = axes(TopologyFigureObj);
			hold(TopologyAxesObj,'on');
			% Node and Label
			r_set = nan(3,obj.AgentQuantity);
			r_set(3,:) = TopologyHeight;
			MarkerRadius = 0.1;
			MarkerSmoothFactor = 20;
			CircleAngleSet = linspace(0,2*pi,MarkerSmoothFactor);
			CircleData = MarkerRadius * [cos(CircleAngleSet);sin(CircleAngleSet)];
			for AgentNr = 1:obj.AgentQuantity
				phi = -(AgentNr - 1) * 2 * pi / obj.AgentQuantity; %让N个节点均匀分布在单位圆上
				r_set(:,AgentNr) = 1 * [sin(phi);cos(phi)];

				
				plot3(TopologyAxesObj, ...
					CircleData(1,:) + r_set(1,AgentNr), ...
					CircleData(2,:) + r_set(2,AgentNr), ...
					r_set(3,:),'r-');

				r_label = [sin(phi);cos(phi)] + [-0.015;0];
				text(TopologyAxesObj,r_label(1),r_label(2),TopologyHeight,num2str(AgentNr));
			end
			% Edge
			for EdgeNr = 1:size(obj.EdgeSet,1)
				EdgeNodeNr1 = obj.EdgeSet(EdgeNr,1);
				EdgeNodeNr2 = obj.EdgeSet(EdgeNr,2);
				NodePosition1 = r_set(:,EdgeNodeNr1);
				NodePosition2 = r_set(:,EdgeNodeNr2);
				DeltaPosition = NodePosition2 - NodePosition1;
				DeltaDirection = DeltaPosition / norm(DeltaPosition);

				ArrowStartPosition = NodePosition1 + ...
					DeltaDirection * MarkerRadius;
				ArrowDeltaPosition = DeltaDirection * (norm(DeltaPosition) - 2 * MarkerRadius);
				quiver(TopologyAxesObj, ...
					ArrowStartPosition(1),ArrowStartPosition(2),ArrowDeltaPosition(1),ArrowDeltaPosition(2),0, ...
					'color','b','AutoScale','off');
			end
			%
			axis(TopologyAxesObj,[-1.2,1.2,-1.2,1.2]);
		end
		%% Check Connectivity
		function isConnected = check_Connectivity(obj,doShowResult)
			if nargin <= 1
				doShowResult = true;
			end
			if rank(obj.LaplacianMatrix) == obj.AgentQuantity - 1 %L*1=0=>奇异矩阵=>降秩
				isConnected = true;
				if doShowResult
					fprintf('Graph is connected!\n');
				end
			else
				isConnected = false;
			end
		end
	end
end
