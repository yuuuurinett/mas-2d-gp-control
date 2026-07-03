%% Plot Data Set
% ET_MAS_GP_Leader_Plot_GPDataSet
InitialDataFigureObj = figure('Name','Data Set');
InitialDataAxesObj = axes(InitialDataFigureObj);
hold(InitialDataAxesObj,'on');
InitialDataLegend_set  = cell(LocalGP_Quantity,1);
for LocalGP_Nr = 1:LocalGP_Quantity
	plot(InitialDataAxesObj,LocalGP_set{LocalGP_Nr}.X(1,:),LocalGP_set{LocalGP_Nr}.X(2,:),'o');
	InitialDataLegend_set{LocalGP_Nr} = ['GP ', num2str(LocalGP_Nr)];
end
axis(InitialDataAxesObj,DomainScale * [-1,1,-1,1]);
legend(InitialDataAxesObj, InitialDataLegend_set);