% long airway tapering analysis

%% Compute taper gradient path of all terminal branches
% Table with all carina-terminal branch data.
% running ComputeAirwayLobes() before this function will group output by lobe.
AllTaperResults = ComputeTaperAll(AQ);

%% Display taper results of a single gradient path
% plot taper grad results from carina to a given end-node.
figure
PlotTaperResults(AQ, AllTaperResults.terminalnode(1), 'inner')
figure
PlotTaperResults(AQ, AllTaperResults.terminalnode(1))

%% Construct single taper gradient path
% this demonstrates more hands on functions to process data as desired.
% list of terminal node, refer to airway graph.
terminalnodelist = ListTerminalNodes(AQ);
% get branch data for carina to terminal node for given node.
% note a positive taper gradient shows that the branch is tapering.
[logtaperrate, cum_arclength, cum_area, path] = ConstructTaperPath(AQ, terminalnodelist(1)); 

%% Display boxplot of tapervalues by lobe
figure
TaperBoxPlot(AQ, 'inner')
figure
TaperBoxPlot(AQ)
