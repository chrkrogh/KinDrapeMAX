function [OrgNode,OrgNodePct] = OrgNodeFromPct(OrgNode,OrgNodePct,Grid)
% This function computes the OrgNode based on the OrgNodePct (origin node
% defined as the percentages of the grid dimensions) and the Grid or the
% OrgNodePct based on the OrgNode and the Grid
if ~isempty(OrgNodePct)
    OrgNode(1) = max(1,min(Grid(1),ceil(Grid(1)*OrgNodePct(1)/100)));
    OrgNode(2) = max(1,min(Grid(2),ceil(Grid(2)*OrgNodePct(2)/100)));
else
    OrgNodePct = (100./(Grid-1)).*OrgNode-(100./(Grid-1));
end
end