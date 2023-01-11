% Example of abs. steering curves to control draping on hemisphere
d = 0.05;			% Discretization distance
Grid = [16 16];		% Grid size [nWidthNodes nLengthNodes]
Org = [inf inf];	% Origin point on the mold surface [xOrg yOrg]. Neither component is used due to the steering curves
Ang = inf;			% Ini. angle of the lengthwise geo generator. Not used due to lengthwise steering curve
OrgNode = [1 1];	% Org. node in the grid, i.e. lower left corner
PreShear = inf;		% Pre-shear angle. Not used due to 2 x steering curves

DraContr = {'Steer','Steer'};									% Generator type, i.e. 2 x steering curve
SteerPtsRef = {'Abs','Abs'};									% Reference for the steer points
SteerPts{1,1} = [0.0 0.0]' + [0.0 0.0 ;  0.1 0.3 ; 0.0 0.9]';	% Points for lengthwise steer. curve
SteerPts{2,1} = [0.0 0.0]' + [0.0 0.0 ; 0.9 -0.05]';			% Points for widthwise steer. curve 