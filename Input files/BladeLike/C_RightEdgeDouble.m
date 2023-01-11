% Follow the right edge of the mold and the root end (2 x steering curve)
% The size of the course is 0.5 x 4 m, i.e. (Grid-1)*d

d = 0.1; 			% Discretization distance
Grid = [6 41]; 		% Grid size [nWidthNodes nLengthNodes]
Org = [inf inf];	% Origin point on the mold surface [xOrg yOrg]. Neither component is used due to the steering curves
Ang = inf;			% Ini. angle of the lengthwise geo generator. Not used due to lengthwise steering curve
OrgNode = [6 1];	% Org. node in the grid, i.e. lower right corner 
PreShear = inf;		% Pre-shear angle. Not used due to 2 x steering curves

DraContr = {'Steer','Steer'};		% Generator type, i.e. 2 x steering curve
SteerPtsRef = {'Right','Bottom'};	% Reference for the steer points 
SteerPts{1,1} = [0 ; 0];			% Points for lengthwise steer. curve (constant offset)
SteerPts{2,1} = [0 ; 0];			% Points for widthwise steer. curve (constant offset)