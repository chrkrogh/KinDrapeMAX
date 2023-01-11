% Follow two steering curves in abs. coordinates
% The size of the course is 0.75 x 4 m , i.e. (Grid-1)*d

d = 0.05;				% Discretization distance		
Grid = [16 81];			% Grid size [nWidthNodes nLengthNodes]
Org = [inf inf];		% Origin point on the mold surface [xOrg yOrg]. Neither component is used due to the steering curves
Ang = inf;				% Ini. angle of the lengthwise geo generator. Not used due to lengthwise steering curve
OrgNodePct = [100 0];	% Org. node % of the grid dim. Defined instead of Org.
PreShear = inf;			% Pre-shear angle. Not used due to 2 x steering curves

DraContr = {'Steer','Steer'};									% Generator type, i.e. 2 x steering curve
SteerPtsRef = {'Abs','Abs'};									% Reference for the steer points
SteerPts{1,1} = [0.5 0.0]' + [0.0 0.0 ; 0.0 2.0 ; -0.2 4.5]';	% Points for lengthwise steer. curve
SteerPts{2,1} = [0.5 0.0]' + [0.0 0.0 ; -1.0 0.0]';				% Points for widthwise steer. curve 