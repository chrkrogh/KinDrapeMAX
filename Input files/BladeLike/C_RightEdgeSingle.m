% Lengthwise/longitudinal generator: offset from the right edge of the mold 
% Widthwise / transverse generator: geodesic line starting at y = 0.5, i.e. Org(2)
% The size of the course is 0.5 x 4 m, i.e. (Grid-1)*d

d = 0.05;			% Discretization distance
Grid = [11 81];		% Grid size [nWidthNodes nLengthNodes]
Org = [inf 0.5];	% Origin point on the mold surface [xOrg yOrg]. 1st component is not used due to the steering curve
Ang = inf;			% Ini. angle of the lengthwise geo generator. Not used due to lengthwise steering curve
OrgNode = [11 10];	% Org. node in the grid, i.e. 10th node on right course edge
PreShear = 0;		% Pre-shear angle (=> angle between generators is 90 deg)

DraContr = {'Steer','Geo'};			% Generator type, i.e. lengthwise steering curve + widthwise geodesic
SteerPtsRef = {'Right','Bottom'};	% Refence for the steer points (2nd cell is not used here)
SteerPts{1,1} = [0 -0.1 0 ; 0 2 4];	% Points for lengthwise steer. curve
SteerPts{2,1} = [];					% Points for widthwise steer. curve (not used)