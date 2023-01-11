% Starting in the lower right corner and drape along two geodesic lines
% The size of the course is 0.5 x 4 m, i.e. (Grid-1)*d

d = 0.05;			% Discretization distance
Grid = [11 81];	    % Grid size [nWidthNodes nLengthNodes]
Org = [0.66 0.0]; 	% Origin point on the mold surface [xOrg yOrg]
Ang = 0;			% Ini. angle of the lengthwise geo generator
OrgNode = [11 1];	% Org. node in the grid, i.e. lower right corner 
PreShear = 0;		% Pre-shear angle (=> angle between generators is 90 deg)

DraContr = {'Geo','Geo'};			% Generator type, i.e. 2 x geodesic.
SteerPtsRef = {'Right','Bottom'};	% Refence for the steer points (not used here)
SteerPts{1,1} = [];					% Points for lengthwise steer. curve (not used)
SteerPts{2,1} = []; 				% Points for widthwise steer. curve (not used)