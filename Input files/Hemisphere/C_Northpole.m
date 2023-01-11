% Draping with 2 x geodesics starting from north pole
d = 0.05;						% Discretization distance
Grid = [53 53];					% Grid size [nWidthNodes nLengthNodes]
Org = [0.0 0.0];				% Origin point on the mold surface [xOrg yOrg].
Ang = 0;						% Ini. angle of the lengthwise geo generator
OrgNodePct = [50 50];			% Org. node % in the grid, i.e. center 
PreShear = 0;					% Pre-shear angle

DraContr = {'Geo','Geo'};		% Generator type, i.e. 2 x geodesics
SteerPtsRef = {'Abs','Abs'};	% Reference for the steer points (not used due to geo.)
SteerPts = {[];[]};				% Steering points (not used due to geo.)