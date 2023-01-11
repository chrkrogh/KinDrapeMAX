% Layer with 4 courses
% First course: follow the right edge of the mold
% Remaining courses: 5 mm gap
% The size of each course is 0.4 x 4 m, i.e. (Grid-1)*d

% Number of courses to be defined in the input file
nCourses = 5;										
% How many of the courses will be used for analysis
nCoursesUse = 4;
% Draping order of the courses
LayerPropagation = 'right-to-left';
% Discretization distance for each course (repeat {0.05} in a 1 x 5 cell array)
d = repmat({0.05},1,5);
% The grid dimensions for each course
Grid = {[9 81],[9 81],[9 81],[9 81],[9 81]};
% The origin for each course (1st components are inactive due to lengthwise steer curve)
Org = {[inf 0],[inf 0],[inf 0],[inf 0],[inf 0]};
% The origin node % of the grid (repeat {[100 0]} in a 1 x 5 cell array)
OrgNodePct = repmat({[100 0]},1,5); 
% Pre-shear for each course (repeat {0.0} in a 1 x 5 cell array)
PreShear = repmat({0.0},1,5);

% Draping control and ref. for steer pts
DraContr = {'Steer','Geo'};
SteerPtsRef = {'Right','Bottom'};

% Define the steer points using the aux variable Gap
Gap = 5e-3;
% Define 5 (relative) lengthwise steer pts for the 1st course, i.e. an offset of 0 
% along x at y-coordinates 0, 1, 2, 3 and 4. Actually one steering point is enough 
% for this constant offset but it is easy to try to manipulate the offsets
SteerPts{1,1} = [0 0 0 0 0 ; ...
                 0 1 2 3 4];
% Define the lengthwise steering points for the remaing courses, i.e. an offset 
% of -Gap along x at y-coordinates 0, 1, 2, 3 and 4. Again, one steering point is 
% enough for this constant offset
SteerPts{1,2} = [-Gap -Gap -Gap -Gap -Gap ; 0 1 2 3 4]; 
SteerPts{1,3} = [-Gap -Gap -Gap -Gap -Gap ; 0 1 2 3 4];
SteerPts{1,4} = [-Gap -Gap -Gap -Gap -Gap ; 0 1 2 3 4];
SteerPts{1,5} = [-Gap -Gap -Gap -Gap -Gap ; 0 1 2 3 4];
% Define the widthwise steering points (not used though unless DraContr{2} is set 
% to 'Steer')
SteerPts(2,1:5) = {[ 0 ; ...
                   5e-3]};