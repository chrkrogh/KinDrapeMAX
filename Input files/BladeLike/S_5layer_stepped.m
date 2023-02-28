% Stack with courses following the right edge of the mold. For each layer, an inset
% of 100 mm is added.

% The number of layers and the number of layers to use
nLayers = 5;
nLayersUse = 5;
% The number of courses in each layer and the number of courses to use
nCourses = [4 4 4 4 4];
nCoursesUse = [4 4 4 3 3];
% Draping order of the courses
LayerPropagation = 'right-to-left';
% Discretization distance for each course (repeat {0.05} in a 5 x 4 cell array)
d = repmat({0.05},5,4);
% The grid dimensions for each course (nLayers x nCourses)
Grid(1,:) = {[9 81],[9 81],[9 81],[9 81]};
Grid(2,:) = {[9 81],[9 81],[9 81],[9 81]};
Grid(3,:) = {[9 81],[9 81],[9 81],[9 81]};
Grid(4,:) = {[9 81],[9 81],[9 81],[9 81]};
Grid(5,:) = {[9 81],[9 81],[9 81],[9 81]};
% The origin for each course (1st components are inactive due to lengthwise steer curve)
Org(1,:) = {[inf 0],[inf 0],[inf 0],[inf 0]};
Org(2,:) = {[inf 0],[inf 0],[inf 0],[inf 0]};
Org(3,:) = {[inf 0],[inf 0],[inf 0],[inf 0]};
Org(4,:) = {[inf 0],[inf 0],[inf 0],[inf 0]};
Org(5,:) = {[inf 0],[inf 0],[inf 0],[inf 0]};
% The origin node % of the grid (repeat {[100 0]} in a 5 x 4 cell array)
OrgNodePct = repmat({[100 0]},5,4);
% Pre-shear for each course (repeat {0.0} in a 5 x 4 cell array)
PreShear = repmat({0.0},5,4);

% Draping control and ref. for steer pts
DraContr = {'Steer','Geo'};
SteerPtsRef = {'Right','Bottom'};

% Define the steer points using the aux variable Gap
Gap = 5e-3;
% Define SteerPts for a course
SteerPts_l{1,1} = [-0.01 -0.01 -0.01 ; ...
                  0  	2     4];
SteerPts_l{1,2} = [-Gap -Gap -Gap ; 0 2 4]; 
SteerPts_l{1,3} = [-Gap -Gap -Gap ; 0 2 4];
SteerPts_l{1,4} = [-Gap -Gap -Gap ; 0 2 4];

SteerPts_l(2,1:4) = {[ 0 ; ...
                   5e-3]};
% Repeat the course steer pt definition for each layer
SteerPts = {SteerPts_l SteerPts_l SteerPts_l SteerPts_l SteerPts_l};

% Make 100 mm offset from edge per layer
SteerPts{2}{1}(1,:) = -0.1;
SteerPts{3}{1}(1,:) = -0.2;
SteerPts{4}{1}(1,:) = -0.3;
SteerPts{5}{1}(1,:) = -0.4;
