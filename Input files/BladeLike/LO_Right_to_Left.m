% Follow the right edge of the mold
nCourses = 4;
nCoursesUse = 4;
LayerPropagation = 'right-to-left';
d = repmat({0.05},1,4);
Grid = {[9 81],[9 81],[9 81],[9 81]};
Org = {[inf 0],[inf 0],[inf 0],[inf 0]};
Ang = inf;
OrgNodePct = repmat({[100 0]},1,4);
PreShear = repmat({0.0},1,4);

DraContr = {'Steer','Geo'};
SteerPtsRef = {'Right','Bottom'};

Gap = 5e-3;
SteerPts{1,1} = [0 0 0 0 0 ; ...
                 0 1 2 3 4];
SteerPts{1,2} = [-Gap -Gap -Gap -Gap -Gap ; 0 1 2 3 4]; 
SteerPts{1,3} = [-Gap -Gap -Gap -Gap -Gap ; 0 1 2 3 4];
SteerPts{1,4} = [-Gap -Gap -Gap -Gap -Gap ; 0 1 2 3 4];

SteerPts(2,1:4) = {[ 0 ; ...
                   5e-3]};

%% Optimization
% Toggle design variables
DesVar.OrgNodeLong = true;
DesVar.FirstTransOffset = false;
DesVar.TransOffset = 3:5;
DesVar.CourseWidth = false;

% Toggle minimization criteria
Obj.Shear = true;
Obj.AngDev = true;
Obj.TrimArea = false;

% Objective function settings
Obj.p_Shear = 6;
Obj.p_AngDev = inf;
Obj.MeanOfLayerPNorms = true;
Obj.TargetFiberAng = 0.0;
Obj.ShearAng_eq = 6;
Obj.TrimArea_eq = 0.4;

% Bounds and constraints
% Set trans offsets
BndAndCon.TransOffset = [-40e-3 40e-3];
BndAndCon.FirstTransOffset = [0 10e-3];
% Set the course widths (grid nodes)
BndAndCon.AvailableCourseWidths = [7 9 11];
% Set stagger dist (empty: N/A)
BndAndCon.MinStaggerDist = [];

% GA settings
GASet = optimoptions('ga',...
    'Display','iter',...
    'MaxGenerations',1,... 
	'MaxStallGenerations',40,...
	'PopulationSize', 50,...
	'EliteCount', 2,...
	'CrossoverFraction', 0.8,...
    'UseParallel',true,...
    'FunctionTolerance',1e-6,...
    'PlotFcn',{@gaplotbestf,@gaplotstopping},...
	'NonlinearConstraintAlgorithm','penalty'...
    );