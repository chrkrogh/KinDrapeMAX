% Stack with courses following the right edge of the mold
% The steering curve of the 1st course in each layer follows the left course edge
% so that it always remains inside the net boundary, which makes smoother results
% For this reason the origin nodes and steer points of these courses are adapted
nLayers = 30;
nLayersUse = 30;
nCourses = repmat(5,1,nLayers);
nCoursesUse = repmat(5,1,nLayers);
LayerPropagation = 'right-to-left';
d = repmat({0.05},nLayers,5);
Grid = repmat({[8 81],[8 81],[8 81],[8 81],[8 81]},nLayers,1);
Org = repmat({[inf 0],[inf 0],[inf 0],[inf 0],[inf 0]},nLayers,1);

OrgNodePct = repmat({[100 0]},nLayers,5);
% Set the origin nodes of the 1st course in each layer to [0 0]
OrgNodePct(:,1) = repmat({[0 0]},nLayers,1);

PreShear = repmat({0.0},nLayers,5);

DraContr = {'Steer','Geo'};
SteerPtsRef = {'Right','Bottom'};

Gap = 5e-3;
SteerPts{1,1} = [-0.35 -0.35 -0.35 -0.35 -0.35 ; ... 					 		
                  0  1  2  3  4];
SteerPts{1,2} = [-Gap -Gap -Gap ; 0 2 4]; 
SteerPts{1,3} = [-Gap -Gap -Gap ; 0 2 4];
SteerPts{1,4} = [-Gap -Gap -Gap ; 0 2 4];
SteerPts{1,5} = [-Gap -Gap -Gap ; 0 2 4];

SteerPts(2,1:5) = {[ 0 ; ...
                   5e-3]};

SteerPts = repmat({SteerPts},1,nLayers);

%% Optimization

% Toggle design variables
DesVar.OrgNodeLong = true;
DesVar.TransOffset = 1:3;
DesVar.FirstTransOffset = 1:5;
DesVar.CourseWidth = true;

% Toggle minimization criteria
Obj.Shear = true;
Obj.AngDev = true;
Obj.TrimArea = true;

% Objective function settings
Obj.p_Shear = inf;
Obj.p_AngDev = 6;
Obj.MeanOfLayerPNorms = true;
Obj.TargetFiberAng = 0;
Obj.ShearAng_eq = 6;
Obj.TrimArea_eq = 0.25*(9-1)*(81-1)*0.05^2;

% Bounds and constraint limits
% Stagger distance constraint (empty for N/A)
BndAndCon.MinStaggerDist = 20e-3;
% Available course widths (grid nodes)
BndAndCon.AvailableCourseWidths = [8 9 10];
% Bounds on trans. offset
BndAndCon.TransOffset = [-12e-3 0e-3];
BndAndCon.FirstTransOffset = [0e-3 150e-3];

% GA settings
% The MaxGenerations, MaxStallGenerations, PopulationSize and EliteCount are 
% applied per layer in use
GASet = optimoptions('ga',...
    'Display','iter',...
    'MaxGenerations',100,...
	'MaxStallGenerations',40,...
	'PopulationSize', 100,...
	'EliteCount', 5,...
	'CrossoverFraction', 0.8,...
    'UseParallel',true,...
    'FunctionTolerance',1e-6,...
    'PlotFcn',{@gaplotbestf,@gaplotstopping},...
	'NonlinearConstraintAlgorithm','penalty'...
    );

% Stored optimization
Inp.x_opt = [];
