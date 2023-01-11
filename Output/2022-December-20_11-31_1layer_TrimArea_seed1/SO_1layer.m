% Stack with courses following the right edge of the mold
nLayers = 5;
nLayersUse = 1;
nCourses = [5 5 5 5 5];
nCoursesUse = [5 5 5 5 5];
LayerPropagation = 'right-to-left';
d = repmat({0.05},5,5);
Grid(1,:) = {[8 81],[8 81],[8 81],[8 81],[8 81]};
Grid(2,:) = {[8 81],[8 81],[8 81],[8 81],[8 81]};
Grid(3,:) = {[8 81],[8 81],[8 81],[8 81],[8 81]};
Grid(4,:) = {[8 81],[8 81],[8 81],[8 81],[8 81]};
Grid(5,:) = {[8 81],[8 81],[8 81],[8 81],[8 81]};

Org(1,:) = {[inf 0],[inf 0],[inf 0],[inf 0],[inf 0]};
Org(2,:) = {[inf 0],[inf 0],[inf 0],[inf 0],[inf 0]};
Org(3,:) = {[inf 0],[inf 0],[inf 0],[inf 0],[inf 0]};
Org(4,:) = {[inf 0],[inf 0],[inf 0],[inf 0],[inf 0]};
Org(5,:) = {[inf 0],[inf 0],[inf 0],[inf 0],[inf 0]};

Ang = inf;

OrgNodePct = repmat({[100 0]},5,5);
OrgNodePct{1,1} = [0 0];
OrgNodePct{2,1} = [0 0];
OrgNodePct{3,1} = [0 0];
OrgNodePct{4,1} = [0 0];
OrgNodePct{5,1} = [0 0];

%OrgNode = {[9 1],[9 1],[9 1],[9 1]};
PreShear = repmat({0.0},5,5);

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

SteerPts = {SteerPts SteerPts SteerPts SteerPts SteerPts};

% Make 100 mm offset from edge per layer
% SteerPts{2}{1}(1,:) = -0.1;
% SteerPts{3}{1}(1,:) = -0.2;
% SteerPts{4}{1}(1,:) = -0.3;
% SteerPts{5}{1}(1,:) = -0.4;

%% Optimization

% Toggle design variables
DesVar.OrgNodeLong = true;
DesVar.TransOffset = 1:3;
DesVar.FirstTransOffset = 1:5;
DesVar.CourseWidth = true;

% Toggle minimization criteria
Obj.Shear = false;
Obj.AngDev = false;
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
BndAndCon.MinStaggerDist = []; %20e-3;
% Available course widths (grid nodes)
BndAndCon.AvailableCourseWidths = [8 9 10];
% Bounds on trans. offset
BndAndCon.TransOffset = [-12e-3 0e-3];
BndAndCon.FirstTransOffset = [0e-3 150e-3];

GASet = optimoptions('ga',...
    'Display','iter',...
    'MaxGenerations',100,...
	'MaxStallGenerations',40,...
	'PopulationSize', 100*nLayersUse,...
	'EliteCount', 5*nLayersUse,...
	'CrossoverFraction', 0.8,...
    'UseParallel',true,...
    'FunctionTolerance',1e-6,...
    'PlotFcn',{@gaplotbestf,@gaplotstopping},...
	'NonlinearConstraintAlgorithm','penalty'...
    );

% Stored optimization
Inp.x_opt = [];


Inp.x_opt = [35           43           72            6           66   0.00297662    0.0766963     0.127765     0.109184    0.0327323  -0.00133057 -0.000303854  -0.00199139  -0.00440982   -0.0025075  -0.00425797  -0.00720439  -0.00218593  -0.00334371  -0.00293212  -0.00138353   -0.0030055            3            3            1            3            2];