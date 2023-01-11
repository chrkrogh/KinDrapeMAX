´´´
   _   ___      ______                     ___  ___  ___  __   __
  | | / (_)     |  _  \                    |  \/  | / _ \ \ \ / /
  | |/ / _ _ __ | | | |_ __ __ _ _ __   ___| .  . |/ /_\ \ \ V / 
  |    \| | '_ \| | | | '__/ _` | '_ \ / _ \ |\/| ||  _  | /   \ 
  | |\  \ | | | | |/ /| | | (_| | |_) |  __/ |  | || | | |/ /^\ \
  \_| \_/_|_| |_|___/ |_|  \__,_| .__/ \___\_|  |_/\_| |_/\/   \/
                                | |                              
			       |_|  

´´´

# User guide

## Introduction
The program can analyze and optimize the draping pattern on a mold using a kinematic draping model and the built-in MATLAB implementation of genetic algorithm, i.e. the function ga. It is intended for wind turbine blade production, i.e. where multiple roll-widths or courses of fabric are necessary to cover the mold and where the mold is more or less rectangular. It can, however, also analyze a single ply on any smooth surface. 

The program was used for generating the results of the journal paper:
[insert reference when published]
The idea with particularly the optimization framework is well documented in the paper.

The program is executed from the script KinDrapeMAX.m. The program reads the input file and the mold data, analyzes or optimizes the draping pattern (depending on the user settings), plots results in figures, plots results in the command window and lastly asks the user if results should be stored. 

When a mold data file is selected in the beginning of the script (Inp.MoldName, see description of settings below), only the input files pertaining to this mold are available.

The program uses mold surfaces of the form z = F(x,y) where F is a function or interpolating object. That is, for each x,y-coordinate there is a unique z-coordinate. For this reason, many variables are defined only with x,y-coordinates and the z-coordinates is retrieved by projection onto the mold surface. 
At a 0 degree angle, rows are in the (positive) x-direction and columns are in the (positive) y-direction.
The (close to) rectangular molds used for wind turbine blade production have the width direction alligned with the x-axis and the length direction alligned with the y-direction. z is positive upwards (away from the mold surface, in the stacking direction).

N.B. the terms lengthwise/longitudinal and widthwise/transverse will be used interchangeably.

## Main input/output
Different structs are used to store input and output to the program:
- Inp: input (including the contents from the input file)
- Set: program settings
- Plt: plot settings and plot data
- Opt: design variable bookkeeping and bounds (plus temporary opt. data)
- OptProb: optimization problem structure that can be input to MATLAB's ga function
- Dra_S: Draping data for the stack
- OptRes_S: Optimization results for the stack

The mold data is stored in different cell arrays (one cell for each layer, i.e. base surface and mold offsets due to the stack build-up)
- Mold: X,Y,Z grid matrices for surface plotting (inside net boundary) plus mold net boundary lines (MoldEdge)
- F: surface interpolation object created with scatteredInterpolant (query any point incl. extrapolation)
- DT: 2D Delaunay triangulation object describing the input triangular mesh in the XY-plane.
- z: z-coordinates pertaining to each 2D node of DT.

The latter two are used for the geodesic curve algorithm. 

## The draping algorithm 
As an introduction to the basic draping algorithm, please see the following educational arcticle:

Krogh, C., Bak, B.L.V., Lindgaard, E. et al. A simple MATLAB draping code for fiber-reinforced composites with application to optimization of manufacturing process parameters. Struct Multidisc Optim 64, 457–471 (2021). https://doi.org/10.1007/s00158-021-02925-z. Free read-only access: https://rdcu.be/c22yN

The implementation described here builds upon the simple code.
The building block of the draping algorithm is a single course, i.e. a rectangular ply stemming from a fabric roll of fixed width. The course dimensions are defined from a discretization distance, d, and a grid size, i.e. the number of nodes in the width and length direction. 
The placement of the course is defined through two "generators", i.e. curves on the mold surface that is populated with nodes and define respectively a row and column of the draping grid. The location on the mold where the generators intersect it the origin. The corresponding node in the draping grid is the origin node.
The generators can either be defined from a steering curve, i.e. a set of points or calculated as a geodesic curve. They can also be mixed. E.g. a steering curve in the length direction can be paired with a geodesic line in the width direction. In this case, the angle between the two generators is defined as the preshear angle.
Furthermore, steering curves can be defined in absolute coordinates on the mold or be defined relative to a net boundary line. That is, the lengthwise or longutidinal generator can e.g. be defined relative to the right mold net boundary line and the widthwise or transverse generator can be defined relative to the bottom/root net boundary line.

When draping a layer, the first course can be draped as described above for a single course, and the remaining courses are draped relative to the previous course. That is, steering points defined for courses #2 to #Ncourses are relative to the edge of the previous course. When draping a stack, each new layer is draped on an offset mold surface.

## Optimization
The optimizer can manipulate the course placement (offsets from net boundary line or other courses), the shear distribution (longituidnal/lengthwise origin node index, i.e. with the lengthwise steering curve) and width, while minimizing criteria such as shear angles, deviations of the UD fiber angles from a nominal value and the trim area/material waste. 

## Settings in the main script
This is an overview of the settings and tweaks in the main script. First general settings:
- rng(\<arg\>): seed for the random number generator used for the genetic algorithm (\<arg\> must be an integer).
- Inp.MoldName: the mold used for draping (string). It must be a MATLAB function located in a subfolder in the "Molds"-folder. Both function and subfolder must have the same name. See the guide in the "Molds"-folder.
- Inp.FileName: [optional] input file name (string). Must correspond to a MATLAB script in a subfolder in the "Input files"-folder. Both script and subfolder must have the same name analogous to the mold system. If Inp.FileName is not provided or it is empty, the program displays the available input files in the command window and the user can choose one.
- Set.Mode: choose how the program runs. The options are (string):
	- 'analysis': run the input file for analysis (i.e. without optimization)
	- 'baseline': run the input file for analysis but through the framework of the optimization (to check) 
	- 'opt': run the input file for optimization (all layers at the same time, "all-layer" formulation)
	- 'opt-seq': run the input file for optimization (one layer at the time, "sequential" formulation)
	- 'storedopt': run an input file with a stored vector of design variables (stored by the program)
	- 'error': run a vector of design variables that was stored during a program error to reproduce it. This storing is active if the variable Set.OptDebug is set to true.
- Set.InpFileType: filter for selecting input file type. The options are (string):
	- 'course': C_*
	- 'layer' : L_*
	- 'stack' : S_*
	- 'layer-optimization': LO_*
	- 'stack-optimization': SO_*

This naming convention must be followed by the input files. I.e. an input file for a single course must start with "C_". The reason why different input file types are defined is to make appropriate checks on data types and values when loading the input file. A single course only needs little information in comparison to a stack for optimization and it is thus a way of making a simpler input file that can still be checked. Notice though that e.g. a stack-optimization input file also can be used to analyze a single course but it must contain all the required information about the entire stack as well as the optimization settings.

Naturally, Set.Mode must match the input file type. If e.g. an optimization input file is not loaded ('layer-optimization' or 'stack-optimization') the only available option is 'analysis'.

Other settings:
- Set.TrimPlyToNetBd: Trim the draped courses to the net boundary (true/false). Must be active if trim area is a criterion in the optimization.
- Set.ExtrapSteerCurves: extrapolate the steering curves input by the user, if they are not long enough for the specified course length (true/false).
- Set.MeshOrCurveWarning: issue warnings when steering curves or geodesic lines are extrapolated or adjusted (true/false).
- Set.ForceGeoLineIniAngle: make the geodesic line straight for one cell discretization such that the first cell has an angle equal to PreShear (true/false).
- Set.OptDebug: store the vector of design variables x_opt if the program encounters an error during optimization (true/false). It is not very easy to debug inside ga...
- Set.Timing: time or profile the code (string: 'time' or 'profile').
- Set.CheckAspect: check the aspect ratio of the triangular mesh (true/false). Distorted triangles can be a problem for the geodesic line algorithm.
- Set.AspectLim: limit on the triangle aspect ratio before a warning is issued (scalar > 0). 

Plot settings:
- Plt.PlyPlt: quantities to plot after draping analysis/optimization (cell array of strings with options: 'shear' for shear angles, 'fibdev' for UD fiber angle deviations and 'flat' for a single color for each course).
- Plt.NodesWithoutTrim: plot all the nodes of the draping grid without trimming operations (true/false). 
- Plt.AbsShearAng: plot the absolute values of the shear angles (true/false).
- Plt.AbaqusColor: use the same three-color (blue/yellow/red) theme as Composites Modeler for Abaqus CAE instead of MATLAB's jet color theme (true/false).
- Plt.ShearLim: Shear limit for Abaqus-style plot (< 50 % of limit: blue, > limit: red, otherwise: yellow) (scalar).
- Plt.TriMesh: plot the triangular mesh (true/false).
- Plt.MoldNetBd: plot the mold net boundary lines from Mold.MoldEdge (true/false).
- Plt.GenCurve: plot the generator curves for each course (true/false).
- Plt.GeoInc: plot the incremental operations in the geodesic line algorithm (true/false). N.B. generates a plot figure per geodesic line in the model.
- Plt.MeshBoundary: plot the mesh boundary, i.e. extended boundary (true/false).
- Plt.DispLegend: display a plot legend (true/false).
- Plt.zShift: shift the plies by this amount in z when plotting. Trick to avoid the ply and mold colliding (scalar).
- Plt.FigWindowSize: specify the figure window size (string with options: 'normal' for MATLAB standard, 'full' for full-screen, 'half' for approximately left half or 'full-other' for full screen on another monitor (define this by first running the aux. script "DefineOtherScreen.m".
- Plt.Display: display output in figures (true/false).

## Input files
While reading the following description of the different input files it is recommended to simultaneously go through the example input files in the "Input files"-folder.

For a single course the following variables must be defined in the input file (filename beginning with "C_"):
- d: discretization (scalar > 0)
- Grid: size of the draping grid ([nWidthNodes nLengthNodes] with each entry being an integer > 0)
- Org: origin point on the mold for geodesics ([xOrg yOrg] with each entry being a scalar, naturally on the mold surface). If a steering curve is used as the lengthwise generator, the x-component of Org becomes redundant and if in addition the widthwise generator is a steering curve, also the y-component of Org becomes redundant. The variable must however always be defined.  
- Ang: the initial angle of a lengthwise geodesic generator (scalar between -90 and 90). If a steering curve is used as the lengthwise generator, the value is not used, but the variable must still be defined.
- OrgNode: the origin node in the draping grid ([WidthOrg LengthOrg] which must be within the draping grid). This variable can also be replaced by the variable OrgNodePct which enables to write the org. node in % of the grid dim ([%OrgWidth %OrgLength] where each component is a scalar between 0 and 100).
- PreShear: pre shear angle between the lengthwise generator (steering curve or geodesic) and a widthwise geodesic generator (scalar < 90 deg). A pre shear angle of 0 means that the two generators will start 90 deg apart. If a steering curve is used as the widthwise generator, the value it not used, but the variable must still be defined.
- DraContr: draping control, i.e. generator type (cell array with a cell for lengthwise and widthwise. Options: 'Steer' for steering curve and 'Geo' for geodesic curve).
- SteerPtsRef: reference for the steering points (cell array. Options for 1st cell (lengthwise): 'Right', 'Left' or 'Abs'. Options for 2nd cell (widthwise): 'Bottom' or 'Abs'. Here 'Abs' means that the steering points are in absolute mold coordinates, i.e. not relative to an edge). The variable is only active if 'Steer' is chosen for any of the DraContr cells but it must always be defined.
- SteerPts: definition of steering points (2 x 1 cell array where the first cell pertains to the lengthwise generator and the second cell pertains to the widthwise generator. Each cell is a numeric array with two rows where the first is the x-coordinates and the second is the y-coordinates of the steering points. There can be as many columns, i.e. steering points as needed. If e.g. relative steering points are used, a single point can be used to specify a constant offset).

For a layer input file (filename beginning with "L_"), some additional variables must be specified:
- nCourses: number of courses defined in the layer input file (integer > 0)
- nCoursesUse [optional]: the number of courses of the ones defined to use (integer > 0 and <= nCourses)
- LayerPropagation: draping order of the courses (string with options: 'right-to-left' or 'left-to-right').
Further, the variables d, Grid, Org, OrgNode/OrgNodePct and PreShear must now be expanded to cell arrays of size 1 x nCourses where each cell is the variable defined for a single course. E.g. Grid for a three course layer will look like: 
Grid = {[nRowNodes1 nColNodes1],[nRowNodes2 nColNodes2],[nRowNodes1 nColNodes3]};
The variable Ang is not active and therefore not necessary to define for layers because steering curves must be used as lengthwise generators.
The variables DraContr and SteerPtsRef are unchanged because these settings apply to all courses in a layer.
The variable SteerPts must now be a 2 x nCourses cell array where the 1st row is still pertaining to the lengthwise steering points and the 2nd row pertaining to the widthwise steering points. Each column now pertains to a course in the layer.

For a stack input file (filename beginning with "S_"), some additional variables must be specified:
- nLayers: number of layers defined in the input file (integer > 0).
- nLayersUse [optional]: the number of layers of the ones defined to use (integer > 0 and <= nLayers)
Further, the variables nCourses and nCoursesUse must now be specified as vectors with the number of elements equal to the number of layers. That is, each entry is the number of courses in the corresponding layer.
Further, the cell arrays of the variables d, Grid, Org, OrgNode/OrgNodePct and PreShear must now be expanded to size nLayers x nCourses. That is, each row of the cell arrays represent a layer as previously defined for a single layer.
The variables DraContr and SteerPtsRef are unchanged because these settings apply to all courses in a stack.
The variable SteerPts must now be a cell array of 1 x nLayers where each cell is a sub cell array of 2 x nCourses as described for a single layer.

For an optimization input file for a layer or a stack (filename beginning with "LO_" or "SO_") the paramterization, criteria and optimization settings must be specified in additional variables. The settings relating to design variables are stored in the struct DesVar:
- DesVar.OrgNodeLong: toggle the longitudinal/lengthwise origin node index of each course as a design variable (true/false).
- DesVar.TransOffset: toggle the transverse/widthwise offset between adjacent courses as design variable. Here, the optimizer will manipulate the x-coordinates of the already defined steering points for courses #2-#Ncourses in each layer (true/false: all steering points or vector of integers specifying the indices of the steering points to be manipulated. N.B: The number of defined steering points for courses #2-#Ncourses must be the same.)
 -DesVar.FirstTransOffset: toggle the transverse/widthwise offset between the first course and the left/right mold net boundary line (i.e. depending on the variable "LayerPropagation") as design variable. Here, the optimizer will manipulate the x-coordinates of the already defined steering points for course #1 in each layer (true/false: all steering points or vector of integers specifying the indices of the steering points to be manipulated. N.B: The number of defined steering points for course #1 must be the same.)
- DesVar.CourseWidth: toggle the width of each course as a design variable. 

The minimization criteria, i.e. objective function settings are stored in the struct Obj:
- Obj.Shear: toggle minimization of shear angles in the objective function (true/false).
- Obj.AngDev: toggle minimization of UD fiber angle deviations in the objective function (true/false).
- Obj.TrimArea: toggle minimization of trim area/waste in the objective function (true/false).
- Obj.p_Shear: value of p in p-norm for shear part of objective function (integer > 1, if inf it corresponds to the maximum value of shear).
- Obj.p_AngDev: value of p in p-norm for UD fiber angle deviation part of objective function (integer > 1, if inf it corresponds to the maximum value of angle deviations).
- Obj.MeanOfLayerPNorms: toggle if the objective function parts should be calculated for each layer and then averaged for the stack. Otherwise global stack values are calculated (true/false).
- Obj.TargetFiberAng: target fiber angle for UD fiber angle deviations (scalar, 0 corresponds to fibers trying to be alligned with the global y-axis).
- Obj.ShearAng_eq: equivalent shear angle used for scaling the trim area term in the objective function (scalar > 0).
- Obj.TrimArea_eq: equivalent trim area used for scaling the trim area term in the objective function (scalar > 0).
Regarding the last two variables, the scaled trim area / waste is calculated as:
W_layer,scaled = Obj.ShearAng_eq * W_layer/Obj.TrimArea_eq

The settings relating to bounds on design variables and constraints are stored in the struct BndAndCon:
- BndAndCon.MinStaggerDist: The minimum stagger distance, i.e. widthwise distance between course edges in adjacent layers (scalar > 0 or empty, i.e. [] for deactivating the constraint).
- BndAndCon.AvailableCourseWidths: The available course widths specified as grid nodes (vector of integers > 0).
- BndAndCon.TransOffset: lower and upper bound on the TransOffset design variable (2-component vector where 1st entry must be smaller than the 2nd).
- BndAndCon.FirstTransOffset: lower and upper bound on the FirstTransOffset design variable (2-component vector where 1st entry must be smaller than the 2nd).

And lastly the genetic algorithm settings in the variable GASet which must call the MATLAB "optimoptions" function to set the values, i.e. (please see the MATLAB Global Optimization Toolbox documentation):
GASet = optimoptions('ga',...
    	'Display','iter',...
    	'MaxGenerations',<maximum number of generations (integer > 0)>,... 
	'MaxStallGenerations',<maximum number of stall generations (integer > 0)>,...
	'PopulationSize', <population size (integer > 0, should at least be > no. of desvar)>,...
	'EliteCount', <number of elite (integer > 0, < PopulationSize)>,...
	'CrossoverFraction', <fraction of the pop. minus elite that reproduces with crossover (scalar > 0, < 1)>,...
    	'UseParallel',<toggle parallel computing (true/false)>,...
    	'FunctionTolerance',<stop if rel. obj. function improvement is less over MaxStallGenerations>,...
    	'PlotFcn',{@gaplotbestf,@gaplotstopping},...
	'NonlinearConstraintAlgorithm','penalty'...
    );

N.B: The MaxGenerations, MaxStallGenerations, PopulationSize and EliteCount are applied per layer in use
