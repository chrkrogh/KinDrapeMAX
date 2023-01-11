%   _   ___      ______                     ___  ___  ___  __   __
%  | | / (_)     |  _  \                    |  \/  | / _ \ \ \ / /
%  | |/ / _ _ __ | | | |_ __ __ _ _ __   ___| .  . |/ /_\ \ \ V / 
%  |    \| | '_ \| | | | '__/ _` | '_ \ / _ \ |\/| ||  _  | /   \ 
%  | |\  \ | | | | |/ /| | | (_| | |_) |  __/ |  | || | | |/ /^\ \
%  \_| \_/_|_| |_|___/ |_|  \__,_| .__/ \___\_|  |_/\_| |_/\/   \/
%                                | |                              
%                                |_|               
%
%
% Written by Christian Krogh, Ph.D 
% Department of Materials and Production, Aalborg University, Denmark.
% Licensed under MIT license. Copyright (c) 2023 Christian Krogh
%
% The program can analyze and optimize the draping pattern on a mold using 
% a kinematic draping model and the built-in MATLAB implementation of 
% genetic algorithm, i.e. the function ga. It is intended for wind turbine 
% blade production, i.e. where multiple roll-widths or courses of fabric 
% are necessary to cover the mold and where the mold is more or less 
% rectangular. It can, however, also analyze a single ply on any smooth 
% surface. Please see the user guide for more information on how to run it.

clc
clear
close all

rng(0); % Seed for random number generator, i.e. 0,1,2,3...

Inp.MoldName = 'BladeLike___FromMesh_5layer';

Set.Mode = 'analysis'; % opt, opt-seq, baseline, analysis, storedopt, error
Set.InpFileType = 'course'; % course, layer, stack, -opt

Set.TrimPlyToNetBd = true;
Set.ExtrapSteerCurves = true;
Set.MeshOrCurveWarning = false;
Set.ForceGeoLineIniAngle = true;
Set.OptDebug = false;
Set.Timing = 'time'; % 'time' or 'profile'
Set.CheckAspect = false;
Set.AspectLim = 20;

Plt.PlyPlt = {'shear','fibdev'}; % 'shear','flat','fibDev'
Plt.NodesWithoutTrim = false;
Plt.AbsShearAng = false;
Plt.AbaqusColor = false;
Plt.ShearLim = 5;
Plt.TriMesh = false;
Plt.MoldNetBd = true;
Plt.GenCurve = true;
Plt.GeoInc = false;
Plt.MeshBoundary = true;
Plt.DispLegend = true;
Plt.zShift = 1.0e-3;
Plt.FigWindowSize = 'normal'; %'normal', 'full', 'full-other', 'half'

Plt.Display = true;

addpath('./Functions/')
addpath('./Aux tools/')

% Initialize and get mold and input files
[Inp,Set,Mold,F,DT,z,FileList,InpFilePath,nCourses,nLayers] = ...
    InitializeAndGetMoldAndInpFile(Inp,Set);

% Process optimization input
Opt = ProcessOptInput(Inp,Set);

% Analyze or optimize draping model
[OptProb,ObjFun,OptRes_S,Dra_S,x_opt,nCourses_trim] = ...
    AnalyzeOrOptimizeDrapingModel(Inp,Set,Opt,Mold,Plt,F,DT,z,nLayers);

% Plot results in figures
Plt = PlotResultsInFigs(Inp,Set,Plt,Dra_S,Mold,F,DT,z,nLayers,...
    nCourses_trim,OptRes_S);

% Display results in command window
DispStackResultsInCommandWindow(Inp,OptRes_S,x_opt)

% Store results
StoreResults(Inp,Plt,InpFilePath,x_opt,Dra_S,DT,F,z,OptRes_S,Mold,...
    ObjFun,Opt,OptProb,Set)
