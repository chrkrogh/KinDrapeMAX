function [OptProb,ObjFun,OptRes_S,Dra_S,x_opt,nCourses_trim] = ...
    AnalyzeOrOptimizeDrapingModel(Inp,Set,Opt,Mold,Plt,F,DT,z,nLayers)
% Optimize, view baseline, view stored optimization result, analyze (w/o
% optimization framework) or handle error from debug mode

% Unpack the substruct with design variable idx
if strcmpi(Set.Mode,'Analysis')
    DesVarIdxSubStruct = [];
else
    DesVarIdxSubStruct = [Opt.DesVarIdx];
end

switch lower(Set.Mode) % Compare only lower case values

    case 'opt'
        [Opt.SeqLayer] = deal(false);

        % Create optimization problem structure
        OptProb.fitnessfcn = @(x)StackObj(x,Mold,Plt,Inp,Set,Opt,F,DT,z);
        OptProb.Aineq = [];
        OptProb.Bineq = [];
        OptProb.Aeq	= [];
        OptProb.Aeq = [];
        OptProb.Beq	= [];
        OptProb.lb = cell2mat(reshape(squeeze(struct2cell([Opt.lb])),1,[]));
        OptProb.ub = cell2mat(reshape(squeeze(struct2cell([Opt.ub])),1,[]));
        OptProb.nvars = length(OptProb.lb);
        
        OptProb.intcon = [DesVarIdxSubStruct.OrgNodeLong ...
            DesVarIdxSubStruct.CourseWidth];
        OptProb.rngstate = [];
        OptProb.solver = 'ga';
        OptProb.options = Inp.GASet;

        % Set the total MaxGenerations, MaxStallGenerations, PopulationSize 
        % and EliteCount, i.e. multiply by nLayers
        OptProb.options.MaxGenerations = ...
            OptProb.options.MaxGenerations * nLayers;
        OptProb.options.MaxStallGenerations = ...
            OptProb.options.MaxStallGenerations * nLayers;
        OptProb.options.PopulationSize = ...
            OptProb.options.PopulationSize * nLayers;
        OptProb.options.EliteCount = ...
            OptProb.options.EliteCount * nLayers;

        if isempty(Inp(1).BndAndCon.MinStaggerDist)
            OptProb.nonlcon	= [];
        else
            OptProb.nonlcon	= @(x)StackNonlCon(x,Mold,Plt,Inp,Set,Opt,F,DT,z);
        end

        [x_opt,fval,exitflag,output] = ga(OptProb);

    case 'opt-seq'

        nVarsTot = sum(cellfun(@numel,squeeze(struct2cell([Opt.lb]))),'all');

        x_opt(1:nVarsTot) = NaN;

        % Put "true" in all rows of the struct in the field SeqLayer
        [Opt.SeqLayer] = deal(true);

        Opt(1).PrevIdx = 0;

        for LNo = 1:nLayers

            fprintf(['\n|=========================== Layer #%d ========'...
                '========================|\n'],LNo)

            % Create optimization problem structure
            OptProb.fitnessfcn = @(x)LayerObj(x,Mold{LNo},Plt,Inp(LNo),...
                Set,Opt(LNo),F{LNo},DT{LNo},z{LNo});
            OptProb.Aineq = [];
            OptProb.Bineq = [];
            OptProb.Aeq	= [];
            OptProb.Aeq = [];
            OptProb.Beq	= [];
            OptProb.lb = cell2mat(struct2cell(Opt(LNo).lb)');
            OptProb.ub = cell2mat(struct2cell(Opt(LNo).ub)');
            OptProb.nvars = length(OptProb.lb);
            
            OptProb.intcon = [Opt(LNo).DesVarIdx.CourseWidth ...
                              Opt(LNo).DesVarIdx.OrgNodeLong] - ...
                              Opt(LNo).PrevIdx;
            OptProb.rngstate = [];
            OptProb.solver = 'ga';
            OptProb.options = Inp.GASet;
            if LNo == 1 || isempty(Inp(1).BndAndCon.MinStaggerDist)
                OptProb.nonlcon	= [];
            else
                OptProb.nonlcon	= @(x)StackNonlCon(x,Mold{LNo},Plt,...
                    Inp(LNo),Set,Opt(LNo),F{LNo},DT{LNo},z{LNo});
            end

            % Optimize!
            [x_opt_iter,fval,exitflag,output] = ga(OptProb);

            % Store in global vector design variables
            x_opt(Opt(LNo).PrevIdx+(1:OptProb.nvars)) = x_opt_iter;

            if LNo < nLayers
                % Update the PrevIdx
                Opt(LNo+1).PrevIdx = Opt(LNo).PrevIdx+OptProb.nvars;

                % Evaluate the model to get the course edges for the
                % stagger distance constraint in next iteration
                [~, OptRes_prev, Dra_L_prev] = LayerObj(x_opt_iter,Mold{LNo},...
                    Plt,Inp(LNo),Set,Opt(LNo),F{LNo},DT{LNo},z{LNo});

                % Store values in Opt struct for next opt. iteration
                Opt(LNo+1).Dra_L_prev = Dra_L_prev;
                Opt(LNo+1).nCourses_prev = OptRes_prev.nCourses_trim;
                Opt(LNo+1).MaxLongNodes_prev = max(cellfun(@(x) x(2),...
                    [Inp(LNo).Grid_L]));
            else
                % Set the SeqLayer field to false;
                % (otherwise the LayerObj function will still offset the
                % vector of design variables in the following evaluations)
                [Opt.SeqLayer] = deal(false);
            end

        end

    case {'baseline','analysis'}
        x_opt = [];
        OptProb = [];
        for k = 1:nLayers
            Opt(k).DesVarIdx.OrgNodeLong = [];
            Opt(k).DesVarIdx.TransOffset = [];
            Opt(k).DesVarIdx.FirstTransOffset = [];
            Opt(k).DesVarIdx.CourseWidth = [];
        end
    case 'storedopt'
        x_opt = Inp.x_opt;
        OptProb = [];
    case 'error'
        load('./Temporary files/ErrorWorkSpace.mat','x'); 
        x_opt = x; 
        clear x;  
        OptProb = [];
        Set.OptDebug = false;
        Set.MeshOrCurveWarning = true;
        dbstop if error
end

if strcmpi(Set.Mode,'Analysis')
    % Compute the layers directly
    % Initialize
    Dra_S = cell(1,nLayers);
    nCourses_trim = zeros(1,nLayers);
    % Loop over all layers
    for LNo = 1:nLayers
        [Dra_S{LNo}, nCourses_trim(LNo)] = ...
            ComputeLayer(Mold{LNo},Plt,Inp(LNo),Set,F{LNo},DT{LNo},z{LNo});
    end
    OptRes_S = [];
    ObjFun = [];
else
    % Run through the objective and possibly constraint function
    [ObjFun, OptRes_S, Dra_S] = StackObj(x_opt,Mold,Plt,Inp,Set,Opt,F,DT,z);

    nCourses_trim = [OptRes_S.nCourses_trim];

    if ~isempty(Inp(1).BndAndCon.MinStaggerDist)
        % Disable warnings and inc. plotting of geodesic line creation
        % (otherwise they will be printed/plotted again with StackNonlCon)
        Set.MeshOrCurveWarning = false;
        Plt.GeoInc = false;
        
        [C,~,StgDist_cell] = StackNonlCon(x_opt,Mold,Plt,Inp,Set,Opt,F,DT,z);

        % Store the stagger distance in the struct
        [OptRes_S.StaggerDist] = StgDist_cell{:};
    end
end

if strcmpi(Set.Timing,'time')
    toc
elseif strcmpi(Set.Timing,'profile')
    profile viewer
    profile off
end
end