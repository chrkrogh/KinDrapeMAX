function [ObjFun, OptRes_S, Dra_S] = StackObj(x,Mold,Plt,Inp,Set,Opt,F,DT,z)
% This function computes the objective function for the stack. It does so
% by calling 'LayerObj' for each layer and processing the output

nLayers = Inp(1).nLayers;

% Initialize
OptRes = cell(1,nLayers);
Dra_S = cell(1,nLayers);
for LNo = 1:nLayers

    if Set.MeshOrCurveWarning
        fprintf(2,'\n|=========== Warning messages for layer #%d ===========|\n',LNo);
    end

    if isfield(Set,'OptDebug') && Set.OptDebug
        % In case of errors in the program during optimization, the
        % complete error message is not always printed and the debugger 
        % will not stop in the function where the error occured. With this
        % try-catch structure, the vector of design variables can be
        % retrieved and the program can be re-run without optimization.
        try
            [~, OptRes{LNo}, Dra_S{LNo}] = LayerObj(x,Mold{LNo},Plt,Inp(LNo),Set,...
                Opt(LNo),F{LNo},DT{LNo},z{LNo});
        catch ME
            x
            save('./Temporary files/ErrorWorkSpace');
            fprintf('\n Workspace saved. Rerun the program in error mode\n')
            rethrow(ME)
        end
    else
        % Just run normally
        [~, OptRes{LNo}, Dra_S{LNo}] = LayerObj(x,Mold{LNo},Plt,Inp(LNo),Set,...
            Opt(LNo),F{LNo},DT{LNo},z{LNo});
    end

end

%% Gather data for optimization

% Combine into one struct
OptRes_S = [OptRes{:}];

% Concatenate vectors from each layer
AllShearAngles_S = [OptRes_S.AllShearAngles];
AllAngDev_S = [OptRes_S.AllAngDev];
TotalTrimArea_S = [OptRes_S.TotalTrimArea];

%% Objective function

if Inp(1).Obj.Shear && ~Inp(1).Obj.MeanOfLayerPNorms
    if Inp(1).Obj.p_Shear == inf
        % Compute the max with all shear angles in the stack
        Shear_obj = max(abs(AllShearAngles_S),[],'omitnan');
    elseif Inp(1).Obj.p_Shear < inf
        % Compute the p-norm with all shear angles in the stack
        p = Inp(1).Obj.p_Shear;
        Shear_obj = sum(abs(AllShearAngles_S).^p,'omitnan')^(1/p);
    end
elseif Inp(1).Obj.Shear && Inp(1).Obj.MeanOfLayerPNorms
    % Take the mean of the p-norms / max computed for each layer
    Shear_obj = mean([OptRes_S.Shear_obj],'omitnan');
else
    Shear_obj = 0;
end

if Inp(1).Obj.AngDev && ~Inp(1).Obj.MeanOfLayerPNorms
    TargetFiberAng = Inp(1).Obj.TargetFiberAng;
    if Inp(1).Obj.p_AngDev == inf
        % Compute the max with all angle deviations in the stack
        AngDev_obj = max(abs(AllAngDev_S-TargetFiberAng),[],'omitnan');
    elseif Inp(1).Obj.p_AngDev < inf
        % Compute the p-norm with all angle deviations in the stack
        p = Inp(1).Obj.p_AngDev;
        AngDev_obj = sum(abs(AllAngDev_S-TargetFiberAng)^p,'omitnan')^(1/p);
    end
elseif Inp(1).Obj.AngDev && Inp(1).Obj.MeanOfLayerPNorms
    % Take the mean of the p-norms / max computed for each layer
    AngDev_obj = mean([OptRes_S.AngDev_obj],'omitnan');
else
    AngDev_obj = 0;
end

if Inp(1).Obj.TrimArea && ~Inp(1).Obj.MeanOfLayerPNorms
    % Scale stack trim area
    ShearAng_eq = Inp(1).Obj.ShearAng_eq;
    TrimArea_eq = Inp(1).Obj.TrimArea_eq;
    TrimArea_obj = ShearAng_eq * sum(TotalTrimArea_S)/TrimArea_eq;
elseif Inp(1).Obj.TrimArea && Inp(1).Obj.MeanOfLayerPNorms
    % Take mean of layer trim areas
    TrimArea_obj = mean([OptRes_S.TrimArea_obj],'omitnan');
else
    TrimArea_obj = 0.0;
end

ObjFun = Shear_obj + AngDev_obj + TrimArea_obj;

% Update the obj parts
[OptRes_S.Shear_obj] = deal(Shear_obj);
[OptRes_S.AngDev_obj] = deal(AngDev_obj);
[OptRes_S.TrimArea_obj] = deal(TrimArea_obj);

end