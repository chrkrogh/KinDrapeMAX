function DispStackResultsInCommandWindow(Inp,OptRes_S,x_opt)
% This function displays results from the design in the command window

if isempty(OptRes_S)
    return
end

nLayers = Inp(1).nLayers;

% Concatenate vectors from each layer
AllShearAngles_S = [OptRes_S.AllShearAngles];
AllAngDev_S = [OptRes_S.AllAngDev];
TotalTrimArea_S = [OptRes_S.TotalTrimArea];

fprintf('\nOptimized result summary\n')

fprintf('Stack shear angles: %.3g deg (min), %.3g deg (max), %.3g deg (avg.)\n\n',...
    min(AllShearAngles_S),max(AllShearAngles_S),mean(AllShearAngles_S))
fprintf('Stack angle dev.: %.3g deg (min), %.3g deg (max), %.3g deg (avg.)\n\n',...
    min(AllAngDev_S),max(AllAngDev_S),mean(AllAngDev_S))
fprintf('Total stack trim area: %.3g m^2 \n\n',sum(TotalTrimArea_S))

if ~isempty(Inp(1).BndAndCon.MinStaggerDist)
    [ViolateLayers,~,~] = find(vertcat(OptRes_S(:).StaggerDist)<Inp(1).BndAndCon.MinStaggerDist);
    ViolateLayers = unique(ViolateLayers);
    if isempty(ViolateLayers)
        fprintf('All stagger dist constraints of %g mm satisfied in stack \n\n',...
            1e3*Inp(1).BndAndCon.MinStaggerDist)
    else
        fprintf('Stagger dist constraint violations in layer(s): ')
        fprintf('%d ',ViolateLayers')
        fprintf('\n\n')
    end
end

%FormatFun = @(x) sprintf('%0.2f', x);

if ~isempty(x_opt)
    T1 = cell(1,nLayers);
    T2 = cell(1,nLayers);
    for LNo = 1:nLayers
        nCourses_curr = OptRes_S(LNo).nCourses_trim;

        fprintf('\n|=========================== Layer #%d ================================|\n',LNo)

        fprintf('\nCriteria\n')
        AllShearAng_L = OptRes_S(LNo).AllShearAngles;
        AllAngDev_L = OptRes_S(LNo).AllAngDev;
        if ~isempty(OptRes_S(LNo).CourseWidths)
            LastCourseArea = ((Inp(LNo).Grid_L{nCourses_curr}(2))-1) * ...
                (Inp(LNo).d_L{nCourses_curr}) * ...
                OptRes_S(LNo).CourseWidths(nCourses_curr);
        else
            LastCourseArea = prod(Inp(LNo).Grid_L{nCourses_curr}-1) * ...
                (Inp(LNo).d_L{nCourses_curr}.^2);
        end
        if ~isempty(Inp(1).BndAndCon.MinStaggerDist)
            StaggerDist_L = OptRes_S(LNo).StaggerDist;
        else
            StaggerDist_L = NaN;
        end

        ShearTabData = [min(AllShearAng_L) ; max(AllShearAng_L) ; ...
            mean(AllShearAng_L)];
        AngDevTabData = [min(AllAngDev_L) ; max(AllAngDev_L) ; ...
            mean(AllAngDev_L)];
        StaggerDistTabData = 1e3*[min(StaggerDist_L,[],'omitnan') ; ...
            max(StaggerDist_L,[],'omitnan') ; mean(StaggerDist_L,'omitnan')];

        % Create table
        T1{LNo} = array2table([ShearTabData, AngDevTabData, ...
            StaggerDistTabData]');
        % Change each entry in the tale to a string to allow number formatting
        T1{LNo} = varfun(@(x) num2str(x, ['%' sprintf('.%df', 2)]), T1{LNo});
        % Set variable names (N.B table array is transposed)
        T1{LNo}.Properties.VariableNames = {'Min','Max','Avg'};
        % Set row names
        T1{LNo}.Properties.RowNames = {'Shear angles [deg]',...
            'Angle dev. [deg]','Stagger dist [mm]'};
        if LNo == 1
            % Remove stagger dist information (not applicable to 1st layer)
            T1{1}(3,:) = [];
        end

        disp(T1{LNo})

        fprintf('\tTotal trim area: %.3g m^2 / %.3g %% of course #%d\n\n',...
            OptRes_S(LNo).TotalTrimArea,100*OptRes_S(LNo).TotalTrimArea/LastCourseArea,...
            nCourses_curr)


        fprintf('\nDesign\n')
        % Extract the information from OptRes_S struct
        [~,OrgNodePct_opt] = OrgNodeFromPct(OptRes_S(LNo).OrgNodeLong,[],...
            cellfun(@(x) x(2),Inp(LNo).Grid_L(1:nCourses_curr)));
        Offsets_opt = OptRes_S(LNo).TransOffsets;
        FirstOffsets_opt = OptRes_S(LNo).FirstTransOffsets;
        nTransOffsets = max(size(Offsets_opt,2),size(FirstOffsets_opt,2));
        if isempty(OptRes_S(LNo).FirstTransOffsets)
            % Add a NaN row for first course (i.e. no des vars) and conv. to mm
            Offsets_opt_mm = [NaN(1,size(Offsets_opt,2)) ; 1e3*Offsets_opt];
        elseif isempty(OptRes_S(LNo).TransOffsets)
            Offsets_opt_mm = [1e3*FirstOffsets_opt ; ...
                NaN(OptRes_S(LNo).nCourses_trim-1,size(FirstOffsets_opt,2))];
        elseif ~isempty(OptRes_S(LNo).FirstTransOffsets) && ...
                ~isempty(OptRes_S(LNo).TransOffsets)
            Offsets_opt_mm = NaN(OptRes_S(LNo).nCourses_trim,nTransOffsets);
            
            Offsets_opt_mm(1,1:size(FirstOffsets_opt,2)) = 1e3*FirstOffsets_opt;
            Offsets_opt_mm(2:end,1:size(Offsets_opt,2)) = 1e3*Offsets_opt;
        end

        Widths_opt_mm = 1e3*OptRes_S(LNo).CourseWidths';

        T2{LNo} = array2table([OrgNodePct_opt' Offsets_opt_mm Widths_opt_mm]);
        % Merge variables from transverse offset into one column
        T2{LNo} = mergevars(T2{LNo},2:2+nTransOffsets-1,...
            'NewVariableName','Trans offset');
        % Change each entry in table to string to allow number formatting
        T2{LNo} = varfun(@(x) num2str(x, ['%' sprintf('.%df ',2)]),T2{LNo});
        % Set variable and row names

        Varname_cell = {};
        if ~isempty(OrgNodePct_opt)
            Varname_cell = [Varname_cell 'OrgNode [%]'];
        end

        if ~isempty(Offsets_opt_mm)
            Varname_cell = [Varname_cell 'Trans. offset [mm]'];
        end

        if ~isempty(Widths_opt_mm)
            Varname_cell = [Varname_cell 'Width [mm]'];
        end

        T2{LNo}.Properties.VariableNames = Varname_cell;

        T2{LNo}.Properties.RowNames = arrayfun(@(x) ['Course #' ...
            num2str(x)],1:length(OrgNodePct_opt),'Uni',0);
        % Remove trans offset information from 1st course (not applicable)
        if ~isempty(OptRes_S(LNo).TransOffsets) && isempty(OptRes_S(LNo).FirstTransOffsets)
            T2{LNo}{1,2} = ' ';
        elseif isempty(OptRes_S(LNo).TransOffsets) && ~isempty(OptRes_S(LNo).FirstTransOffsets)
            T2{LNo}{2:OptRes_S(LNo).nCourses_trim,2} = ' ';
        end
        disp(T2{LNo})
    end
end
end