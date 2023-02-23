function FigHan = Plot2DFlattenedStack(Inp,Mold,Set,OptRes_S)
% This function takes the output from the draping analysis or optimization
% and plots the flattened stack at the root end. That is, the flattened 
% mold surface is a straight line with length equal to the widthwise mold 
% arc length at the root. The courses are plotted as straight lines. Course
% ends and parts of courses that have been trimmed are also plotted. This
% function could be extended with flattened stack plots at other sections /
% y-coordinates.
% It is assumed that the 1st steering point is at the root end. The
% implementation is rather quick and dirty. There could be other bugs if
% not used as in the paper.

% First compute length of bottom mold edge (root) for the ref plane
MoldRootLength_ref = sum(sqrt(...
    diff(Mold{1}.MoldEdge{3}(:,1)).^2 + ...
    diff(Mold{1}.MoldEdge{3}(:,2)).^2 + ...
    diff(Mold{1}.MoldEdge{3}(:,3)).^2));

FigHan = figure;
hold on

PlotOffset = 50e-3;

% Plot the flattened mold surface
if strcmpi(Inp(1).LayerPropagation,'right-to-left')
    plot([0 -MoldRootLength_ref],[0 0],'k-','LineWidth',2)
elseif strcmpi(Inp(1).LayerPropagation,'left-to-right')
    plot([0 MoldRootLength_ref],[0 0],'k-','LineWidth',2)
end
for LNo = 1:Inp(1).nLayers

    MoldRootLength_curr = sum(sqrt(...
        diff(Mold{LNo}.MoldEdge{3}(:,1)).^2 + ...
        diff(Mold{LNo}.MoldEdge{3}(:,2)).^2 + ...
        diff(Mold{LNo}.MoldEdge{3}(:,3)).^2));

    NewXCoor = 0.0;
    for CNo = 1:OptRes_S(LNo).nCourses_trim
        % Course widths
        if ~isempty(OptRes_S(LNo).CourseWidths)
            % From optimization
            CourseWidth_curr = OptRes_S(LNo).CourseWidths(CNo);
        else
            % From input file
            CourseWidth_curr = (Inp(LNo).Grid_L{CNo}(1)-1) * ...
                Inp(LNo).d_L{CNo};
        end

        if CNo == 1 && strcmpi(Inp(LNo).LayerPropagation,'right-to-left') && ...
                Inp(LNo).OrgNodePct_L{CNo}(1) == 0
            % If the 1st course has the steering curve along the left edge and the
            % layer propagates from right to left, introduce a course width
            % offset
            CourseWidthOffset = CourseWidth_curr;
        elseif CNo == 1 && strcmpi(Inp(LNo).LayerPropagation,'left-to-right') && ...
                Inp(LNo).OrgNodePct_L{CNo}(1) == 100
            % If the 1st course has the steering curve along the right edge and the
            % layer propagates from left tor right, introduce a course width
            % offset
            CourseWidthOffset = -CourseWidth_curr;
        else
            CourseWidthOffset = 0.0;
        end

        % Trans offsets incl. FirstTransOffset
        if ~isempty(OptRes_S(LNo).FirstTransOffsets) && CNo == 1
            MinGap_curr = min(OptRes_S(LNo).FirstTransOffsets(1,:));
            MaxOverlap_curr = max(OptRes_S(LNo).FirstTransOffsets(1,:));
            RootPt = OptRes_S(LNo).FirstTransOffsets(1,1);
        elseif isempty(OptRes_S(LNo).TransOffsets) || CNo == 1
            MinGap_curr = min(Inp(LNo).SteerPts_L{1,CNo}(1,:)) + CourseWidthOffset;
            MaxOverlap_curr = max(Inp(LNo).SteerPts_L{1,CNo}(1,:)) + CourseWidthOffset;
            RootPt = Inp(LNo).SteerPts_L{1,CNo}(1,1) + CourseWidthOffset;
        else
            MinGap_curr = min(OptRes_S(LNo).TransOffsets(CNo-1,:));
            MaxOverlap_curr = max(OptRes_S(LNo).TransOffsets(CNo-1,:));
            RootPt = OptRes_S(LNo).TransOffsets(CNo-1,1);
        end

        % Calculate the min gap and max overlap
        if ~isempty(OptRes_S(LNo).FirstTransOffsets) && CNo == 1 ...
                && strcmpi(Inp(LNo).LayerPropagation,'right-to-left')
            MinGap_curr = MaxOverlap_curr;
            MaxOverlap_curr = 0;
        elseif strcmpi(Inp(LNo).LayerPropagation,'right-to-left')
            MinGap_curr = min(0,MinGap_curr);
            MaxOverlap_curr = max(0,MaxOverlap_curr);
        elseif strcmpi(Inp(LNo).LayerPropagation,'left-to-right')
            MinGap_curr = max(0,MaxOverlap_curr);
            MaxOverlap_curr = min(0,MinGap_curr);
        end

        % If the 1st course uses the left edge as steering curve
        if CNo == 1 && Inp(LNo).OrgNodePct_L{1}(1) == 0 ...
                && strcmpi(Inp(LNo).LayerPropagation,'right-to-left')
            NewXCoor = CourseWidth_curr;
            MinGap_curr = MinGap_curr - CourseWidth_curr;
            RootPt = RootPt - CourseWidth_curr;
        end

        PrevXCoor = NewXCoor + RootPt;
        if strcmpi(Inp(LNo).LayerPropagation,'right-to-left')
            NewXCoor = PrevXCoor - CourseWidth_curr;
        elseif strcmpi(Inp(LNo).LayerPropagation,'left-to-right')
            NewXCoor = PrevXCoor + CourseWidth_curr;
        end

        % Check trimming
        if strcmpi(Inp(LNo).LayerPropagation,'right-to-left') && ...
                Set.TrimPlyToNetBd
            % Check if trimming was necessary in right side (Q&D)
            if CNo == 1 && PrevXCoor > 0
                WasteCoor_r = PrevXCoor;
                PrevXCoor = 0;
                plot([NewXCoor WasteCoor_r],PlotOffset*[LNo LNo],'r-')
            end
            % Check if trimming was necessary in left side (Q&D)
            if NewXCoor < -MoldRootLength_curr
                WasteCoor_l = NewXCoor;
                NewXCoor = -MoldRootLength_curr;
                plot([NewXCoor WasteCoor_l],PlotOffset*[LNo LNo],'r-')
                % Trick to avoid issues after the net boundary has been
                % passed
                if PrevXCoor < NewXCoor
                    PrevXCoor = NewXCoor;
                end
            end
        elseif strcmpi(Inp(LNo).LayerPropagation,'left-to-right') && ...
                Set.TrimPlyToNetBd
            % Check if trimming was necessary in left side (Q&D)
            if CNo == 1 && PrevXCoor < 0
                WasteCoor_r = PrevXCoor;
                PrevXCoor = 0;
                plot([NewXCoor WasteCoor_r],PlotOffset*[LNo LNo],'r-')
            end
            % Check if trimming was necessary in left side (Q&D)
            if NewXCoor > MoldRootLength_curr
                WasteCoor_r = NewXCoor;
                NewXCoor = MoldRootLength_curr;
                plot([NewXCoor WasteCoor_r],PlotOffset*[LNo LNo],'r-')
                % Trick to avoid issues after the net boundary has been
                % passed
                if PrevXCoor > NewXCoor
                    PrevXCoor = NewXCoor;
                end
            end
        end

        % Plot the course at the root
        plot([PrevXCoor NewXCoor],PlotOffset*[LNo LNo],'b-')
        plot([PrevXCoor NewXCoor],PlotOffset*[LNo LNo],'b|',...
            'MarkerFaceColor','b','MarkerSize',12)
        % Plot the maximum (abs) overlap value
        if CNo > 1 && MaxOverlap_curr > 0
            plot([PrevXCoor- MinGap_curr PrevXCoor- ...
                MinGap_curr+MaxOverlap_curr],PlotOffset*[LNo LNo],'g-')
            plot(PrevXCoor- MinGap_curr+MaxOverlap_curr*[1 1],...
                PlotOffset*[LNo-0.3 LNo+0.3],'g-')
        end
    end
end

axis tight
ylim([-PlotOffset PlotOffset*(Inp(1).nLayers+1)])
pbaspect([3 1 1])
set(gca,'clipping','off')
set(gcf,'color','w');
set(gca,'position',[.05 .05 .9 .9])
axis off
title('Flattened stack plot at root')
pause(1e-6)

end