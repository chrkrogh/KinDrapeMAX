function [kg,kn,k] = CheckGeodesicCurvature(Curve,DT,z,PlotRes,Step)
% Check the geodesic curvature of a curve defined by points on a mold. It
% is assumed that the points are all located on the mold.
% The idea is to compute a unit normal to the surface, a unit tangent
% vector to the curve (i.e. the velocity vector since the curve is unit 
% speed) and from those to a mutually perpendicular vector.
% The acceleration vector is then computed. The acceleration vector 
% projected onto the normal vector is the normal curvature and the 
% acceleration vector projected onto the mutually perpendicular vector is
% the geodesic curvature (see Section 7.3 (p. 166) in Pressly (2010)).

if nargin <= 4
    Step = 1;
end

PlotIntermediateResult = false;

if PlotIntermediateResult
    figure
    hold on
    trimesh(DT.ConnectivityList,DT.Points(:,1),...
        DT.Points(:,2),z(:),'LineWidth',0.1,'FaceColor','none','EdgeColor','k')
    axis equal
    view(3)
end
    
kg = cell(1,size(Curve,2));
kn = cell(1,size(Curve,2));
k = cell(1,size(Curve,2));

if PlotRes
   figure 
end

for CurveNo = 1:size(Curve,2)

    if Step > 1 && ~isempty(Curve{CurveNo})
        Curve{CurveNo}(1:Step:end,:) = [];
    end
    
    if PlotIntermediateResult
        plot3(Curve{CurveNo}(:,1),Curve{CurveNo}(:,2),Curve{CurveNo}(:,3),'r-o')
    end
        
    Ctr = 0;
    for i = 2:1:size(Curve{CurveNo},1)-1
       
        Ctr = Ctr + 1;
        
        % The current and two next points
        Pt_curr = Curve{CurveNo}(i,:);
        Pt_p1 = Curve{CurveNo}(i+1,:); 
        Pt_m1 = Curve{CurveNo}(i-1,:); 
        
        % Corresponding triangle of the mold mesh
        % Find the triangle ID
        TriID = pointLocation(DT,Pt_curr(1:2));
        if isnan(TriID)
            break
        end
        % Find the vertex indices
        VertIdx = DT.ConnectivityList(TriID,:);
        % Calculate the current vertex coordinates
        VertCoor(1:3,1:2) = DT.Points(VertIdx,:);
        VertCoor(1:3,3) = F_DT(DT,z,VertCoor(1:3,1),VertCoor(1:3,2));
        
        % Calculate the normal vector to the surface at the point
        % First get two edge vectors
        EdgeVec(1,1:3) = VertCoor(2,1:3) - VertCoor(1,1:3);
        EdgeVec(2,1:3) = VertCoor(3,1:3) - VertCoor(2,1:3);
        % Compute the cross product
        NormVec = cross(EdgeVec(1,1:3),EdgeVec(2,1:3));
        NormVec = NormVec / norm(NormVec);
        
        % Calculate a unit tangent (velocity) vector to the curve
        % using forward finite difference
        TangVec = (Pt_p1 - Pt_curr) / norm (Pt_p1 - Pt_curr);

        % The mutually perpendicular vector
        PerpVec = cross(NormVec,TangVec);
        
        % Calculate the acceleration as the 2nd order central finite
        % difference approximation. It is based on the FD of two 1st order
        % central FD's
        delta1 = norm(Pt_p1 - Pt_curr);
        delta2 = norm(Pt_curr - Pt_m1);
        
        FD_1st_1 = (Pt_p1 - Pt_curr)/delta1;
        FD_1st_2 = (Pt_curr - Pt_m1)/delta2;
        
        % Calculate the 2nd order FD with a step size as the average of the
        % two
        delta3 = 0.5*(delta1 + delta2);
        Acc = (FD_1st_1 - FD_1st_2) / delta3;
        
        % Calculate the curvatures
        k{CurveNo}(Ctr) = norm(Acc);
        kg{CurveNo}(Ctr) = abs(dot(Acc,PerpVec));
        kn{CurveNo}(Ctr) = abs(dot(Acc,NormVec));
        
        if PlotIntermediateResult
            Sc = 0.25;
            
           quiver3(Pt_curr(1),Pt_curr(2),Pt_curr(3),...
               NormVec(1),NormVec(2),NormVec(3),Sc,'b-') 
           
           quiver3(Pt_curr(1),Pt_curr(2),Pt_curr(3),...
               TangVec(1),TangVec(2),TangVec(3),Sc,'g-')
           
           quiver3(Pt_curr(1),Pt_curr(2),Pt_curr(3),...
               PerpVec(1),PerpVec(2),PerpVec(3),Sc,'y-')
           
           quiver3(Pt_curr(1),Pt_curr(2),Pt_curr(3),...
               Acc(1),Acc(2),Acc(3),Sc,'m-')
           
           keyboard
        end
    end 
    if PlotRes
        
        subplot(2,2,CurveNo)
        hold on
        plot(kg{CurveNo},'b-o')
        
        plot(kn{CurveNo},'r-s')
        
        plot(k{CurveNo},'g-x')
        
        xlabel('Vector index')
        ylabel('Curvature')
        
        legend('Geo. curvature, computed',...
            'Normal curvature, computed',...
            'Curvature, computed',...
            'Location','Best')
    end
end
end
