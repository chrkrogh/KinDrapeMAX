function Plt = PlyPlot(Dra,Plt,Org,F)
% Plot the ply using various color options
% Update the min/max shear and fib. dev. values for the color limits
Plt.MinShear = min(Plt.MinShear,min(min(Dra.P(:,:,4))));
Plt.MaxShear = max(Plt.MaxShear,max(max(Dra.P(:,:,4))));

Plt.MinFibDev = min(Plt.MinFibDev,min(min(Dra.P(:,:,5))));
Plt.MaxFibDev = max(Plt.MaxFibDev,max(max(Dra.P(:,:,5))));

for i = 1:length(Plt.PlyPlt)
    if isfield(Dra,'P')
        if strcmpi(Plt.PlyPlt{i},'Shear')
            figure(Plt.ShearFig);
            if Plt.AbsShearAng || Plt.AbaqusColor
                Dra.P(:,:,4) = abs(Dra.P(:,:,4));
            end
            patch(Dra.P(:,:,1)',Dra.P(:,:,2)',Dra.P(:,:,3)'+Plt.zShift,...
                Dra.P(:,:,4)','EdgeColor','k');
            
            % Set the color limits
            if Plt.AbaqusColor
                NewCmap = [ones(50,1)*[0 0 1] ; ones(49,1)*[1 1 0] ; [1 0 0]];
                colormap(NewCmap)
                caxis([0 Plt.ShearLim])
                view([-6.3458 10.5085])
            elseif ~Plt.AbaqusColor
                colormap('jet');
                if ~all(all(isnan(Dra.P(:,:,4)))) && ~Plt.AbsShearAng 
                    caxis([Plt.MinShear Plt.MaxShear])
                elseif ~all(all(isnan(Dra.P(:,:,4)))) && Plt.AbsShearAng
                    caxis([0 max(abs(Plt.MinShear),Plt.MaxShear)])
                end
            end
            
        elseif strcmpi(Plt.PlyPlt{i},'Flat')
            figure(Plt.FlatFig);
            
            FlatPlyColor = {'k','r','g','b','y','c','m'};
            
            patch(Dra.P(:,:,1)',Dra.P(:,:,2)',Dra.P(:,:,3)'+Plt.zShift,0,...
                'EdgeColor',FlatPlyColor{Plt.PlyColorNo},'FaceColor','none');
            
        elseif strcmpi(Plt.PlyPlt{i},'FibDev')
            figure(Plt.FibDevFig);
            
            patch(Dra.P(:,:,1)',Dra.P(:,:,2)',Dra.P(:,:,3)'+Plt.zShift,...
                Dra.P(:,:,5)');
            
            colormap('jet');
            if ~all(all(isnan(Dra.P(:,:,5))))
                caxis([Plt.MinFibDev Plt.MaxFibDev])
            end
            
        end
    end
    
    if ~isempty(Org)
        % Plot the (adjusted) origin
        scatter3(Org(1),Org(2),F(Org(1),Org(2)),'kx','LineWidth',5);
    end
    
    if isfield(Plt,'NodesWithoutTrim') && Plt.NodesWithoutTrim
        plot3(Dra.Node(:,:,1),Dra.Node(:,:,2),Dra.Node(:,:,3),'mo')
    end

    % Generator curves
    if isfield(Dra,'GenCurve') && Plt.GenCurve
        for GNo = find(~cellfun(@isempty,Dra.GenCurve))
            plot3(Dra.GenCurve{1,GNo}(:,1),...
                Dra.GenCurve{1,GNo}(:,2),...
                Dra.GenCurve{1,GNo}(:,3),'r--','LineWidth',3);
            if (GNo == 1 || GNo == 3) && ~isempty(Dra.SteerPtsAbs{2})
                SteerPts_curr = Dra.SteerPtsAbs{2};
            elseif (GNo == 2 || GNo == 4) && ~isempty(Dra.SteerPtsAbs{1})
                SteerPts_curr = Dra.SteerPtsAbs{1};
            else
                continue
            end
            scatter3(SteerPts_curr(1,:),SteerPts_curr(2,:),...
                F(SteerPts_curr(1,:),SteerPts_curr(2,:)),100,'ro');
        end
    end

end
end