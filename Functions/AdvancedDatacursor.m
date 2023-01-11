function output_txt = AdvancedDatacursor( obj,event_obj) 
dataIndex = get(event_obj,'DataIndex'); 
pos = get(event_obj,'Position'); 
h=get(event_obj,'Target'); 
X=get(h,'XData'); 
Y=get(h,'YData'); 
Z=get(h,'ZData'); 
Clr = get(h,'FaceVertexCData');
idx_x=find(X==pos(1)); 
idx_y=find(Y==pos(2)); 
output_txt = {[ 'X ',num2str(pos(1),4)],... 
[ 'Y ',num2str(pos(2),4)],... 
[ 'Z ',num2str(pos(3),4)],...
strjoin((['C ' string(num2str(Clr(idx_x),4))']))};
%['C ',num2str(mean(Clr(idx_x)),4)]}; 
end