function [IntersecPt,t,u] = IntersecOfTwo2DLines(Pt11,Pt12,Pt21,Pt22)
% Find intersection point using vectors
% See "Intersection of 2 2D lines" from Stack Overflow

p = Pt11;
r = Pt12 - Pt11;

q = Pt21;
s = Pt22 - Pt21;

% The expression implemented using Cross2D function
%t = Cross2D((q-p),s) / Cross2D(r,s);
%u = Cross2D((p-q),r) / Cross2D(s,r);

% Inline
qmp = q-p;
pmq = p-q;
t = (qmp(1)*s(2) - qmp(2)*s(1)) / (r(1)*s(2) - r(2)*s(1));
u = (pmq(1)*r(2) - pmq(2)*r(1)) / (s(1)*r(2) - s(2)*r(1));

IntersecPt = p+t*r;

end

% function result = Cross2D(a,b)
% % 2D cross product
% % if length(a) > 2 || length(b) > 2
% %     error('Vectors must have length 2')
% % end
% result = a(1)*b(2) - a(2)*b(1);
% end

