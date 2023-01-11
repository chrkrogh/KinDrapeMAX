function [Idx, CellNo] = LinearCellIdx(Grid,Row,Col,Dir1,Dir2,No)
% Output: linear indices in Node of cell (4 x 3, i.e. vertices x coord.)
Rows = Row + [0 Dir2(No,1) Dir1(No,1)+Dir2(No,1) Dir1(No,1)]';
Cols = Col + [0 Dir2(No,2) Dir1(No,2)+Dir2(No,2) Dir1(No,2)]';
Idx = Rows + (Cols-1)*Grid(1) + ([1 2 3]-1)*Grid(1)*Grid(2);
CellNo = Rows(No) + (Cols(No)-1)*(Grid(1)-1);
end