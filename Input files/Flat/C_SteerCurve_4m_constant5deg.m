% -------------------------------------------------------------------
%  Generated by MATLAB on 26-May-2022 09:25:22
%  MATLAB version: 9.12.0.1927505 (R2022a) Update 1
% -------------------------------------------------------------------
                                                       

d = 0.1;

Grid = [6 41];

Org = [1.5 0];

Ang = 0;

OrgNode = [1 1];

PreShear = 5;

DraContr = cell(1, 2);
DraContr{1} = 'Steer';
DraContr{2} = 'Geo';

SteerPtsRef = cell(2, 1);
SteerPtsRef{1} = 'Abs';
SteerPtsRef{2} = 'Bottom';

SteerPts = cell(2, 1);
SteerPts{1,1} = ...
  [1.5 1.5 1.5;
   0 2 4];
SteerPts{2,1} = [0; 0];

