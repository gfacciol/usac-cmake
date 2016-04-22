%
% Created by Gabriele Facciolo on 22/04/16.
%
pts = importdata('../data/fundmatrix/test1/orig_pts.txt', ' ', 1);
pts = pts.data;
sorting = importdata('../data/fundmatrix/test1/sorting.txt');
[F,inliers] = MEX_usac(0, '../data/fundmatrix/example.cfg', pts, sorting);
F
