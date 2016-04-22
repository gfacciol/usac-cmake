%
% Created by Gabriele Facciolo on 22/04/16.
%
pts = importdata('../data/fundmatrix/test1/orig_pts.txt', ' ', 1);
pts = pts.data;
sorting = importdata('../data/fundmatrix/test1/sorting.txt');
[F,inliers] = MEX_usac(0, '../data/fundmatrix/example.cfg', pts, sorting);
F


% display some matches and epipolar lines as in
% http://fr.mathworks.com/help/vision/ref/epipolarline.html?refresh=true
I1=imread('../data/fundmatrix/test1/im1.jpg');
I2=imread('../data/fundmatrix/test1/im2.jpg');

inl = find(inliers);
inl=inl(1:50:end);

subplot(1,2,1)
hold on
imagesc(I1); 
plot(pts(inl,1), pts(inl,2), 'go')
epiLines = epipolarLine(F', pts(inl,3:4));
points = lineToBorderPoints(epiLines, size(I1));
line(points(:, [1,3])', points(:, [2,4])');
axis equal

subplot(1,2,2)
hold on
imagesc(I2); 
plot(pts(inl,3), pts(inl,4), 'go')
epiLines = epipolarLine(F, pts(inl,1:2));
points = lineToBorderPoints(epiLines, size(I2));
line(points(:, [1,3])', points(:, [2,4])');
axis equal
