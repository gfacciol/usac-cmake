%
% Created by Gabriele Facciolo on 22/04/16.
%
pts = importdata('../data/fundmatrix/test1/orig_pts.txt', ' ', 1);
pts = pts.data;
sorting = importdata('../data/fundmatrix/test1/sorting.txt');
[F,inliers] = MEX_usac(0, '../data/fundmatrix/example.cfg', true, pts, sorting, 1.5, 10000);
F

pts = importdata('../data/essential/test1/orig_pts.txt', ' ', 1);
pts = pts.data;
sorting = importdata('../data/essential/test1/sorting.txt');
calib = importdata('../data/essential/test1/calib_matrices.txt');
K1 = ParseCalib(calib(1,:));
K2 = ParseCalib(calib(2,:));
X1 = pts(:,1:2)';
X1(3,:) = 1;
X2  = pts(:,3:4)';
X2(3,:) = 1;
nX1 = NormalizedCoordinates(X1,K1);
nX2 = NormalizedCoordinates(X2,K2);
ptsNormalized = [nX1(1:2,:)' nX2(1:2,:)'];
[F,inliers] = MEX_usac(0, '../data/essential/example_unix.cfg', true, ptsNormalized, sorting, 1.5, 10000);

E = K2' * F * K1;
[u, s, v] = svd(E);
if (size(s,1) == size(s,2))
    s = diag(s);
end
avg = (s(1) + s(2))/2;
d = diag([avg avg 0]);
E = u * d * v';


%[E,inliers] = MEX_usac(0, '../data/essential/example.cfg', true, ptsNormalized, sorting, 1.5, 10000);

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
axis equal ij

subplot(1,2,2)
hold on
imagesc(I2); 
plot(pts(inl,3), pts(inl,4), 'go')
epiLines = epipolarLine(F, pts(inl,1:2));
points = lineToBorderPoints(epiLines, size(I2));
line(points(:, [1,3])', points(:, [2,4])');
axis equal ij