function [K] = ParseCalib(calib)
    K = eye(3);
    K(1,1:3) = calib(1:3);
    K(2,2:3) = calib(4:5);
end