%&%% DIGITAL IMAGING PROCESSING Exercise 3
% Karam Mawas	    2946939 	

clear all
close all
clc

O_coord = [512868.940 512982.899 513047.601 512781.152 512866.691 512990.398 512995.967
     5427723.833 5427801.281 5427704.443 5427721.195 5427556.801 5427681.815 5427683.424
     280.885 323.376 326.998 229.286 229.433 330.030 326.758
     1 1 1 1 1 1 1];
im_coord = [3182.89734794033 735.934408729377 626.040646560719 3883.00857598429 4277.44761332331 1950.62337050776 1831.15016856533
      2309.95058488979 2663.28816090427 236.234572089787 2654.66673210226 235.689027402227 533.919832052842 502.625318577091
      1 1 1 1 1 1 1];

NPts = length(im_coord); % Number of points

A = zeros(2*NPts,12);
o = [0 0 0 0];
for n = 1:NPts
    X = O_coord(:,n)';
    x = im_coord(1,n);
    y = im_coord(2,n);
    w = im_coord(3,n);
    % building the A matrix
    A(2*n-1,:) = [ o -w*X y*X];
    A(2*n,:) = [ w*X o x*X];
end

% SVD Singular Value Decomposition
[U,D,V] = svd(A,0);
% Extract Homography
P = reshape(V(:,12),4,3)';

% finding the null vector X_0
X0 = P(:,4); % the t matrix wich contains the translations (t=-R.X~0)

cal_cam = P(:,1:3); % the camera calibration  K & R matrices


%MATLAB QR factorization
%returns an upper triangular matrix R and
%a unitary matrix Q, such that A = Q*R
cal_cam_inv = inv(cal_cam); % inverting the P matrix in order to get the right values
[R,K] = qr(cal_cam_inv); % Q=R, R=K % 

% the output would be inverted, so we have to inverse the matrix 
R = R'; % the inverse of the rotation matrix is the same as the transpose
K = inv(K); K = K./K(3,3); % K matrix should be normalized in order to be homogeneous
K(1,:) = -K(1,:);

% Ploting the points depending on the computed P matrix
im_point = P*O_coord;
im_point_x = -im_point(1,:)./im_point(3,:);
im_point_y = im_point(2,:)./im_point(3,:);
im_point_z = im_point(3,:)./im_point(3,:);

im_point = [im_point_x; im_point_y; im_point_z]

%%
data = dlmread('Signalized_Points_R0020851_test.txt','\t');
label = cellstr(num2str(data(:,1)));

c = 21.63348509; % focal length [mm]
delta_pix = 0.0055; % the pixel size [mm]
mx = 1/delta_pix;
px = 2143.5; % principal point in x-direction [pixel]
py = -1423.5; % principal point in y-direction [pixel]


num = data(:,1); % the points number
X = data(:,2:4); % the points coordinates
X = [X ones(7,1)]; % the homogeneous coordinates of the points
X = X';

% the camera matrix
K = [-3933.363636 0       2143.5
     0        3933.363636 1423.5
     0            0            1];
% K = [mx 0 0; 0 -mx 0; 0 0 1]*[c 0 px; 0 c py; 0 0 1]

% the rotation matrix
R = [-0.82757075  -0.557530412  0.065471324
      0.549857636 -0.82856977  -0.10549273
      0.113062965 -0.05130279   0.99226246];

% the translation vector (X~0)
t = [512980.9951 5427701.527 514.79429];


t = t';
% the exterior orientation matrix of the camera
R_t = [R -R*t];

% the projection matrix
P = K*R_t; 

% the image coord.
x = P*X;
% normalizing the coord. to be homogen.
x_x = x(1,:)./x(3,:);
x_y = x(2,:)./x(3,:);
x_z = x(3,:)./x(3,:);

x = [x_x;x_y;x_z];

%%
% computting the difference between the image points
diff_x = im_point_x - x_x
diff_y = im_point_y - x_y


figure

im = imread('R0020851.jpg');
hold on
imshow(im)
hold on
plot (x_x,x_y,'r+','MarkerSize',12);
text(x_x,x_y,label, 'VerticalAlignment','bottom','HorizontalAlignment','right','color','g');
hold on
plot (im_point_x,im_point_y,'bo','MarkerSize',12);
title('Projecting points on Image')
