function [r] = draw_ellipse(C,m,K)
% 2D case
% (r - m)' inv(C) (r-m) = K
%

M = 100;

% Calculate the eigenvectors and eigenvalues
covariance = C;
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

phi = angle;
theta_grid = linspace(0,2*pi);

a=sqrt(K)*sqrt(largest_eigenval);
b=sqrt(K)*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) -sin(phi); sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r = R * [ellipse_x_r;ellipse_y_r] + m;

end

