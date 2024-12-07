
function [ya, yg, ym, e3, T] = xsens_read(file_name)
% Read data from a Xsens DOT
%   Input
%       file_name : csv format
%   Output
%       ya : acceleration in m/s^2 unit
%       yg : angular velocity in rad/s unit
%       ym : magnetic field in a.u. unit (1 a.u. ~ 40uT)
%       e3 : euler angles in degree unit
%       T : time interval in second unit

% set the parameters
sample_time_col = 1;

% Read data starting from cell B10
A = readmatrix(file_name, 'Range', 'B10');

e3 = A(:, 2:4)';
ya = A(:, 5:7)';
yg = A(:, 8:10)' * pi / 180; % convert from deg/s -> rad/s
ym = A(:, 11:13)';

t_us = A(:, sample_time_col)'; % The unit of sample time on xsens dot is us
t  = t_us / 10^6; % Convert them to unit of second
T = mean(t(2:end) - t(1:end-1)); % Compute the change in time in second

% todo: check status at the last column 14th

end % function
