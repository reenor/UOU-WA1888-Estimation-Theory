%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Student name: CHUNG QUANG KHANH
% Student ID:   20245360
% Homework #4
% Prof. SUH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clearvars, close all

% measurement: m = 1
n = 2; % state variables

A = [ 1   0;
      0   1 ]; % n x n
d = [ 5   5 ]';
Q = [ 1     0.5;
      0.5   1   ]; % n x n
H = [ 1   0 ]; % m x n
R = 1; % m x m

xh0 = [ 0   0]';
P0  = eye(n);

z = [4  9   16];
N = length(z);

xh_result = zeros(n, N);
P_result  = zeros(n, n, N); % array of N matrices, each matrix is n x n

% Run Kalman filter
xh = xh0;
P  = P0;
for k = 1:N
    % Time update stage
    xh_minus = A*xh + d;             % Prediction of state
    P_minus = A*P*A' + Q;            % Prediction of error covariance
    
    % Measurement update stage
    K = P_minus*H'*inv(H*P_minus*H' + R);           % Computation of Kalman gain
    xh = xh_minus + K*(z(k) - H*xh_minus);          % Computation of estimate
    P = (eye(n) - K*H)*P_minus;                     % Computation of error covariance
    
    % Save results after each round
    xh_result(:, k) = xh;
    P_result(:, :, k) = P;
end

% Show the results
disp('xh_3 ='); disp(xh_result(:, 3));
disp('P_3 = '); disp(P_result(:, :, 3));

% Draw the results
% figure;
% subplot 211;    plot(tt, x, tt, xh_result, 'r');
%                 xlabel('Time'); ylabel('x vs xhat'); legend('x', 'xhat');
% subplot 212;    plot(tt, P_result); xlabel('Time'); ylabel('P(k)');

% figure;
% plot(tt, e_result, tt, sqrt(P_result), 'r', tt, -sqrt(P_result), 'r');
% xlabel('Time'); ylabel('e(k) and SQRT of P(k)');
% legend('e(k)', 'SQRT of P(k)');
