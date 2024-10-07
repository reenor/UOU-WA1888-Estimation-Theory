%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Student name: CHUNG QUANG KHANH
% Student ID:   20245360
% Homework #3
% Prof. SUH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clearvars, close all

A = 0.8;
H = 1;
T = 0.02;
N = 100; 
Q = 0.01;
R = 1;
x0 = 0;
P0 = 1;

tt = 0:T:T*(N-1);

% Create 100 samples
w = sqrtm(Q) * randn(N,1);
v = sqrtm(R) * randn(N,1);

x = zeros(N,1);
z = zeros(N,1);

x(1) = sqrtm(P0) * randn;
for i = 2:N
    x(i) = A*x(i-1) + w(i-1);   % True values
    z(i) = H*x(i) + v(i);       % Measurement samples
end

% Run Kalman filter
xh_result = zeros(1, N);
P_result = zeros(1, N);
e_result = zeros(1, N);
xh = x0;
P = P0;
for k = 1:N
    % Time update stage
    xh_minus = A*xh;                % Prediction of state
    P_minus = A*A*P + Q;            % Prediction of error covariance
    
    % Measurement update stage
    K = H*P_minus*inv(H*H*P_minus + R);             % Computation of Kalman gain    
    xh = xh_minus + K*(z(k) - H*xh_minus);          % Computation of estimate
    P = (1 - K*H)*P_minus;                          % Computation of error covariance
    
    % Save result after each round
    xh_result(k) = xh;
    P_result(k) = P;
    e_result(k) = x(k) - xh;
end

% Draw the results
figure;
subplot 211;    plot(tt, x, tt, xh_result, 'r');
                xlabel('Time'); ylabel('x vs xhat'); legend('x', 'xhat');
subplot 212;    plot(tt, P_result); xlabel('Time'); ylabel('P(k)');

figure;
plot(tt, e_result, tt, sqrt(P_result), 'r', tt, -sqrt(P_result), 'r');
xlabel('Time'); ylabel('e(k) and SQRT of P(k)');
legend('e(k)', 'SQRT of P(k)');
