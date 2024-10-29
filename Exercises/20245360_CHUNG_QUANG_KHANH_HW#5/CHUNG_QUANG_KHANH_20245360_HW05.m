%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Student name: CHUNG QUANG KHANH
% Student ID:   20245360
% Homework #5
% Prof. SUH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clearvars, close all

n = 3; % state variables
m = 1; % measurement

T = 0.01;
F = [ 1   T   T^2/2 ;
      0   1   T ;
      0   0   1       ]; % n x n
Q = 0.05 * [ T^5/20     T^4/8      T^3/6  ;
             T^4/8      T^3/3      T^2/2  ;
             T^3/6      T^2/2      T        ]; % n x n
H = [ 0 0 1 ]; % m x n
R = 0.01; % m x m
I = eye(n);

%%%%%%%%%%%%%%%%%%%% (1) Do simulation to generate xk and zk with initial position 0.
N = 500;
tt = 0:T:T*(N-1);

w = sqrtm(Q) * randn(n, N);
v = sqrtm(R) * randn(m, N);

x = zeros(n, N);
z = zeros(m, N);

% Given that the initial position is 0
x(:, 1) = [0 0 0]';
z(1) = 0;

for k = 2:N
    x(:, k) = F*x(:, k-1) + w(:, k-1);    % True values
    z(:, k) = H*x(:, k) + v(:, k);        % Measurement samples
end

%%%%%%%%%%%%%%%%%%%% (2) With a Kalman filter, find xˆk from zk.
%%%%%%%%%%%%%%%%%%%%     Assume that you know the exact initial position 0.

xh0 = x(:, 1);
P0  = zeros(n, n);

xh_result = zeros(n, N);
P_result  = zeros(n, n, N); % array of N matrices, each matrix is n x n
e_result = zeros(n, N);

% Run Kalman filter
xh = xh0;
P  = P0;
for k = 2:N
    % Time update stage
    xh_minus = F*xh;                 % Prediction of state
    P_minus = F*P*F' + Q;            % Prediction of error covariance
    
    % Measurement update stage
    K = P_minus*H'*inv(H*P_minus*H' + R);           % Computation of Kalman gain
    xh = xh_minus + K*(z(k) - H*xh_minus);          % Computation of estimate
    P = (I - K*H)*P_minus;                          % Computation of error covariance
    
    % Save results after each round
    xh_result(:, k) = xh;
    P_result(:, :, k) = P;
    e_result(:, k) = x(:, k) - xh;
end

%%%%%%%%%%%%%%%%%%%% (3) Draw plot of ek, Pk(1, 1), Pk(2, 2), Pk(3, 3)
figure;
plot(tt, e_result);
xlabel('Time'); ylabel('e(k)');
legend('acc', 'vel', 'pos');

% Extract values from P_result
P_result_1_1 = zeros(1, N);
P_result_2_2 = zeros(1, N);
P_result_3_3 = zeros(1, N);
for k = 1:N
      linearInd1 = sub2ind(size(P_result), 1, 1, k);
      linearInd2 = sub2ind(size(P_result), 2, 2, k);
      linearInd3 = sub2ind(size(P_result), 3, 3, k);
      P_result_1_1(k) = P_result(linearInd1);
      P_result_2_2(k) = P_result(linearInd2);
      P_result_3_3(k) = P_result(linearInd3);
end
figure;
plot(tt, P_result_1_1);
xlabel('Time'); ylabel('Pk(1,1)');
figure;
plot(tt, P_result_2_2);
xlabel('Time'); ylabel('Pk(2,2)');
figure;
plot(tt, P_result_3_3);
xlabel('Time'); ylabel('Pk(3,3)');
