%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Student name: CHUNG QUANG KHANH
% Student ID:   20245360
% Homework #6
% Prof. SUH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clearvars, close all

%%%%%%%%%%%%%%%%%%%% (2) Find theta_hat using the indirect KF without considering the gyroscope
%%%%%%%%%%%%%%%%%%%%        bias.
load attitude1.mat

r = var(z1);
q = var(z2);

load attitude2.mat

N = length(z1);
thetahat = zeros(N, 1);

P = r;
thetahat(1) = z1(1);

alpha = 2; % numerical integration error??
for k = 2:N
    T = t(k) - t(k-1);
    thetahat(k) = thetahat(k-1) + T * z2(k-1);
    P = P + alpha * q * T^2;
    K = P * inv(P + r);
    delta_theta = K * (z1(k) - thetahat(k));
    thetahat(k) = thetahat(k) + delta_theta;
    P = (1 - K) * P;
end

n = norm( theta - thetahat)^2;
plot(t, theta, 'r', t, thetahat, 'b');
xlabel('Time'); ylabel('Theta and ThetaHat');
legend('Theta', 'ThetaHat');

