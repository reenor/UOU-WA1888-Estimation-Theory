clc, clearvars, close all

load('3dsim.mat');

N = size(ya,2);
T = 0.01;

qhat = zeros(4,N);
eulerhat = zeros(3,N);
qhat(:,1) = q(:,1);
eulerhat(:,1) = quaternion2euler(q(:,1));
oldOmega = zeros(4,4);

for i = 2:N
    wx = yg(1,i-1);
    wy = yg(1,i-1);
    wz = yg(1,i-1);
    Omega = [ 0 , -wx, -wy, -wz ; wx , 0 , wz , -wy ; wy , -wz, 0, wx ; wz , wy , -wx , 0 ];
    qhat(:,i) = (eye(4) + (1/2) * Omega * T) * qhat(:,i-1);

    qhat(:,i) = (eye(4) + (3/4) * Omega * T - (1/4) * oldOmega * T - (1/8) * T^2 * norm(yg(:,i-1))^2) * qhat(:,i-1);
    qhat(:,i) = qhat(:,i) / norm(qhat(:,i));

    qhat(:,i) = qhat(:,i) / norm(qhat(:,i));
    eulerhat(:,i) = quaternion2euler(qhat(:,i));
end    

plotsensor(eulerhat);
plotsensor(eulertrue - eulerhat);
sum(sum(abs(eulertrue - eulerhat))) %%239.4734