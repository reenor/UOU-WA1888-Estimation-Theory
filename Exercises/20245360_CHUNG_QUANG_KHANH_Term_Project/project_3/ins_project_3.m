
clc, clearvars, close all

% Read data
%file_name = '7_680A_D422CD00680A_20241207_205743.csv'; %29.125139410202475
file_name = '7_680A_D422CD00680A_20241207_210050.csv'; %29.030281208085523
[ya, yg, ym, e3, T] = xsens_read(file_name);
%plotsensor(ya);
%plotsensor(yg);
N = size(ya, 2);
tt = 0:T:T*(N-1);

% Zero velocity detection
ya_threshold = 1.6;
ya_length = 12;
yg_threshold = 0.5;
yg_length = 5;

ya_diff_norm = compute_sensor_norm(ya, 2);
ya_zero_vel = find_interval_lessthan(ya_diff_norm, ya_threshold, ya_length);
% figure; plot(tt, ya_diff_norm, '-', Color="#FF0000"); hold on;
%         plot(tt, ya_zero_vel, '*', Color="#0000FF"); hold off

yg_norm = compute_sensor_norm(yg, 1);
yg_zero_vel = find_interval_lessthan(yg_norm, yg_threshold, yg_length);
% figure; plot(tt, yg_norm, '-', Color="#FF00FF"); hold on;
%         plot(tt, yg_zero_vel, '*', Color="#0000FF"); hold off

zero_vel = ya_zero_vel.*yg_zero_vel;
figure; plot(tt, ya_diff_norm, '-', Color="#FF0000"); hold on;
        plot(tt, yg_norm, '-', Color="#FF00FF");
        plot(tt, zero_vel, '*', Color="#0000FF"); hold off
legend('ya_diff_norm','yg_norm','zero_vel')

ra = 0.005; %?
rg = 0.001; %?
qhat = zeros(4,N);
vhat = zeros(3,N);
rhat = zeros(3,N);

Q = diag([0.25 * rg,0.25 * rg,0.25 * rg,0,0,0,ra,ra,ra]);
A = zeros(9,9);
A(4:6,7:9) = eye(3);

qhat(:,1) = quaternionya(ya(:,1));
gtilde = [0 ; 0 ; 9.8];
H = zeros(4,9); % project 1: insert position
H(1,6) = 1; % project 1: position
H(2:4,7:9) = eye(3); % project 1: velocity
R = zeros(4,4); % project 1: insert position
R(1,1) = 0.01; % project 1: position
R(2:4,2:4) = 0.001 * eye(3); % project 1: velocity
P = diag([0.001 0.001 0.001 0 0 0 0 0 0 ]);
oldomega4 = zeros(4,4);
for i = 2:N
    Cq = quaternion2dcm( qhat(:,i-1) );
    
    A(1:3,1:3) = -vec2product(yg(:,i-1));
    A(7:9,1:3) = - 2 * Cq' * vec2product(ya(:,i-1));
    Qd = Q * T + (T^2 / 2) * A * Q + (T^2/2) * Q * A';

    dA = eye(9) + A * T + A * A * T^2 / 2;
    P = dA * P * dA' + Qd;
    omega4 = compute_44(yg(:,i-1));
    qhat(:,i) = ( eye(4) + 0.75 * omega4 * T - 0.25 * oldomega4 * T - (1/6) * norm(yg(:,i-1))^2 * T^2 * eye(4) - (1/24) * omega4 * oldomega4 * T^2  - (1/48)*norm(yg(:,i-1))^2 *omega4 * T^3) * qhat(:,i-1);
    qhat(:,i) = qhat(:,i) / norm(qhat(:,i));
    oldomega4 = omega4;
    vhat(:,i) = vhat(:,i-1) +  0.5* T * ( quaternion2dcm(qhat(:,i-1))' * ya(:,i-1) + quaternion2dcm(qhat(:,i))' * ya(:,i) )  -  T * [ 0 ; 0 ; 9.8];
    rhat(:,i) = rhat(:,i-1) + 0.5* T * (vhat(:,i) + vhat(:,i-1));
    
    if ( zero_vel(i) == 1 )
        z = zeros(3,1) - vhat(:,i);
        z = [0 - rhat(3,i); z]; % project 1: insert position into z
        K = P * H' * inv(H * P * H' + R);
        x = K * z;
    
        P = ( eye(9) - K * H) * P;
        P = 0.5 * (P + P');
        
        vhat(:,i) = vhat(:,i) + x(7:9);
        rhat(:,i) = rhat(:,i) + x(4:6);

        qe = [ 1 ; x(1:3) ];
        qhat(:,i) = quaternionmul(qhat(:,i),qe);
        qhat(:,i)=  qhat(:,i) / norm(qhat(:,i));
    end
end    
    
plotsensor(vhat);
plotsensor(rhat);

figure; plot(rhat(1,:),rhat(2,:)); axis equal
