clc, clearvars, close all

% 25.7 cm lifting
%load('lift.mat');

% 5 step walking
%load('walking.mat');

% 50 m walking
% load('longwalking1.mat');
% load('longwalking2.mat');

% Running
load('running.mat');
%plotsensor(ya);

N = size(ya,2);
ra = 0.005;
rg = 0.001;
T = 0.01;
tt = 0:T:T*(N-1);
length_ya = 3;
length_yg = 4;
threshold_ya = 35;
threshold_yg = 8.5;

% zero velocity detection for ya
yanorm = compute_sensor_norm(ya,1);
zero_vel_ya = find_interval_lessthan(yanorm,threshold_ya,length_ya);

% zero velocity detection for yg
ygnorm = compute_sensor_norm(yg,1);
zero_vel_yg = find_interval_lessthan(ygnorm,threshold_yg,length_yg);

% combine ya and yg
zerovel = zero_vel_ya.*zero_vel_yg;
figure; plot(tt,yanorm,'r-',tt,ygnorm,'b-',tt,zerovel,'b*')
legend('yanorm','ygnorm','zerovel')

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
    
    if ( zerovel(i) == 1 )
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
