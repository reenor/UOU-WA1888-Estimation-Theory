% 25.7 cm lifting
% load('lift.mat');

% 5 step walking
% load('walking.mat');

% 50 m walking
 load('longwalking1.mat');
% load('longwalking2.mat');

N = size(ya,2);
ra = 0.005;
rg = 0.001;
T = 0.01;

% zero velocity detection
yanorm = zeros(1,N);
zerovel2 = zeros(1,N);
for i = 1:N
    yanorm(i) = norm(ya(:,i));
    if ( (yanorm(i) < 9.8+0.5) && (yanorm(i) > 9.8 -0.5) )
        zerovel2(i) = 1;
    end
end


zerovel = zerovel2;
M = 10;
for i = M+1:N-M
    if ( sum(zerovel2(i-M:i+M)) == 2*M+1 )
        zerovel(i) = 1;
    else
        zerovel(i) = 0;
    end      
end

% plot(zerovel,'*');
% return;

qhat = zeros(4,N);
vhat = zeros(3,N);
rhat = zeros(3,N);

Q = diag([0.25 * rg,0.25 * rg,0.25 * rg,0,0,0,ra,ra,ra]); 
A = zeros(9,9);
A(4:6,7:9) = eye(3);

qhat(:,1) = quaternionya(ya(:,1));
gtilde = [0 ; 0 ; 9.8];
H = zeros(3,9);
H(:,7:9) = eye(3);
R = 0.001 * eye(3);
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
