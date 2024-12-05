function [q,Pqe] = quaternionya(ya,ra)
% function [q,Pqe] = quaternionya(ya)
% ya : 3x1 accelerometer output
% q : quaternion

if ( nargin < 2 )
  ra = 0;
end

w1 = ya;
v1 = [ 0 ; 0 ; 9.8];
w2 = [ 1 ; 0 ; 0 ];
v2 = w2;

C = triad([w1 w2],[v1 v2]);
q = dcm2quaternion(C);

w1 = w1 / norm(w1);

s1 = w1;
s2 = cross(w1,w2);
s3 = cross(s1,s2);
s4 = cross(w2,s2);


s1hat = s1 / norm(s1);
s2hat = s2 / norm(s2);
s4hat = s4 / norm(s4);

ra = ra / (9.8^2);

sigma2 = 0.00000000001;

P = ra * (eye(3) + 1/(norm(s2)*norm(s2)) * ( s4hat * s4hat' + (w1' * w2)^2 * s2hat * s2hat'  )  ) + ...
    sigma2 * ( 1 / ( norm(s2)*norm(s2) ) ) * (eye(3) - s1hat * s1hat');
Pqe = 0.5 * trace(0.25*P) * eye(3) - 0.25*P;
        