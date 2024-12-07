function [C] = triad(a,b)
% function [C] = triad(a,b)
% a = (3 x 2 ) 
% b = (3 x 2 )
% 
% find the rotation matrix using TRIAD algorithm
%   a1 = C b1 
%   a2 = C b2

foo1 = size(a);
foo2 = size(b);

if ( (foo1(1) ~= 3) || (foo1(2) ~= 2) || (foo2(1) ~= 3) || (foo2(2) ~= 2) )
    disp('The dimension is wrong (TRIAD method)');
    return;
end

a1 = a(:,1);
a2 = a(:,2);
b1 = b(:,1);
b2 = b(:,2);

b1 = b1 / norm(b1);
a1 = a1 / norm(a1);
b2 = b2 / norm(b2);
a2 = a2 / norm(a2);

r1 = b1;
foo = cross(b1,b2);
r2 = foo / norm(foo);
r3 = cross(b1,r2);

s1 = a1;
foo = cross(a1,a2);
s2 = foo / norm(foo);
s3 = cross(a1,s2);

C = s1*r1' + s2*r2' + s3*r3';