function [C] = dcmfromyaym(ya,ym)
% function [C] = dcmfromyaym(ya,ym)


ya = ya / norm(ya);
ym = ym / norm(ym);

foo = cross(ya,ym);
foo = foo / norm(foo);
C = [ -cross(ya,foo) , foo, ya ];

