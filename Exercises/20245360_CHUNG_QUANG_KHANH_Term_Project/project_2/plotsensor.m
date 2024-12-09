function [ynorm] = plotsensor(y,mode)
% function [ynorm] = plotsensor(y,mode)
% 
%  mode = 1  3 axis
%         2  norm

if nargin < 2
    mode = 1;
end

N = size(y,2);

ynorm = zeros(1,N);
for i = 1:N
   ynorm(i) = norm(y(:,i));
end

tt = 0:0.01:(N-1)*0.01;

if ( mode == 1 )
  figure;
  subplot(3,1,1);
  plot(tt,y(1,:));
  subplot(3,1,2);
  plot(tt,y(2,:));
  subplot(3,1,3);
  plot(tt,y(3,:));
elseif ( mode == 2)
  figure;
  plot(tt,ynorm);
end