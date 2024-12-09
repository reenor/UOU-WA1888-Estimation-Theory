function [result] = find_interval_lessthan(y,threshold,minlen)
% [result] = find_interval_lessthan(y,threshold,minlen)
% minlen should be even number

N = max(size(y));

hminlen = round(minlen/2);

foo = zeros(size(y));
result = foo;

for i = 1:N
    if ( y(i) < threshold )
        foo(i) = 1;
    else
        foo(i) = 0;
    end
end

if ( minlen == 1 )
  result = foo;
  return;
end

for i = 1:N
    starti = max(1,i-hminlen);
    endi = min(N,i+hminlen);
    if ( sum(foo(starti:endi)) == (endi - starti + 1) )
        result(i) = 1;
    end
end

