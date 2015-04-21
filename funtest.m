
function [y, dy] = funtest(x)
    y = exp(-x'*x);
    dy = -2*y*x;
end