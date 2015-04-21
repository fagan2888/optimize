% lienar newtonw CG method
function [x, f]= linear_newtonCG(A, b)
    n = length(b);
    x = zeros(n,1); xx = x;
    f = 0;
    % CG-direction
    k = 0;
    while k<n 
        p = conjugateGradient(A, b, min(0.5, sqrt(norm(b)))*norm(b))
        if p'*A*p<=0
            break;
        else
            x = xx;
        end
        
    end
    
end
