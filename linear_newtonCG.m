% lienar newtonw CG method
function [x, f]= linear_newtonCG(A, b)
    n = length(b);
    x = zeros(n,1); xx = x;
    f = 0;
    % CG-direction
    r = A*x-b;
    p = -r;
    k = 0;
    while k<n || norm(r)<min(0.5, sqrt(norm(b)))*norm(b)
        if p'*A*p<=0
            break;
        else
            x = xx;
        end
        ak = r'*r/(p'*A*p);
        xx = x + ak*p;
        rr = r + ak*A*p;
        beta = rr'*rr/(r'*r);
        p = -rr + beta*p;
        k = k +1;
    end
    
end
