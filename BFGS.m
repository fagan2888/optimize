
function [x, fmin] = BFGS(f,x0)
    n = length(x0);
    H = eye(n);
    [fmin, dy] = f(x0);
    x = x0; alpha = 1; k = 0;
    while norm(dy)>1e-10
        
        p = - H * dy;
        x = x + alpha * p;
        
        dy0 = dy;
        [fmin, dy] = f(x);
        % calculate H_k+1
        s = alpha * p;
        yk = dy - dy0;
        rh = 1/(yk'*s);
        H = (eye(n) - rh*s*yk')*H*(eye(n) - rh*yk*s') + rh*s*s';
        k = k + 1;
        fprintf('%d dy=%g\n', k, norm(dy));
    end
end