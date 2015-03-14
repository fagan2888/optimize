
function [x, P] = conjugateGradient(A, b, epsilon)
    if ~exist('epsilon','var')
        epsilon = 1e-10;
    end
    x0 = zeros(size(b));
    r0 = A*x0 - b;
    p0 = - r0;
    k = 1;
    rk = r0;
    xk = x0;
    pk = p0;
    n = length(b);
    P = zeros(n,n);
    P(:,1) = p0;
    while norm(rk)>epsilon
        
        ak = rk'*rk/(pk'*A*pk);
        xk = xk + ak*pk;
        rkk = rk+ak*A*pk;
        beta = rkk'*rkk/(rk'* rk);
        pk = -rkk + beta*pk;
        P(:,k+1) = pk;
        k = k+1;
        rk = rkk;
        fprintf('rk=%g\n', norm(rk));
        
    end
    x = xk;
end