function [xf, rhof, x, rho] = IOM_tGMRES(A, b, k)
 
    e1 = eye(k, 1);
    ek = zeros(k, 1); ek(k) = 1;
    theta = norm(b);
    v = b/theta;
    n = length(v);

    [V, H, f] = ArnoldiC(A, k, v);

    beta = norm(f);
    if (k < n)
       [Q,R] = qr([H; beta*ek']);
       rho = abs(Q(1, k+1)*theta);
    else
       [Q, R] = qr(H);
       rho = 0;
    end
    q = Q(1, 1:k)'; 
    y = R(1:k, 1:k)\q;
    x = V*y*theta;

    yf = (H\e1);
    xf = V*yf*theta;
    rhof = beta*abs(yf(k));
    
