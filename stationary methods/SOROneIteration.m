function x = SOROneIteration(A, b, x0, omega)

%% initialization
n = length(b);
x = zeros(n,1);
D = diag(diag(A));
U = triu(A - D);
L = tril(A - D);

%% SOR
x = (inv(D + omega*L)) * (((1-omega) * D - omega*U) * x0 + omega*b);
% end