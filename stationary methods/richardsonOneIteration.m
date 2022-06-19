function x = richardsonOneIteration(A, b, x0, omega)

%% initialization
n = length(b);
x = ones(size(x0));
I = eye(n);

%% richardson step
x = (I - omega*A) * x0 + omega*b;
% end
