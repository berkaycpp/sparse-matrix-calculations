function x = gaussSeidelOneIteration(A, b, x0, omega)

%% initialization
n = length(b);
x = zeros(n,1);

%% gauss-seidel iteration
for i = 1 : n
    temp = (-A(i, 1:(i-1)) * x(1:(i-1))) - (A(i, (i+1):n) * x0((i+1):n)) + b(i); % -Uxk + b
    x(i) = temp / A(i, i);
end

%% inserting constant omega
x = omega * x + (1 - omega) * x0;
% end