function x = jacobiOneIteration(A, b, x0, omega)

%% initialization
n = length(b);
x = ones(size(x0));

%% jacobi iteration
for i = 1 : n
     x(i) = b(i);
     for j = 1 : n
         if j ~= i
             x(i) = x(i) - A(i, j) * x0(j);
         end
     end
     x(i) = x(i) / A(i, i);
end

%% inserting constant omega
x = omega * x + (1 - omega) * x0;
% end
