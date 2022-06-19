function [x, iteration_count, relative_residuals, time] = Minimum_Residual(A, x0, b, tolerance, max_iterations)

iteration_count = 0;
x = x0;
relative_residual = 1;
relative_residuals = zeros(1, max_iterations);
time = zeros(1, max_iterations);

r = b - A*x;
p = A*r;

tic;
while tolerance < relative_residual
    if iteration_count > max_iterations
        error('max iteration limit is exceeded');
    end

    alpha = r'*p / (p'*p);
    x = x + alpha*r;
    r = r - alpha*p;
    p = A*r;                %% only one matrix vector operation

    relative_residual = norm(b - A*x) / norm(b - A*x0);
    iteration_count = iteration_count + 1;
    relative_residuals(iteration_count) = relative_residual;
    time(iteration_count) = toc;
end