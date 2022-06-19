function [x, iteration_count, relative_residuals, time] = Residual_Norm_Steepest_Descent(A, x0, b, tolerance, max_iterations)

iteration_count = 0;
x = x0;
relative_residual = 1;
relative_residuals = zeros(1, max_iterations);
time = zeros(1, max_iterations);

r = b - A*x;

tic;
while tolerance < relative_residual
    if iteration_count > max_iterations
        error('max iteration limit is exceeded');
    end

    v = A'*r;           %% we have 2 matrix vector 
    z = A*v;            %% operations here
    alpha = norm(v)^2 / norm(z)^2;
    x = x + alpha*v;
    r = r - alpha*z;

    relative_residual = norm(b - A*x) / norm(b - A*x0);
    iteration_count = iteration_count + 1;
    relative_residuals(iteration_count) = relative_residual;
    time(iteration_count) = toc;
end