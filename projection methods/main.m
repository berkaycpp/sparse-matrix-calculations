close all;
clear all; 

%% matrix market format converter
A = mmread("inputs_A/pde900.mtx");

%% checking spectral radius is less than 1
[sizeA1, sizeA2] = size(A);
if sizeA1 ~= sizeA2
    error('A is not square');
end
factorized = (inv(diag(diag(A)))) * (A - (diag(diag(A))));
spectralRadius = max(abs(eigs(factorized)));
if spectralRadius >= 0.99   % safer
    error('Spectral Radius is = %s',spectralRadius);
end

b = ones(length(A), 1);
x0 = zeros(length(A), 1);
tolerance = 1e-6;
max_iterations = 10000000;
[x, iteration_count, relative_residuals, time] = Steepest_Descent(A, x0, b, tolerance, max_iterations);
subplot(3,1,1)
semilogy(time, relative_residuals, '-o', 'MarkerSize', 4, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
title('Steepest Descent Method for PDE-900 Matrix');
subtitle('Time vs Relative Residual Graph');
xlabel('Time (sec)')
ylabel('Relative Residual')
grid on
ylim([relative_residuals(iteration_count) / 4, relative_residuals(1) * 4])
xlim([0 time(iteration_count)*1.1])
legendText1 = sprintf('Iteration Count = %d', iteration_count);
legend(legendText1);

b = ones(length(A), 1);
x0 = zeros(length(A), 1);
tolerance = 1e-6;
max_iterations = 10000000;
[x, iteration_count, relative_residuals, time] = Minimum_Residual(A, x0, b, tolerance, max_iterations);
subplot(3,1,2)
semilogy(time, relative_residuals, '-o', 'MarkerSize', 4, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow');
title('Minimum Residual Method for PDE-900 Matrix');
subtitle('Time vs Relative Residual Graph');
xlabel('Time (sec)')
ylabel('Relative Residual')
grid on
ylim([relative_residuals(iteration_count) / 4, relative_residuals(1) * 4])
xlim([0 time(iteration_count)*1.1])
legendText1 = sprintf('Iteration Count = %d', iteration_count);
legend(legendText1);

b = ones(length(A), 1);
x0 = zeros(length(A), 1);
tolerance = 1e-6;
max_iterations = 10000000;
[x, iteration_count, relative_residuals, time] = Residual_Norm_Steepest_Descent(A, x0, b, tolerance, max_iterations);
subplot(3,1,3)
semilogy(time, relative_residuals, '-o', 'MarkerSize', 4, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green');
title('Residual Norm Steepest Descent Method for PDE-900 Matrix');
subtitle('Time vs Relative Residual Graph');
xlabel('Time (sec)')
ylabel('Relative Residual')
grid on
ylim([relative_residuals(iteration_count) / 4, relative_residuals(1) * 4])
xlim([0 time(iteration_count)*1.1])
legendText1 = sprintf('Iteration Count = %d', iteration_count);
legend(legendText1);
