clear all;
close all;

A = mmread('inputs/jpwh_991.mtx');
b = ones(length(A), 1);

%% matlab built-in GMRES function
tol = 1e-12;
restart = linspace(2,75,74); % inner loop, iterations(2), varying parameter
time_gmres = zeros(1, length(restart));
relres_final = zeros(1, length(restart));
max_iteration = 1; % outer loop, iterations(1)
for i = 1:length(restart)
    tic;
    [x_gmres, flag_gmres, relres_final(i)]  = gmres(A, b, restart(i), tol, max_iteration);
    time_gmres(i) = toc;
end

%% IOM and Truncated GMRES implementations
k = linspace(2,75,74); % truncation parameter, varying
time_implemented = zeros(1, length(k));
rho_iom = zeros(1, length(k));
rho_tgmres = zeros(1, length(k));
for j = 1:length(k)
    tic;
    [x_iom, rho_iom(j), x_tgmres, rho_tgmres(j)] = IOM_tGMRES(A, b, k(j));
    time_implemented(j) = toc;
    % no need to calculate IOM and truncated GMRES timings independetly,
    % dominant factor is arnoldi process that is common in both of them.
end

%% figures
subplot(1,2,1)
colororder({'b','r'})
yyaxis left
semilogy(restart,time_gmres, '-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue', 'LineWidth', 1.5);
ylabel('Time (sec)')
yyaxis right
f = semilogy(restart,relres_final, 'r', 'LineWidth', 1.5);
title('Built-in GMRES Method on JPWH-991');
subtitle('Timing and Relative Residual of Varying Restart Parameters');
xlabel('Restart Parameter (m)')
ylabel('Relative Residual')
grid minor

subplot(1,2,2)
colororder({'b','g'})
yyaxis left
semilogy(k,time_implemented, '-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue', 'LineWidth', 1.5);
title('Implemented FOM and Truncated GMRES Method on JPWH-991');
subtitle('Timing and Relative Residual of Varying Truncation Parameters');
xlabel('Truncation Parameter')
ylabel('Time (sec)')
yyaxis right
semilogy(k,rho_iom, '-', 'LineWidth', 1.5);
hold on
f = semilogy(k,rho_tgmres, '--', 'MarkerSize', 3, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green', 'LineWidth', 1.5);
ylabel('Relative Residual')
legend('', 'Rel. Res. of IOM', 'Rel. Res. of tGMRES', 'location', 'best')
hold off
grid minor