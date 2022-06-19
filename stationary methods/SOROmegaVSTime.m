clear all
close all

%% matrix market format converter
A = mmread("inputs_A/jpwh_991.mtx");

%% checking spectral radius is less than 1
factorized = (inv(diag(diag(A)))) * (A - (diag(diag(A))));
spectralRadius = max(abs(eigs(factorized)));
if spectralRadius >= 0.99   % safer
    error('Spectral Radius is = %s',spectralRadius);
end

%% initializing the jacobi inputs
omega = linspace(0.7, 1, 30);
time = zeros(1, length(omega));

for i = 1 : length(omega)
    b = ones(length(A), 1);
    x0 = ones(length(A), 1);
    tic;

    %% start of the iteration
    x = SOROneIteration(A, b, x0, omega(i));
    relativeResidual2Norm = norm(b - A*x) / norm(b - A*x0);

    %% continue to iterate
    while relativeResidual2Norm > 10^(-6)
        x = SOROneIteration(A, b, x, omega(i));
        relativeResidual2Norm = norm(b - A*x) / norm(b - A*x0);
    end
    time(i) = toc;
end

plot(omega, time);
% end
