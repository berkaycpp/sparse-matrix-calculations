clear all;
close all;

A = mmread('inputs/bfw62b.mtx');

%% parameters
k = 25;
X = ones(length(A), k);
max_iter = 100000;
tol = 1e-12;

%% matlab eigs functions for finding largest and smallest eigenpairs
tic;
[V1, D1] = eigs(A, [], k, 'largestabs');
eigs_max_time = toc;
V1 = abs(V1);
tic;
[V2, D2] = eigs(A, [], k, 'smallestabs');
V2 = abs(V2);
eigs_min_time = toc;

%% ssi(subspace iterations) and sii(subspace inverse iterations) algorithms
i = 0;
tic;
while (i < max_iter) & tol < norm(abs(V1) - abs(X))
    Z = A*X;
    [X, R] = qr(Z, 0);
    i = i+1;
end
iteration_ssi = i;
ssivec = abs(X);
for i = 1:k
    ssival_(:, 1) = A * X(:, i);
    ssival(i,i) = ssival_(1, 1) / X(1, i);
end
ssi_time = toc;

X = ones(length(A), k);
i = 0;
tic;
A = inv(A);
while (i < max_iter) & tol < norm(abs(V2) - abs(X))
    Z = A*X;
    [X, R] = qr(Z, 0);
    i = i+1;
end
iteration_sii = i;
siivec = abs(X);
for i = 1:k
    siival_(:, 1) = A * X(:, i);
    siival(i,i) = 1 / (siival_(1, 1) / X(1, i));
end
sii_time = toc;

%% figures
xa = categorical({'eigs()_l_a_r_g_e_s_t Time Spent','SSI Time Spent','eigs()_s_m_a_l_l_e_s_t Time Spent','SII Time Spent'});
xa = reordercats(xa, {'eigs()_l_a_r_g_e_s_t Time Spent','SSI Time Spent','eigs()_s_m_a_l_l_e_s_t Time Spent','SII Time Spent'});
ya = [eigs_max_time, ssi_time, eigs_min_time, sii_time];
b = bar(xa,ya, 0.1);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Time Consumption Comparison');
subtitle('Built-in eigs() with SSI and SII, tolerance value is 1e-12');
ylabel('Time (sec)')
ylim([0, max(ya)+max(ya)/10])
grid minor








