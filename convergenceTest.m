function [] = convergenceTest()
%Convergence test for nBodySolver

tmax = 40;
N = 2;
tracefreq = 1;
masses = [1; 1];
pos0 = [1 1 0;-1 -1 0];
v0 = [-0.375 0 0; 0.375 0 0];

%Call function with level = 6,7,8
[t6, pos6] = nBodySolver(tmax, 6, N, masses, pos0, v0, tracefreq);
[t7, pos7] = nBodySolver(tmax, 7, N, masses, pos0, v0, tracefreq);
[t8, pos8] = nBodySolver(tmax, 8, N, masses, pos0, v0, tracefreq);

%Use x-coordinates for test and downSample the 7 and 8 arrays
pos6 = squeeze(pos6(1,1,:));
pos7 = squeeze(pos7(1,1,1:2:end));
pos8 = squeeze(pos8(1,1,1:4:end));

%Compute differences in grid function from level to level
dpos67 = pos6 - pos7;
dpos78 = 4*(pos7 - pos8);

clf; hold on; plot(t6, dpos67, 'r-.o'); plot(t6, dpos78, 'g-.+');

% plot(t6, squeeze(pos6(1,1,:)), 'r-.o');
% plot(t7, squeeze(pos7(1,1,:)), 'g-.+');
% plot(t8, squeeze(pos8(1,1,:)), 'b-.*');





end

