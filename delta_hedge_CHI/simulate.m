clear

% Initialize TC state
c = 5; heading = 45;
R = 100; V = 120; % following initial values in CCY
nsample = 4; lambda = 100;

% Initialize parameters
params = initialize_parameters(nsample, c, heading, R, V);

%%
% Run four-path ensemble simulation
[results, path_results] = simulate_tc_ensemble(lambda, params);

% Plot TC tracks
figure;
hold on;

for i = 1:params.nsample
    plot(path_results(i).lon_history, path_results(i).lat_history, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d', i));
end

xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Four Simulated TC Tracks');
%legend('Location', 'best');
grid on;
axis equal;
hold off;

%%
% Run four-path ensemble simulation
[results_hedge, path_results] = simulate_tc_ensemble_hedge(lambda, params);