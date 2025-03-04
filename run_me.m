%% Design 1
ons = (3:5:5*26)'; % every 5 seconds a new trial sart
execution = max(ons) - min(ons); % subjects do 135s of experimental trials
rest_duration = 12; % subject take 12s of rest between every 2 execution period
num_blocks = 4; % they do experimental trials in four mini blocks.
block_offsets = (0:num_blocks-1) * (execution + rest_duration);
% Expand onset times for all mini blocks:
onsets = [];
for i = 1:num_blocks
    onsets = [onsets; ons + block_offsets(i)];
end

nscan = 549; % total 549 volumes recorded

onsets = [onsets(1:52),onsets(53:end)];

%% Design 2


%% Run Simulation
% cannonical hrf: 
hrf_params_canon = [6 16 1 1 6 0 32];
% an example young man hrf: 
hrf_params_youngMan = [4.5 11 1 1 6 0 32];

% SIMULATE DATA:
Y = simulate_GLM('simulate_GLM', nscan, onsets, hrf_params_youngMan, ...
                'add_noise', 1, 'mu', 0, 'sigma', 0.003);

% Design Matrix:
X = simulate_GLM('simulate_GLM', nscan, onsets, hrf_params_canon);

% GLM:
beta = (X' * X)^-1 * X' * Y;
Y_pred = X * beta;

figure;
hrf_y = spm_hrf(1, hrf_params_youngMan);
hrf_x = spm_hrf(1, hrf_params_canon);
plot(hrf_y, 'k', 'LineWidth', 2); hold on;
plot(hrf_x, '--r', 'LineWidth', 2);
legend('data hrf', 'glm hrf')

figure;
plot(Y_pred, '--r', 'LineWidth', 2); hold on;
plot(Y, 'k', 'LineWidth', 2); 
legend('data', 'glm')


