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

% Make onsets into multiple neuronal causes/conditions:
half1 = onsets(1:52);
half2 = onsets(53:end);

half1 = half1(randperm(length(half1)));
half2 = half2(randperm(length(half2)));
onsets = [half1';half2'];

%% Design 2


%% Run Simulation
% an example young man hrf: 
hrf_params_y = [4.5 11 1 1 1 0 32];
% cannonical hrf: 
hrf_params_x = [6 16 1 1 6 0 32];

figure;
hrf_y = spm_hrf(1, hrf_params_y);
hrf_x = spm_hrf(1, hrf_params_x);
plot(hrf_y/max(hrf_y), 'k', 'LineWidth', 2); hold on;
plot(hrf_x/max(hrf_x), '--r', 'LineWidth', 2);
legend('data hrf', 'glm hrf')

% SIMULATE DATA:
Y = simulate_GLM('simulate_GLM', nscan, onsets, hrf_params_y, ...
                'add_noise', 0, 'mu', 0, 'sigma', 0.003);
Y = sum(Y,2);
mu = 0;
sigma = max(Y)/7;
Y = Y + mu + sigma * randn(size(Y));

% Design Matrix:
X = simulate_GLM('simulate_GLM', nscan, onsets, hrf_params_x);

% GLM:
% mean rmv:
X = X - mean(X,1);
Y = Y - mean(Y);

% OLS:
beta = (X' * X)^-1 * X' * Y;
Y_pred = X * beta;

figure;
plot(Y, 'k', 'LineWidth', 2); hold on;
plot(Y_pred, ':r', 'LineWidth', 2); hold on;
legend('data', 'glm')


