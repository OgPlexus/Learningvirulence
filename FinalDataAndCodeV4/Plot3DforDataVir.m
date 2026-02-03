% Parameters
num_points = 10000; 
mu = 0.8; 
sigma = 0.08; 
iterations = [1, 10, 20]; 
logistic_r = 3.9999; 

% Generate initial points with truncated normal distribution
pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
t = truncate(pd, 0, 1);
initial_points = random(t, num_points, 1);

% Initialize matrix to store data at specified iterations
data = zeros(num_points, length(iterations));
data(:, 1) = initial_points;

% Update points using the logistic map and store snapshots
for iter = 2:max(iterations)
    initial_points = logistic_r * initial_points .* (1 - initial_points);
    if ismember(iter, iterations)
        data(:, find(iter == iterations)) = initial_points;
    end
end

% Plotting
figure;
hold on;
colors = lines(length(iterations));
for i = 1:length(iterations)
    [counts, edges] = histcounts(data(:, i), 100, 'Normalization', 'pdf');
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    n = iterations(i) * ones(size(bin_centers));
    plot3(n, bin_centers, counts, 'Color', colors(i, :), 'LineWidth', 1.5);

    
    min_y = min(bin_centers);
    max_y = max(bin_centers);
    max_z = max(counts);
    plot3([iterations(i) iterations(i)], [min_y max_y], [0 0], 'k', 'LineWidth', 1.5); % Vertical line at the front
    plot3([iterations(i) iterations(i)], [min_y max_y], [max_z max_z], 'k', 'LineWidth', 1.5); % Vertical line at the back
    plot3([iterations(i) iterations(i)], [min_y min_y], [0 max_z], 'k', 'LineWidth', 1.5); % Left vertical line
    plot3([iterations(i) iterations(i)], [max_y max_y], [0 max_z], 'k', 'LineWidth', 1.5); % Right vertical line
end
xlabel('Generation($n$)');
ylabel('Virulence ($\mu$)');
zlabel('density ($p(\mu_n)$)');
title('Distribution of virulence at different generations');
grid on;
view(3);
hold off;
hfig = gcf;
fname = 'muDEnsityplot';
