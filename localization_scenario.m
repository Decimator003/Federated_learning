% Main script: federated_localization.m
clc; clear; close all;

%% Initialize Parameters
% Grid and system parameters
N = 7;                     % Grid size (NxN)
NP = ((N-1)*5)^2;         % Number of pixels
v_d = 0.01;               % Distance from focus to vertex
beta = 0.75;              % Regularization parameter
mu = 2.35;                % Path loss exponent

% Federated Learning parameters
num_nodes = 4;            % Number of distributed nodes
num_rounds = 10;          % Number of federated learning rounds
local_epochs = 5;         % Local training epochs per round

% Measurement parameters
SNRdb = 10;               % Signal-to-noise ratio in dB
Ts = 0.1;                % Sampling period
sigma = 0.4;              % Pixel variance
delta = 2.05;             % Pixel correlation constant

%% Generate Node Positions
[L, node_positions] = generate_node_positions(N);

%% Generate Object Position (Target)
ov = [2 3 1 1];  % Object parameters [x y width height]
[o1, o2] = generate_object(ov);

%% Split Data Among Nodes
node_data = split_data(L, num_nodes);

%% Initialize Global Variables
global_weights = ones(num_nodes, 1) / num_nodes;  % Equal node weights
pixel_grid = generate_pixel_grid(N, NP);
[pixmat, pixm_inv] = generate_correlation_matrix(pixel_grid, sigma, delta);

%% Main Federated Learning Loop
global_model = [];
local_models = cell(num_nodes, 1);

for round = 1:num_rounds
    fprintf('Processing Federated Learning Round %d/%d\n', round, num_rounds);
    
    % Local Computation at Each Node
    parfor node = 1:num_nodes
        local_data = node_data{node};
        [W_M_local, yq_local] = compute_local_measurements(local_data, ov, o1, o2, SNRdb, Ts);
        
        % Local Training
        local_model = train_local_model(local_data, W_M_local, yq_local, pixm_inv, beta, local_epochs);
        local_models{node} = local_model;
    end
    
    % Global Aggregation
    global_model = federated_averaging(local_models, global_weights);
    
    % Visualize Current Round Results
    visualize_results(global_model, N, round);
end

%% Supporting Functions

function [L, node_positions] = generate_node_positions(N)
    % Generate grid points
    x = linspace(1, N, N);
    y = linspace(1, N, N);
    
    % Initialize arrays
    l1 = zeros(N, 2);
    l2 = zeros(N-2, 2);
    l3 = zeros(N, 2);
    l4 = zeros(N-2, 2);
    
    % Generate positions
    for i = 1:N
        l1(i,:) = [x(1), y(i)];
        l3(i,:) = [x(N), y(N-i+1)];
    end
    
    for j = 1:N-2
        l2(j,:) = [x(N-j), y(1)];
        l4(j,:) = [x(j+1), y(N)];
    end
    
    L = [l1; l4; l3; l2];
    node_positions = L;
end

function [o1, o2] = generate_object(ov)
    o1 = [ov(1) ov(1) ov(1)+ov(3) ov(1)+ov(3)];
    o2 = [ov(2) ov(2)+ov(4) ov(2)+ov(4) ov(2)];
end

function node_data = split_data(L, num_nodes)
    total_links = size(L, 1);
    links_per_node = floor(total_links / num_nodes);
    node_data = cell(num_nodes, 1);
    
    for i = 1:num_nodes
        start_idx = (i-1) * links_per_node + 1;
        if i == num_nodes
            end_idx = total_links;
        else
            end_idx = i * links_per_node;
        end
        node_data{i} = L(start_idx:end_idx, :);
    end
end

function pixel_grid = generate_pixel_grid(N, NP)
    pixel_grid = zeros(NP, 2);
    idx = 1;
    for i = 1:sqrt(NP)
        for j = 1:sqrt(NP)
            pixel_grid(idx,:) = [j*0.2, i*0.2];
            idx = idx + 1;
        end
    end
end

function [pixmat, pixm_inv] = generate_correlation_matrix(pixel_grid, sigma, delta)
    NP = size(pixel_grid, 1);
    dist_matrix = zeros(NP, NP);
    
    % Compute distance matrix
    for i = 1:NP
        for j = 1:NP
            dist_matrix(i,j) = norm(pixel_grid(i,:) - pixel_grid(j,:));
        end
    end
    
    % Generate correlation matrix
    fact = sigma/delta;
    pixmat = fact * exp(-dist_matrix/delta);
    pixm_inv = pinv(sqrtm(pixmat));
end

function [W_M, yq] = compute_local_measurements(node_data, ov, o1, o2, SNRdb, Ts)
    num_links = size(node_data, 1);
    NP = ((ov(1)+ov(3)-1)*5)^2;
    
    % Initialize weight matrix
    W_M = zeros(num_links, NP);
    
    % Initialize measurements
    y = zeros(num_links, 1);
    
    % Compute measurements and weights for each link
    for i = 1:num_links
        tx = node_data(i,1:2);
        rx = node_data(i,3:4);
        
        % Compute intersection with object
        dist = compute_intersection(tx, rx, o1, o2);
        
        % Add noise
        noise = sqrt(1/(Ts*(10^(SNRdb/10)))) * randn();
        y(i) = dist + noise;
        
        % Compute weights
        W_M(i,:) = compute_link_weights(tx, rx, NP);
    end
    
    % Quantize measurements
    q = quantizer('fixed', 'convergent', 'wrap', [8 2]);
    yq = quantize(q, y);
end

function dist = compute_intersection(tx, rx, o1, o2)
    % Compute line-polygon intersection
    line_vec = rx - tx;
    t = linspace(0, 1, 100);
    
    points = zeros(length(t), 2);
    for i = 1:length(t)
        points(i,:) = tx + t(i)*line_vec;
    end
    
    % Check if points are inside object
    in_obj = inpolygon(points(:,1), points(:,2), o1, o2);
    
    if any(in_obj)
        % Compute intersection distance
        intersect_points = points(in_obj,:);
        first_intersect = intersect_points(1,:);
        dist = norm(first_intersect - tx);
    else
        dist = 0;
    end
end

function weights = compute_link_weights(tx, rx, NP)
    weights = zeros(1, NP);
    c = norm(rx - tx)/2;  % Half distance
    
    % Compute ellipse parameters
    center = (tx + rx)/2;
    major_axis = norm(rx - tx)/2;
    minor_axis = sqrt(major_axis^2 - c^2);
    angle = atan2(rx(2)-tx(2), rx(1)-tx(1));
    
    % Generate ellipse points
    t = linspace(0, 2*pi, NP);
    ellipse_x = center(1) + major_axis*cos(t)*cos(angle) - minor_axis*sin(t)*sin(angle);
    ellipse_y = center(2) + major_axis*cos(t)*sin(angle) + minor_axis*sin(t)*cos(angle);
    
    % Compute weights based on distance to ellipse
    for i = 1:NP
        dist_to_ellipse = min(sqrt((ellipse_x-tx(1)).^2 + (ellipse_y-tx(2)).^2));
        weights(i) = exp(-dist_to_ellipse);
    end
    
    % Normalize weights
    weights = weights / sum(weights);
end

function local_model = train_local_model(local_data, W_M, yq, pixm_inv, beta, epochs)
    for epoch = 1:epochs
        % Compute local gradients
        W_MT = W_M';
        W_MF = W_MT * W_M;
        
        % Update local model with regularization
        W_regco = pinv(W_MF + beta*pixm_inv);
        local_model = W_regco * W_MT * yq;
    end
end

function global_model = federated_averaging(local_models, weights)
    model_size = size(local_models{1}, 1);
    global_model = zeros(model_size, 1);
    
    for i = 1:length(local_models)
        global_model = global_model + weights(i) * local_models{i};
    end
end

function visualize_results(model, N, round)
    % Reshape model to image
    I_FL = reshape(model, sqrt(length(model)), []);
    I_FL = I_FL / max(max(I_FL));
    
    % Create figure
    figure(round);
    x1 = [1 N];
    y1 = [1 N];
    imagesc(x1, y1, I_FL);
    colormap('jet');
    colorbar;
    
    % Set properties
    set(gca, 'YDir', 'normal');
    xlabel('Distance X(meters)');
    ylabel('Distance Y(meters)');
    title(sprintf('Federated Learning Round %d - Localization Result', round));
    
    drawnow;
end