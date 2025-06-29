%% This function use BVDPSO to solve weighted single object for GPO problem
% BVDPSO cited from
% M. Shen, Z. H. Zhan, W. N. Chen, Y. J. Gong, J. Zhang, and Y. Li. 2014. Bi-velocity discrete particle 
% swarm optimization and its application to multicast routing problem in communication networks. 
% IEEE Transactions on Industrial Electronics 61, 12 (2014), 7141--7151.

clc;
NP = 30;      % Population size
G = 100;       % Number of iterations
D = 250;        % Dimension of decision variables
c1 = 2;      % Learning factor
c2 = 2;
w_max = 0.8;   % Maximum inertia weight
w_min = 0.4;   % Minimum inertia weight

penality = 2;
% Initialize population individuals
x = (rand(NP,D)>0.5);                        % Generate uniformly distributed binary strings  randn generates normally distributed random numbers
v = rand(NP,D,2);      % Initialize velocity

result_matrix = population_print();
gateway_matrix = test_worldmap();

% Initialize personal best
individual_best = x;       % Historical best of each individual
pbest = zeros(NP, 1);      % Fitness value corresponding to personal best position
global_best_fit = 2000;
for k=1:NP
    pbest(k, 1) = func(individual_best(k, :), result_matrix, gateway_matrix);
    if pbest(k, 1) < global_best_fit
        global_best = individual_best(k, :);
        global_best_fit = pbest(k, 1);
    end
end
 
% Initialize global best
for k=1:NP
    temp = func(individual_best(k, :), result_matrix, gateway_matrix);
    if temp < global_best_fit
        global_best = individual_best(k, :);
        global_best_fit = temp;
    end
end
 
% Start iterations
for gen = 1:G

    disp(gen)
    % Calculate dynamic inertia weight
    w = w_max - (w_max-w_min) * gen / G;
   for k=1:NP
        for j=1:D
            % Calculate c1*r1_j*(pij - xij)
            if individual_best(k,j)==x(k,j)
                vp_0 = 0;
                vp_1 = 0;
            elseif individual_best(k,j)==1
                vp_0 = 0;
                vp_1 = 1;
            elseif individual_best(k,j)==0
                vp_0 = 1;
                vp_1 = 0;
            end
    
            c11 = c1 * rand();
            vp_0 = min(c11 * vp_0,1);
            vp_1 = min(c11 * vp_1,1);
            % Calculate c2*r2_j*(lij - xij)
            if global_best(j) == x(k, j)
                vl_0 = 0;
                vl_1 = 0; 
            elseif global_best(j) == 1
                vl_0 = 0;
                vl_1 = 1; 
            elseif global_best(j) == 0
                vl_0 = 1;
                vl_1 = 0;  
            end
    
            c22 = c2 * rand();
            vl_0 = min(c22 * vl_0,1);
            vl_1 = min(c22 * vl_1,1);
            % Step 3: Determine the final velocity
            v(k, j, 1) = max([w * v(k, j, 1), vp_0, vl_0]);  % Velocity for position 0
            v(k, j, 2) = max([w * v(k, j, 2), vp_1, vl_1]);  % Velocity for position 1
        end
        % Update particle position
        for t=1:D
            a=rand();
            if v(k, t ,1)>a&&v(k, t ,2)>a
                x(k, t) = rand()>0.5;
            elseif v(k, t ,1)>a&&v(k, t ,2)<=a
                x(k, t) = 0;
            elseif v(k, t ,1)<=a&&v(k, t ,2)>a
                x(k, t) = 1;
            elseif v(k, t ,1)<=a&&v(k, t ,2)<=a
                x(k, t) = x(k, t);
            end
        end
    end

    % Calculate personal historical best and global best
    % Personal historical best
    for k=1:NP
        old_fitness = func(individual_best(k, :), result_matrix, gateway_matrix);
        new_fitness = func(x(k, :), result_matrix, gateway_matrix);
        if new_fitness < old_fitness
            individual_best(k, :) = x(k, :);
            pbest(k, 1) = new_fitness;
        end
    end
    % Global best
    for k=1:NP
        temp = func(individual_best(k, :), result_matrix, gateway_matrix);
        if temp < global_best_fit
            global_best = individual_best(k, :);
            global_best_fit = temp;
        end
    end
    disp(sum(global_best==1))
    disp(global_best_fit-0.03*sum(global_best==1))
    global_optimal(gen) = global_best_fit;      % Record in each iteration
    
end

global_best_verse= global_best';

figure(1)
plot(global_optimal);

% Create figure window
figure;
% Read Excel file
data = readtable('aggregated_data_with_lat_lon3.xlsx');

% Extract lat, lon, avg columns and convert to numeric type
lat = data.lat;
lon = data.lon;
avg = data.avg;

lon=lon(global_best_verse==1);
lat=lat(global_best_verse==1);
avg=avg(global_best_verse==1);
% Set figure window size to ensure 2:1 ratio
set(gcf, 'Position', [100, 80, 1200, 600]); % e.g. 1200x600 pixels

% Draw world map
axesm('MapProjection', 'mercator', 'Frame', 'off', 'Grid', 'off', ...
      'MapLatLimit', [-60 90], 'MapLonLimit', [-170 170]); % Limit displayed latitude range

% Set background color to white
set(gcf, 'Color', 'w');

% Add world coastline data
load coastlines
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'k'); % Draw black coastlines

% Get min and max values of avg
avg_min = min(avg);
avg_max = max(avg);


% Replace zero values in avg with a small positive number to avoid log10(0) issues
avg(avg == 0) = 1e-10;

% Scale using log10
avg_log = log10(avg);

% Plot
scatterm(lat, lon, 20, avg_log, 's', 'filled');

% Check caxis range to ensure min and max are not equal
clim_min = min(avg_log);
clim_max = max(avg_log);

if clim_min == clim_max
    clim_min = clim_min - 1; % If equal, adjust the range
    clim_max = clim_max + 1;
end

% Set color bar range
caxis([clim_min clim_max]);

% Set color bar
h = colorbar;
ylabel(h, 'Log10 of Deviation');

% Set axis labels
xlabel('Longitude (°)','FontSize', 18);
ylabel('Latitude (°)','FontSize', 18);

% Add title
title('IXP Locations', 'FontSize', 18);


% Adjust map display range and scale
axis tight;
axis equal;

% Add title
title('IXP Locations', 'FontSize', 18);

function res = func(x, result_matrix, gateway_matrix)
    res=0;
    pop_num = size(result_matrix, 2);
    gate_num = size(gateway_matrix, 2);
    
    % Initialize variables
    Tj = zeros(1, gate_num); % Store the supply traffic for each IXP
    
    % Calculate distances between regions and enabled IXPs, and allocate traffic
    for i = 1:pop_num
        % Calculate traffic demand for region i
        Ti = result_matrix(3, i)*0.1;
        
        % Calculate distances from region i to all IXPs
        distances = sqrt((result_matrix(1, i) - gateway_matrix(1, :)).^2 + ...
                         (result_matrix(2, i) - gateway_matrix(2, :)).^2);
                     
        % Sort distances to find nearest enabled IXPs
        [sorted_distances, idx] = sort(distances);
        
        % Allocate traffic to enabled IXPs
        for k = 1:gate_num
            j = idx(k);
            if x(j) == 1 % If IXP is enabled
                % Allocate traffic until IXP reaches its capacity
                allocated = min(Ti, gateway_matrix(3, j) - Tj(j));
                Tj(j) = Tj(j) + allocated;
                Ti = Ti - allocated;
                
                % If traffic demand is satisfied, stop allocation
                if Ti <= 0
                    break;
                end
            end
        end
        if Ti>0
           res=res+Ti*0.05;
        end
           
    end
    % Calculate average traffic T0
    T0 = sum(Tj(x == 1)) / sum(x);
    
    % Calculate DT
    DT = sum(abs(Tj(x == 1) - T0)) / sum(x);
    
    % Weighted multi object to single object
    res =res+ DT / T0;
    res =res+0.03*(sum(x==1));
end
 
