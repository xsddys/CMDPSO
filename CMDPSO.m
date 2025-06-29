% our paper Jianyu Li, Zhida He, et al. "Coevolutionary Multi-objective Discrete Particle Swarm Optimization for Gateway Placement Optimization Problem" ICACI 2025
% CMDPSO Set velocity_type to 'bi-velocity' to activate CMDPSO dual-velocity mechanism, 
% use 'sigmoid' to activate sigmoid velocity mapping mechanism.

clc;
clear;
close all;

global NA Tmax na T DIMENSION OBJECTIVE_NUM PARTICLE_NUM c1 c2 c3 obj_bunds gBest Swarms archives pmutate mu pexplore D_min archives_accumulate velocity_type func_call

NA = 100;  % Archive size limit
Tmax = 100;  % Maximum iteration count
na = 0; % Current archive position
T = 1; % Current iteration
DIMENSION = 250; % Variable dimension
OBJECTIVE_NUM = 2; % Number of objectives
PARTICLE_NUM = 10; % Number of particles
pmutate = 1;
mu = 0.02;
c1 = 2;
c2 = 2;
c3 = 1.5;
pexplore = 0;
penality = 2;
D_min = 3;
velocity_type = 'bi-velocity'; % Options: 'bi-velocity' or 'sigmoid'
func_call = 0; % Count of fitness function calls
global result_matrix gateway_matrix;

% get the population and gateway data
result_matrix = population_print();
gateway_matrix = test_worldmap();

single_time();

%% this function used to run multiple time inpendently
function run_multiple_times()
    global archives archives_accumulate na func_call
    NUM_RUNS = 20;  % Set number of independent runs
    func_call = 0;
    archives_accumulate = cell(1, NUM_RUNS);  % Initialize a cell array to store archives of each run
    
    for run_id = 1:NUM_RUNS
        fprintf('Starting run %d...\n', run_id);
        na = 0;
        
        % Call main function
        single_time();
        
        % After each run, save the current archives into archives_accumulate
        archives_accumulate{run_id} = cell(1, na);  % Allocate space for current run's archives
        for i = 1:na
            archives_accumulate{run_id}{i} = archives{i};  % Copy each non-dominated solution
        end
        
        fprintf('Run %d completed\n', run_id);
    end
    
    fprintf('Run %d completed, CalFitness function called %d times\n', run_id, func_call);
    % Display information after completion
    fprintf('All %d runs completed\n', NUM_RUNS);
end
%% The main code
% Define structures
function empty_particle = create_empty_particle()
    empty_particle.v = [];
    empty_particle.pos = [];
    empty_particle.fitness = [];
    empty_particle.pBest = [];
    empty_particle.pBest_fitness = [];
    empty_particle.self_archive_pos = [];
    empty_particle.self_archive_fitness=[];
end

function empty_archive = create_empty_archive()
    empty_archive.pos = [];
    empty_archive.fitness = [];
end

function empty_gbest = create_empty_gbest()
    empty_gbest.pos = [];
    empty_gbest.fitness = 0;
end

function empty_boundary = create_empty_boundary()
    empty_boundary.pos_max = [];
    empty_boundary.pos_min = [];
    empty_boundary.v_max = [];
    empty_boundary.v_min = [];
end

% Initialize parameters
function initialize()
    
    global DIMENSION obj_bunds PARTICLE_NUM OBJECTIVE_NUM Swarms gBest archives NA result_matrix gateway_matrix D_min velocity_type func_call
  
    % Initialize boundaries and particles
    obj_bunds = arrayfun(@(~) create_empty_boundary(), 1:DIMENSION, 'UniformOutput', false);
    for i = 1:DIMENSION
        obj_bunds{i}.pos_max = 1.0;
        obj_bunds{i}.pos_min = 0.0;
        
        if strcmp(velocity_type, 'bi-velocity')
            obj_bunds{i}.v_max = 0.2 * (obj_bunds{i}.pos_max - obj_bunds{i}.pos_min);
            obj_bunds{i}.v_min = -obj_bunds{i}.v_max;
        else % sigmoid
            obj_bunds{i}.v_max = 5;
            obj_bunds{i}.v_min = -5;
        end
    end
    
    archives = repmat(create_empty_archive(), NA, 1);
    Swarms = arrayfun(@(x) repmat(create_empty_particle(), PARTICLE_NUM, 1), 1:OBJECTIVE_NUM, 'UniformOutput', false);
    gBest = arrayfun(@(~) create_empty_gbest(), 1:OBJECTIVE_NUM, 'UniformOutput', false);
    
    % Initialize particles
    for i = 1:OBJECTIVE_NUM
        gBest{i}.fitness = inf;  % Initialize global best fitness as infinity
        for j = 1:PARTICLE_NUM
            
            Swarms{i}(j).pos = (rand(1,DIMENSION)>0.5);  % Initialize with random 0-1 string
            
            % Initialize velocity based on velocity_type
            if strcmp(velocity_type, 'bi-velocity')
                Swarms{i}(j).v = rand(DIMENSION,2);  % Dual velocity
            else % sigmoid
                Swarms{i}(j).v = rand(1,DIMENSION)*10-5;  % One-dimensional velocity
            end
            
            % Calculate fitness
            Swarms{i}(j).fitness = CalFitness(Swarms{i}(j).pos, result_matrix, gateway_matrix, D_min);
            func_call = func_call + 1;
            
            Swarms{i}(j).pBest = Swarms{i}(j).pos;  % Initialize personal best
            Swarms{i}(j).pBest_fitness = Swarms{i}(j).fitness;
    
            % Update global best
            if Swarms{i}(j).pBest_fitness < gBest{i}.fitness
                gBest{i}.pos = Swarms{i}(j).pBest;  % Update global best position
                gBest{i}.fitness = Swarms{i}(j).pBest_fitness;
            end
        end
    end

    %disp(Swarms);
    disp(gBest{1}.fitness);
    disp(gBest{2}.fitness);
end

% Calculate fitness，this is our multi object 
function res = CalFitness(x, result_matrix, gateway_matrix, D_min)
    
    res(1)=0;
    res(2)=0;
    pop_num = size(result_matrix, 2);
    gate_num = size(gateway_matrix, 2);
    weight = gateway_matrix(4,:);
    % Initialize variables
    Tj = zeros(1, gate_num); % Store traffic supplied by each IXP
    
    % Calculate distances between regions and enabled IXPs, and allocate traffic
    for i = 1:pop_num
        % Calculate traffic demand for region i
        Ti = result_matrix(3, i)*0.05;
        
        % Calculate distances from region i to all IXPs
        distances = sqrt((result_matrix(1, i) - gateway_matrix(1, :)).^2 + ...
                         (result_matrix(2, i) - gateway_matrix(2, :)).^2);
                     
        % Sort distances to find nearest enabled IXPs
        [~, idx] = sort(distances);
        
        % Allocate traffic to enabled IXPs
        for k = 1:gate_num
            j = idx(k);
            if x(j) == 1 % If IXP is enabled
                % Allocate traffic up to IXP capacity
                allocated = min(Ti/weight(j), gateway_matrix(3, j) - Tj(j));
                Tj(j) = Tj(j) + allocated;
                Ti = Ti - allocated;
                
                % If traffic demand is satisfied, stop allocation
                if Ti <= 0
                    break;
                end
            end
        end
        if Ti>0
           res(1)=res(1)+Ti*0.0005;
           res(2)=res(2)+Ti*0.0005;
        end
           
    end
    % Calculate average traffic T0
    %x
    %size(Tj)
    T0 = sum(Tj(x == 1)) / sum(x);
    
    % Calculate DT
    DT = sum(abs(Tj(x == 1) - T0)) / sum(x);

    % Minimum distance penalty g1(x)
    active_indices = find(x == 1); % Indices of enabled gateways
    active_coords = gateway_matrix(1:2, active_indices); % Coordinates of enabled gateways
    dist_matrix = squareform(pdist(active_coords')); % Calculate pairwise distances
    % Find gateway pairs with distance less than D_min (excluding self-distance)
    [row, col] = find(triu(dist_matrix, 1) < D_min & triu(dist_matrix, 1) > 0); 
    penalty = sum(0.001 * (D_min - dist_matrix(sub2ind(size(dist_matrix), row, col))));
    % Return final result
    res(1) = res(1) + DT / T0 + penalty;
    res(2) = res(2) + sum(x==1)/250 + penalty;
end

% main loap for one time
function single_time()
    
    global OBJECTIVE_NUM PARTICLE_NUM Swarms archives gBest T Tmax na obj_bunds result_matrix gateway_matrix pexplore D_min
    % Start timing
    tic;
    initialize();  % Initialize global variables
    archives = {};  % Initialize archive
    % Start evolution
    while T <= Tmax

        for m= 1:OBJECTIVE_NUM
            for j = 1:PARTICLE_NUM
                ptc = Swarms{m}(j);
                %disp(na)
                % If archive is not empty
                if na ~= 0
                    select = randi([1, na]);  % Randomly select from archive
                    
                    ptc.self_archive_pos = archives{select}.pos;  % Select from archive
                   
                else
                    select = m;
                    while select == m
                        select = randi([1, OBJECTIVE_NUM]);  % Randomly select other objective
                    end
                    ptc.self_archive_pos = gBest{select}.pos;
                    %ptc.self_archive_fitness(m) = gBest{select}.fitness;
                   
                end
                
                % Update velocity and position
                update_V_PSO(ptc, m, j, result_matrix, gateway_matrix, D_min);
               
                %Swarms{m}(j).v = ptc.v;
                %Swarms{m}(j).pos = ptc.pos;
                %Swarms{m}(j).fitness = ptc.fitness;
                %Swarms{m}(j).pBest = ptc.pBest;
                %Swarms{m}(j).pBest = ptc.pBest;
                %Swarms{m}(j).self_archive_pos=ptc.self_archive_pos;

                % Save updated particle
            end
        end
        %r1 = CalFitness(gBest{1}.pos, result_matrix, gateway_matrix, D_min);
        %r2 = CalFitness(gBest{2}.pos, result_matrix, gateway_matrix, D_min);
        %fprintf('gBest1 objective values: %f,%d\n', r1(1), r1(2))
        %fprintf('gBest2 objective values: %f,%d\n', r2(1), r2(2))
        % Print gBest information
        

        UpdateArchive();  % Update archive
        
        % Plot F1 Costs
        figure(1);

        PlotArchives(archives);
        pause(0.01);
        T = T + 1  % Increment iteration count
    end

    % Print final non-dominated solutions
    disp("Run completed")
    T=0;
end


function w=GetWeight(T,Tmax)
    w=0.9-0.5*T/Tmax;
end

% velocity to update
function update_V_PSO(ptc, swarm_num, j, result_matrix, gateway_matrix, D_min)
    global DIMENSION T Tmax Swarms pexplore gBest velocity_type func_call

    w = GetWeight(T, Tmax);  % Get inertia weight
    c1 = 2;
    c2 = 2;
    c3 = 1.5;
    
    if strcmp(velocity_type, 'bi-velocity')
        % Dual velocity update logic  
        for i = 1:DIMENSION
            %c1 - Personal best
            if ptc.pBest(i) == ptc.pos(i)
                vp_0 = 0;
                vp_1 = 0;
            elseif ptc.pBest(i) == 1
                vp_0 = 0;
                vp_1 = 1;
            elseif ptc.pBest(i) == 0
                vp_0 = 1;
                vp_1 = 0;
            end
            
            c11 = c1 * rand();
            vp_0 = min(c11 * vp_0, 1);
            vp_1 = min(c11 * vp_1, 1);
            
            %c2 - Global best
            if gBest{swarm_num}.pos(i) == ptc.pos(i)
                vl_0 = 0;
                vl_1 = 0; 
            elseif gBest{swarm_num}.pos(i) == 1
                vl_0 = 0;
                vl_1 = 1; 
            elseif gBest{swarm_num}.pos(i) == 0
                vl_0 = 1;
                vl_1 = 0;  
            end
            
            c22 = c2 * rand();
            vl_0 = min(c22 * vl_0, 1);
            vl_1 = min(c22 * vl_1, 1);
            
            %c3 - Archive
            if ptc.self_archive_pos(i) == ptc.pos(i)
                va_0 = 0;
                va_1 = 0; 
            elseif ptc.self_archive_pos(i) == 1
                va_0 = 0;
                va_1 = 1; 
            elseif ptc.self_archive_pos(i) == 0
                va_0 = 1;
                va_1 = 0;  
            end
            
            c33 = c3 * rand();
            va_0 = min(c33 * va_0, 1);
            va_1 = min(c33 * va_1, 1);
            
            ptc.v(i, 1) = max([w * ptc.v(i, 1), vp_0, vl_0, va_0]);  % Velocity for position 0
            ptc.v(i, 2) = max([w * ptc.v(i, 2), vp_1, vl_1, va_1]);  % Velocity for position 1
        end   

        
        % Update position based on dual velocity
        for t = 1:DIMENSION
            a = rand();
            if ptc.v(t, 1) > a && ptc.v(t, 2) > a
                ptc.pos(t) = rand() > 0.5;
            elseif ptc.v(t, 1) > a && ptc.v(t, 2) <= a
                ptc.pos(t) = 0;
            elseif ptc.v(t, 1) <= a && ptc.v(t, 2) > a
                ptc.pos(t) = 1;
            elseif ptc.v(t, 1) <= a && ptc.v(t, 2) <= a
                ptc.pos(t) = ptc.pos(t);
            end
        end
    else
        % Sigmoid velocity mapping logic
        ptc.v = w * ptc.v + c1 * rand() * (ptc.pBest - ptc.pos)...
                          + c2 * rand() * (gBest{swarm_num}.pos - ptc.pos)...
                          + c3 * rand() * (ptc.self_archive_pos - ptc.pos);
        
        % Limit velocity range
        for t = 1:DIMENSION
            if ptc.v(t) > 5
                ptc.v(t) = 5;
            end
            if ptc.v(t) < -5
                ptc.v(t) = -5;
            end
        end
        
        % Map velocity using sigmoid function
        sig_v = 1./(1 + exp(-ptc.v));
        
        % Update position
        for t = 1:DIMENSION
            if sig_v(t) > rand()
                ptc.pos(t) = 1;
            else
                ptc.pos(t) = 0;
            end
        end
    end

    % Calculate fitness
    ptc.fitness = CalFitness(ptc.pos, result_matrix, gateway_matrix, D_min);
    func_call = func_call + 1;
    
    % Compare and update personal best
    if ptc.fitness(swarm_num) < ptc.pBest_fitness(swarm_num)
        %fprintf('Update personal best')
        ptc.pBest = ptc.pos;
        ptc.pBest_fitness = ptc.fitness;
    end
    
    % Update global best
    if ptc.pBest_fitness(swarm_num) < gBest{swarm_num}.fitness(swarm_num)
        gBest{swarm_num}.pos = ptc.pBest;
        gBest{swarm_num}.fitness = ptc.pBest_fitness;
        %fprintf('Objective %d global best: %f\n', swarm_num, gBest{swarm_num}.fitness(swarm_num))
    end

    % Explicitly update Swarms
    Swarms{swarm_num}(j).v = ptc.v;
    Swarms{swarm_num}(j).pos = ptc.pos;
    Swarms{swarm_num}(j).fitness = ptc.fitness;
    Swarms{swarm_num}(j).pBest = ptc.pBest;
    Swarms{swarm_num}(j).pBest_fitness = ptc.pBest_fitness;
end


function UpdateArchive()
    global gBest archives obj_bunds OBJECTIVE_NUM Swarms na NA pmutate mu
    % Archive S
    S = {};
    
    % Add each population's personal best to archive S
    for m = 1:OBJECTIVE_NUM
        for j = 1:length(Swarms{m})
  
            
            temp = create_empty_archive();  % Create a new structure
            temp.pos = Swarms{m}(j).pBest;
            temp.fitness = Swarms{m}(j).pBest_fitness;
           
            S{end + 1} = temp;  % Add to S
        end
    end

    % Add previous non-dominated solutions to S
    old_archives = archives;
    for i = 1:na
        S{end + 1} = old_archives{i};  % Add to S
    end
    new_archives = Elitist_learning_strategy(archives, pmutate, mu);


    len = round(na*pmutate);
    for i = 1:len
        S{end + 1} = new_archives{i};  % Add to S
    end
    % Remove duplicate non-dominated solutions
    S = ArchiveFilter(S);
    R = [];

    % Determine non-dominated solutions
    R = Nondominated_solution_determining(S);

    % Process archive R
    size_R = length(R);  % Calculate number of non-dominated solutions
    fprintf('Archive size: %d\n', size_R);  % Output number of non-dominated solutions
    % If exceeding NA, apply density selection
    if size_R > NA
        Density_based_selection(R);
        na = NA;
    else
        na = size_R;
        archives = R(1:na); % Update archive
        
        
    end
    
    
end

function S = ArchiveFilter(S)
    % Remove duplicate non-dominated solutions
    temp = {};
    for i = 1:length(S)
        flag = false;
        for j = 1:length(temp)
            if isequal(temp{j}.pos, S{i}.pos)  % Check if positions are identical
                flag = true;
                break;
            end
        end
        if ~flag
            temp{end + 1} = S{i};  % Add to temporary archive
        end
    end
    S = temp;
end

function flag = dominates(U,W)
    global OBJECTIVE_NUM
    for i=1:OBJECTIVE_NUM
        if U.fitness(i)>W.fitness(i)
            flag=false;
            return;
        end
    end
    flag=true;
end

function R = Nondominated_solution_determining(S)
    % Determine non-dominated solutions
    R = {};
    for i = 1:length(S)
        flag = true;
        for j = 1:length(S)
            if j ~= i && dominates(S{j}, S{i})  % Check dominance relationship
                flag = false;
                break;
            end
        end
        if flag
            R{end + 1} = S{i};  % Add non-dominated solution to R
        end
    end
end

function Density_based_selection(R)
    
    global OBJECTIVE_NUM archives
    L = length(R);
    d = zeros(L, 1);

    % Calculate density vector
    for m = 1:OBJECTIVE_NUM
        max_ftns = GetMaxFtns(m);
        min_ftns = GetMinFtns(m);
        SortRwithObjVal(R, m, d);
        d(1) = inf; d(L) = inf;  % Set ends to infinity

        for i = 2:L-1
            d(i) = d(i) + (R(i + 1).fitness(m) - R(i - 1).fitness(m)) / (max_ftns - min_ftns);
        end
    end

    SortRwithD(R, d);
    archives(1:NA) = R(1:NA);  % Update archive
end

% the modified MELS
function mutated = Elitist_learning_strategy(archives, pmutate, mu)
    global na DIMENSION result_matrix gateway_matrix D_min func_call
    
    len = round(na*pmutate);
    
    % Apply small perturbation to non-dominated solutions
    mutated = num2cell(repmat(create_empty_archive(), len, 1));
    for k=1:len
        i = randi([1, na]);  % Randomly select an archive
        
        archive = archives{i};
        nMu = ceil(mu*DIMENSION);
        j = randsample(DIMENSION, nMu);
        y = archive.pos;
        y(j) = ~y(j);
        mutated{k}.pos = y;
        mutated{k}.fitness = CalFitness(mutated{k}.pos, result_matrix, gateway_matrix, D_min);
        func_call = func_call + 1;
        % Apply small perturbation to dimensions and ensure within bounds
        %E.pos(d) = E.pos(d) + (obj_bunds{d}.pos_max - obj_bunds{d}.pos_min) * normrnd(0, 1);
        %if E.pos(d) > obj_bunds{d}.pos_max
        %   E.pos(d) = obj_bunds{d}.pos_max;
        %elseif E.pos(d) < obj_bunds{d}.pos_min
        %    E.pos(d) = obj_bunds{d}.pos_min;
        %end

        % Recalculate fitness

        % Update the non-dominated solution in the archive
        %archives{i} = E;
    end
end
