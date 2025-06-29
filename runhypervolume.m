
% 初始化一个数组来存储每个 hypervolume 的计算结果
hypervolume_results = zeros(5, length(current_pareto_front));  % 5次独立运行，每次运行存储所有 hypervolume

for run_id = 1:5
    % 获取当前运行的 Pareto 前沿
    current_pareto_front = pareto_fronts{run_id};
    
    % 初始化一个数组来存储当前运行的所有 hypervolume
    current_hypervolume = zeros(1, length(current_pareto_front));
    
    for i = 1:length(current_pareto_front)
        % 获取当前的 Pareto 前沿矩阵
        pareto_front_CMPSO = vpa(current_pareto_front{i});
        
        % 处理矩阵，进行归一化
        pareto_front_CMPSO(:, 1) = vpa(pareto_front_CMPSO(:, 1) / 0.9774);
        pareto_front_CMPSO(:, 2) = vpa(pareto_front_CMPSO(:, 2) / 0.3333);
        pareto_front_CMPSO=pareto_front_CMPSO(1:2);
        pareto_front_CMPSO = double(pareto_front_CMPSO);
        % 计算 hypervolume
        current_hypervolume(i) = hypervolume(pareto_front_CMPSO, [1.1, 1.1], 5000);
    end
    
    % 存储当前运行的 hypervolume 结果
    hypervolume_results(run_id, :) = current_hypervolume;  % 存储每个 hypervolume
end
% ... existing code ...

% 保存结果到文件
%save('hypervolume_results.mat', 'hypervolume_results');

