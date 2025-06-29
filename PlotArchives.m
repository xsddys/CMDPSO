% This function used to plot the Archives

function PlotArchives(archives)
    % Initialize array to store fitness values
    CostArray = [];
    
    % Iterate through each struct in archives and extract fitness values
    for i = 1:length(archives)
        CostArray = [CostArray; archives{i}.fitness];
    end
    
    % Plot scatter points
    plot(CostArray(:,1), CostArray(:,2), 'r*', 'MarkerSize', 8);
    % Set axis labels and title
    xlabel('Relative Deviation', 'FontSize', 14);
    ylabel('Activation Ratio of Gateways', 'FontSize', 14);
    title('Comparison of CMPSO', 'FontSize', 16);
    % Show grid
    grid on;
    
    % Add legend
    legend('Pareto Front of CMDPSO');
    
    hold off;

end