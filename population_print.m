% gpw_v4_population_count_rev11_2020_1_deg.txt
% Cited from NASA, “Gridded population of the world, version 4 (gpwv4): population count, revision 11 | NASA earthdata,". Earthdata,
% Here is the URL:https://earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-popcount-r11-4.11
% this function used to obtain the population density per 1°×1°

function result_matrix = population_print()
    
    
    fileID = fopen('gpw_v4_population_count_rev11_2020_1_deg.txt', 'r');
    
    % skip the first six lines
    for i = 1:6
        fgetl(fileID);
    end
    
    % Initialize matrix
    population_matrix = zeros(180, 360);
    
    % Read data line by line and fill the matrix
    for i = 1:180
        line = fgetl(fileID); % 读取一行
        data = str2num(line); % 将字符串转化为数值
        population_matrix(i, :) = data; % 填充矩阵
    end
    
    % Close file
    fclose(fileID);
    
    % Display matrix
    %disp(population_matrix);
    
    
    
    
    % Replace non-numeric and infinite data with a very small positive number
    population_matrix(population_matrix == -9999) = 0;
    population_matrix=flipud(population_matrix);
    %disp(population_matrix);
    
    %disp(max(population_matrix));
    
    
    lon = linspace(-180, 180, 360);
    lat = linspace(-90, 90, 180);
    
    % Create longitude and latitude grid
    [Lon, Lat] = meshgrid(lon, lat);
    
    % Flatten population matrix into a column vector
    population_flat = population_matrix(:);
    lon_flat = Lon(:);
    lat_flat = Lat(:);
    % Find indices of non-zero elements in population matrix
    non_zero_indices = population_flat ~= 0;
    
    % Extract non-zero elements
    population_non_zero = population_flat(non_zero_indices);
    lon_non_zero = lon_flat(non_zero_indices);
    lat_non_zero = lat_flat(non_zero_indices);
    
    % Generate 3×N array, where N is the number of non-zero elements
    result_matrix = [lon_non_zero'; lat_non_zero'; population_non_zero'];
    
    % Display the size of the result matrix
    %disp(size(result_matrix));  % 应该显示 [3 N]
    
    % Create world map
    figure;
    axesm('pcarree', 'MapLatLimit', [-90 90], 'MapLonLimit', [-180 180]);
    axis off; % Close axis
    
    % Display population density data
    geoshow(Lat, Lon, log10(population_matrix), 'DisplayType', 'texturemap');
    
    % Add coastlines
    load coastlines;
    geoshow(coastlat, coastlon, 'Color', 'cyan');
    
    % Add color bar, adjust color mapping
    colorbar;
    colormap(jet);
    caxis([3 7]); % Adjust color axis range
    
    % Add title
    title('World Population Density');
    
    % Display graph
    set(gcf, 'Color', 'w');
    colormap(jet);
    caxis([3 7]); % Adjust color axis range
    
    % Add title
    title('World Population Density');
    
    % Display graph
    set(gcf, 'Color', 'w');
end