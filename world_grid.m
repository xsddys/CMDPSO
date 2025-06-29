% This function used to plot the world grid per 1°×1°, as an example

% Set the latitude and longitude range of Japan
lat_range = [24, 46];  % 纬度范围（日本南北）
lon_range = [122, 150]; % 经度范围（日本东西）


figure;
set(gcf, 'Position', [100, 80, 600, 600]); % 正方形画布
% Plot the map of Japan
axesm('pcarree', 'MapLatLimit', lat_range, 'MapLonLimit', lon_range, ...
    'Frame', 'on', 'FFaceColor', 'white'); % 关闭外部框架

set(gcf, 'Color', 'w');


load coastlines;
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'k', 'LineWidth', 2); % 绘制黑色的海岸线


title('Japan Region with 1° × 1° Grid', 'FontSize', 24,'FontWeight','bold','FontName','Times New Roman');

% plot 1° × 1° grid
hold on;
% set the style of the grid line
setm(gca, 'Grid', 'on', 'GLineStyle', '-', 'Gcolor', [0.5 0.5 0.5], 'GLineWidth', 1);
setm(gca, 'MLineLocation', 1, 'PLineLocation', 1); % 经线和纬线间隔为 1°

% only label the integer longitude and latitude
setm(gca, 'MeridianLabel', 'on', 'ParallelLabel', 'on', ...
    'MLabelLocation', 10, 'PLabelLocation', 10, ... % 仅整十标注
    'MLabelParallel', 'south', 'PLabelMeridian', 'west', ...
    'FontSize', 12, 'FontWeight', 'bold'); % 标签字体大小和粗细


set(gca, 'XColor', 'none', 'YColor', 'none');

hold off;