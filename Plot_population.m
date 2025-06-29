% This function used to plot the population map, whose function is similar to the population_print.m

fileID = fopen('gpw_v4_population_count_rev11_2020_1_deg.txt', 'r');


for i = 1:6
    fgetl(fileID);
end

% 初始化矩阵
population_matrix = zeros(180, 360);

% 逐行读取数据并填充矩阵
for i = 1:180
    line = fgetl(fileID); % 读取一行
    data = str2num(line); % 将字符串转化为数值
    population_matrix(i, :) = data; % 填充矩阵
end

% 关闭文件
fclose(fileID);

% 处理非数值和无穷大数据，将其替换为一个很小的正数
population_matrix(population_matrix == -9999) = 0;
population_matrix = flipud(population_matrix);

% 设置经纬度
lon = linspace(-180, 180, 360);
lat = linspace(-90, 90, 180);

% 创建经纬度网格
[Lon, Lat] = meshgrid(lon, lat);

% 将人口矩阵拉平为一个列向量
population_flat = population_matrix(:);
lon_flat = Lon(:);
lat_flat = Lat(:);

% 找出人口矩阵中不为零的索引
non_zero_indices = population_flat ~= 0;

% 提取不为零的元素
population_non_zero = population_flat(non_zero_indices);
lon_non_zero = lon_flat(non_zero_indices);
lat_non_zero = lat_flat(non_zero_indices);

% 生成 3×N 数组，其中 N 为不为零的元素个数
result_matrix = [lon_non_zero'; lat_non_zero'; population_non_zero'];

% 显示结果矩阵的大小
disp(size(result_matrix));  % 应该显示 [3 N]

% 创建世界地图并设置图窗属性
figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.6], 'Color', 'white'); % 添加白色背景

% 设置图像分辨率
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 16 9]); % 设置输出图像的尺寸比例为16:9
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'Renderer', 'painters'); % 使用矢量渲染器提高质量

% 使用 'pcarree' 投影绘制地图
axesm('pcarree', 'MapLatLimit', [-90 90], 'MapLonLimit', [-180 180], ...
    'Frame', 'on', 'FFaceColor', 'white', ...
    'FLineWidth', 1); % 增加边框线宽

% 设置坐标轴背景为白色
set(gca, 'Color', 'white');

axis off; % 关闭轴线

% 显示人口密度数据
geoshow(Lat, Lon, log10(population_matrix), 'DisplayType', 'texturemap');

% 添加海岸线并设置颜色
load coastlines;
geoshow(coastlat, coastlon, 'Color', 'white', 'LineWidth', 0.5); % 将海岸线设置为黑色并加粗

% 设置颜色条
colorbar;
colormap(jet);
caxis([3 7]); % 调整颜色轴范围
% 在绘图后添加以下代码
c = colorbar;
c.Label.String = 'log10 of population';
% 添加标题


% 设置轴标签
xlabel('Longitude (°)','FontSize', 18);
ylabel('Latitude (°)','FontSize', 18);
% 添加纬度和经度标签
textm(-80, 0, 'Latitude (°)', 'FontSize', 18, 'HorizontalAlignment', 'center','FontName','Times New Roman');
textm(-80, 0, 'World Population Density', 'FontSize', 18, 'HorizontalAlignment', 'center','FontName','Times New Roman');
textm(0, 190, 'Longitude (°)', 'FontSize', 18, 'Rotation', 90, 'HorizontalAlignment', 'center','FontName','Times New Roman');


% 提高图片分辨率和清晰度
set(gca, 'FontSize', 16); % 设置坐标轴的字体大小

set(gcf, 'PaperPositionMode', 'auto'); % 保持图形的宽高比例

% 在绘图设置的最后添加
set(gca, 'Position', [0.05 0.05 0.9 0.9]); % 调整绘图区域，减少边距
