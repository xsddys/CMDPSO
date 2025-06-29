% This function used to plot the gateway map, whose function is similar to the test_worldmap.m

data = readtable('aggregated_data_with_lat_lon3.xlsx');

lat = data.lat;
lon = data.lon;
avg = data.traf;
city = data.cit;
contrary=data.ctry;
%xx = pareto_latitude(7,:);
%lon=lon(xx==1);
%lat=lat(xx==1);
%avg=avg(xx==1);
%sum(xx==1)
lat_t=lat';
lon_t=lon';
avg_t=avg';
city_t = data.cit;
contrary_t=contrary';




format long g



figure;

% set the size of the figure window to ensure 2:1 ratio
set(gcf, 'Position', [100, 80, 1200, 600]); % for example, 1200x600 pixels

% plot the world map
axesm('pcarree', 'MapLatLimit', [-90 90], 'MapLonLimit', [-170 170], ...
    'Frame', 'off', 'FFaceColor', 'white');


set(gcf, 'Color', 'w');


load coastlines
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'k', 'LineWidth', 1); % 绘制黑色的海岸线


avg_min = min(avg);
avg_max = max(avg);

disp(avg_min);
disp(avg_max);

% Replace zero values in avg with a very small positive number to avoid log10(0) problem
avg(avg == 0) = 1e-10;

% Use log10 scaling
avg_log = log10(avg);


scatterm(lat, lon, 100, avg_log, 'filled', 'filled');


clim_min = min(avg_log);
clim_max = max(avg_log);

if clim_min == clim_max
    clim_min = clim_min - 1; % 如果相等，调整范围
    clim_max = clim_max + 1;
end

caxis([clim_min clim_max]);


h = colorbar;
ylabel(h, 'Log10 of Capacity','FontSize', 14.4,'FontName','Times New Roman');


xlabel('Longitude (°)','FontSize', 22,'FontName','Times New Roman');
ylabel('Latitude (°)','FontSize', 22,'FontName','Times New Roman');

% Add title
%title('Gateway Locations', 'FontSize', 26,'FontName','Times New Roman','FontWeight','bold');

