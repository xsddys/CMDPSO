% PCH, "Internet exchange point datasets," Packet Clearing House, 2024
% Here is the URL: http://www.pch.net/ixp/data

function gateway_matrix = test_worldmap()
    % Read gateway Excel 
    data = readtable('aggregated_data_with_lat_lon3.xlsx');
    
    % get lat, lon, avg
    lat = data.lat;
    lon = data.lon;
    traf = data.traf;
    weight = data.overlap;
    
    gateway_matrix=[lon';lat';traf';weight'];
    
    
    c = colorbar;
    c.Label.String = 'log10 of population';
    
end
