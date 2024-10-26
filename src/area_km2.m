function area_km2 = extract_area(lat, lon)
    % Constants
    R = 6371;  % Earth's radius in kilometers

    % Convert degrees to radians
    lat = deg2rad(lat);
    lon = deg2rad(lon);

    % Number of vertices in the polygon
    n = length(lat);

    % Initialize the sum
    sumExcess = 0;

    % Loop through the vertices
    for i = 1:n
        % Get the next index, wrapping around to 1 if necessary
        j = mod(i, n) + 1;

        % Use spherical excess formula
        sumExcess = sumExcess + (lon(j) - lon(i)) * (2 + sin(lat(i)) + sin(lat(j)));
    end

    % Compute the absolute area using the spherical excess formula
    area_radians2 = abs(sumExcess) / 2;

    % Convert the area to square kilometers
    area_km2 = area_radians2 * R^2;
end

% % Example usage:
% lat = [10, 20, 20, 10, 10];  % Latitude of polygon vertices
% lon = [-100, -100, -90, -90, -100];  % Longitude of polygon vertices
% 
% % Compute the area in km^2
% area = polygon_area_geographic(lat, lon);
% 
% % Display the result
% fprintf('The area of the polygon is %.2f km^2\n', area);