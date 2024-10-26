% This script allows to successively map the plastics for all time steps
% and then perform 30 days or spatial averages that are then smoothed before
% extracting the contours

% Start and end dates
date_start = date_num(2008,01,01);
date_end = datenum(2024,07,25);

% Map mesh size and boundaries
dx = 0.08;
gpgp_boundaries = [20 45 -160 -125];
gpgp_boundaries_p = [20-dx/2 45+dx/2 -160-dx/2 -125+dx/2];

% Number of temporal averages
isteps = 30;

% Grids generation
binEdgesX = linspace(gpgp_boundaries(3),gpgp_boundaries(4),(gpgp_boundaries(4)-gpgp_boundaries(3))/dx+1);
binEdgesY = linspace(gpgp_boundaries(1),gpgp_boundaries(2),(gpgp_boundaries(2)-gpgp_boundaries(1))/dx+1);
binEdgesX_p = linspace(gpgp_boundaries_p(3),gpgp_boundaries_p(4),(gpgp_boundaries_p(4)-gpgp_boundaries_p(3))/dx+1);
binEdgesY_p = linspace(gpgp_boundaries_p(1),gpgp_boundaries_p(2),(gpgp_boundaries_p(2)-gpgp_boundaries_p(1))/dx+1);
[YGrid, XGrid] = meshgrid(binEdgesY(1:end),binEdgesX(1:end));

% Allows to avoid recomputing everything
cod_case = 4;

% Initialize counters
iyear = 0;
icurrent = 0;
icurrentyear = 0;
%%
for date = date_start : date_end
    tic
    icurrent = icurrent + 1;
    time(icurrent) = date;
    time_s(icurrent) = (date - datenum(1970,01,01))*3600*24;
    if cod_case <= 0
        datevect = datevec(date);
        if datevect(1) == iyear + 1
            netcdf.close(ncid);
            icurrentyear = 0;
            ncid = netcdf.open(['ADVECTOR_2D_output_' num2str(datevect(1)) '.nc']);
            iyear = datevect(1);
        else
            if icurrent == 1
                ncid = netcdf.open(['ADVECTOR_2D_output_' num2str(datevect(1)) '.nc']);
                release_date = netcdf.getVar(ncid,1);
                iyear = datevect(1);
            else
            end
        end
        if datevect(1)<=2022
            X = netcdf.getVar(ncid,4,[icurrentyear 0],[1 1144442]);
            Y = netcdf.getVar(ncid,5,[icurrentyear 0],[1 1144442]);
        else
            X = netcdf.getVar(ncid,3,[0 icurrentyear],[1144442 1]);
            Y = netcdf.getVar(ncid,4,[0 icurrentyear],[1144442 1]);
        end
        X(X>=180) = X(X>=180) - 360;
        mask1 = (X >= gpgp_boundaries(3)) & (X <= gpgp_boundaries(4)) & (Y >= gpgp_boundaries(1)) & (Y <= gpgp_boundaries(2));
        mask2 = (release_date <= time_s(icurrent));
        if datevect(1)<=2022
            Xtrimmed = X(mask1' & mask2);
            Ytrimmed = Y(mask1' & mask2);
            cut = sum(mask1'.*mask2);
        else
            Xtrimmed = X(mask1 & mask2);
            Ytrimmed = Y(mask1 & mask2);
            cut = sum(mask1.*mask2);
        end        
        if cut == 0
            counts = zeros(size(binEdgesX,2),size(binEdgesY,2));
        else
            [counts, edgesX, edgesY] = histcounts2(Xtrimmed, Ytrimmed, binEdgesX_p, binEdgesY_p);
        end
        save(['raw_heatmaps/heatmap_' num2str(icurrent) '.mat'], 'counts');
        figure(1)
        pcolor(XGrid,YGrid,counts); shading flat;
        title(datestr(date));
        drawnow()
        icurrentyear = icurrentyear + 1;
    end
    toc
end
%% Spatial averaging
if cod_case <= 1
    for k = 1:length(date_start:date_end)
        load(['raw_heatmaps/heatmap_' num2str(k) '.mat']);
        % Compute the block size in fine grid coordinates
        block_size_x = 438 / 71;
        block_size_y = 313 / 51;

        % Loop over the coarse grid
        Zcoarse=zeros(71,51);
        for i = 1:71
            for j = 1:51
                % Define the region in the fine grid that corresponds to this coarse grid cell
                x_start = round((i-1)*block_size_x + 1);
                x_end = min(round(i*block_size_x), 438);  % Ensure we don't exceed the grid size
                y_start = round((j-1)*block_size_y + 1);
                y_end = min(round(j*block_size_y), 313);  % Ensure we don't exceed the grid size

                % Handle the fine grid points that correspond to this coarse grid cell
                fine_grid_block = counts(x_start:x_end, y_start:y_end);

                % Compute the average of the block
                Zcoarse(i,j) = mean(fine_grid_block(:))*(0.5/0.081)^2;
            end
        end
        counts = Zcoarse;
        save(['raw_heatmaps2/heatmap_' num2str(k) '.mat'], 'counts');
    end
end
%%
%% Temporal averaging
if cod_case <= 2

    % Loading projected heatmaps
    for i = 1:length(date_start:date_end)
        load(['raw_heatmaps/heatmap_' num2str(i) '.mat']);
        counts_tab{i} = counts;
    end
    istart = isteps/2+1; iend = length(date_start:date_end) - isteps/2;
    %%
    for i = 1:length(date_start:date_end)
        tic
        i
        if cod_case <=1
            if i ==1
                count_tot = zeros(size(XGrid));
                for j=1:30
                    count_tot = count_tot + counts_tab{i-15+j};
                end
            else
                count_tot = count_tot + counts_tab{i+isteps/2} - counts_tab{i-isteps/2-1};
            end
            %figure(2)
            %pcolor(XGrid,YGrid,count_tot); shading flat;
            %title(datestr(time(i)));
            %drawnow()
            save(['30days_average_heatmaps/heatmap_average' num2str(i) '.mat'], 'count_tot');
        end
        toc
    end
end

%%
for i = length(date_start:date_end)
    tic
    i
    load(['30days_average_heatmaps/heatmap_average' num2str(i) '.mat']);
    counts_tab_average{i} = count_tot;
    toc
end
%%
areaMap = (1.853*60)^2.*(cos(YGrid*pi/180))*dx*dx;
areaMapvect = areaMap(:);
for i = istart:iend
    i
    smoothedCounts = counts_tab_average{i}./areaMap;
    smoothedCounts_vect = smoothedCounts(:);
    [smoothedCount_vect_sorted, I_sorted] = sort(smoothedCounts_vect);
    smoothedCount_vect_sorted_cum_sum = cumsum(smoothedCount_vect_sorted);
    mask1 = smoothedCount_vect_sorted_cum_sum > 0.25*smoothedCount_vect_sorted_cum_sum(end);
    threshold = min(smoothedCount_vect_sorted(mask1));
    mask2 = smoothedCount_vect_sorted_cum_sum > 0.5*smoothedCount_vect_sorted_cum_sum(end);
    threshold = min(smoothedCount_vect_sorted(mask2));
    mask3 = smoothedCount_vect_sorted_cum_sum > 0.75*smoothedCount_vect_sorted_cum_sum(end);
    threshold = min(smoothedCount_vect_sorted(mask3));
    for i_case = 1:3
        if i_case == 1
            mask = mask1;
            mass = 75000000;
        elseif i_case == 2
            mask = mask2;
            mass = 50000000;
        else
            mask = mask3;
            mass = 25000000;
        end
        sum_part = sum(areaMapvect(I_sorted(mask)).*smoothedCounts_vect(I_sorted(mask)));
        fact = mass/sum_part;
        average = fact*mean(smoothedCounts_vect(I_sorted(mask)));
        mediant = fact*median(smoothedCounts_vect(I_sorted(mask)));
        if i_case == 1
            average_1(i) = average;
            median_1(i) = mediant;
        elseif i_case == 2
            average_2(i) = average;
            median_2(i) = mediant;
        else
            average_3(i) = average;
            median_3(i) = mediant;
        end
    end
    gpgp_area1(i)=sum(areaMapvect(I_sorted(mask1)));
    gpgp_area2(i)=sum(areaMapvect(I_sorted(mask2)));
    gpgp_area3(i)=sum(areaMapvect(I_sorted(mask3)));
end
%%
figure
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),gpgp_area1(istart:iend));hold on;
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),gpgp_area2(istart:iend))
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),gpgp_area3(istart:iend))
legend('GPGP','HST','HS');
grid on;
figure
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),average_1(istart:iend));hold on;
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),average_2(istart:iend))
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),average_3(istart:iend))
legend('GPGP','HST','HS');
grid on;
figure
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),median_1(istart:iend));hold on;
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),median_2(istart:iend))
plot(datetime(time(istart:iend),'ConvertFrom','datenum'),median_3(istart:iend))
legend('GPGP','HST','HS');
grid on;
figure(3)
title('Median densities')
figure(2)
title('Averaged densities')
figure(1)
title('GPGP surface')
res_tab = table(datetime(time(istart:iend)','ConvertFrom','datenum'),gpgp_area1(istart:iend)',gpgp_area2(istart:iend)',gpgp_area3(istart:iend)',average_1(istart:iend)',average_2(istart:iend)',average_3(istart:iend)',median_1(istart:iend)',median_2(istart:iend)',median_3(istart:iend)','VariableNames',{'Time','surf_GPGP','surf_HST','surf_HS','mean_density_GPGP','mean_density_HST','mean_density_HS','median_density_GPGP','median_density_HST','median_density_HS'});
writetable(res_tab,'output.csv');
%%
for i = istart : iend
    tic
    if cod_case <=2
        hold off
        i
        sigma = 10;  % Standard deviation for Gaussian kernel
        smoothedCounts = imgaussfilt(counts_tab_average{i}, sigma);
        % smoothedCounts_vect = smoothedCounts(:);
        % smoothedCount_vect_sorted = sort(smoothedCounts_vect);
        % smoothedCount_vect_sorted_cum_sum = cumsum(smoothedCount_vect_sorted);
        % mask1 = smoothedCount_vect_sorted_cum_sum > 0.25*smoothedCount_vect_sorted_cum_sum(end);
        % threshold = min(smoothedCount_vect_sorted(mask1));
        % [C1, h1] = contour(XGrid, YGrid, smoothedCounts, [-1 threshold ], 'LineColor', 'b', 'Showtext', 'on'); hold on;
        % mask2 = smoothedCount_vect_sorted_cum_sum > 0.5*smoothedCount_vect_sorted_cum_sum(end);
        % threshold = min(smoothedCount_vect_sorted(mask2));
        % [C2, h2] = contour(XGrid, YGrid, smoothedCounts, [-1 threshold], 'LineColor', 'g', 'Showtext', 'on');  % 'k' for black lines
        % mask3 = smoothedCount_vect_sorted_cum_sum > 0.75*smoothedCount_vect_sorted_cum_sum(end);
        % threshold = min(smoothedCount_vect_sorted(mask3));
        % [C3, h3] = contour(XGrid, YGrid, smoothedCounts, [-1 threshold], 'LineColor', 'r', 'Showtext', 'on');
        %pcolor(XGrid,YGrid,smoothedCounts); shading flat;
        %title(datestr(time(i)));
        %drawnow()
        save(['30days_heatmaps_gaussian/heatmap_average_gaussian' num2str(i) '.mat'], 'smoothedCounts');
        %caxis([0 50])
    end
    toc
end
%%
for i = istart:iend
    tic
    i
    load(['30days_heatmaps_gaussian/heatmap_average_gaussian' num2str(i) '.mat']);
    counts_tab_smoothed{i} = smoothedCounts;
    toc
end
%%
for i = istart:iend
    tic
    i
    if cod_case <=4
        smoothedCounts_vect = counts_tab_smoothed{i}(:);
        smoothedCount_vect_sorted = sort(smoothedCounts_vect);
        smoothedCount_vect_sorted_cum_sum = cumsum(smoothedCount_vect_sorted);
        mask1 = smoothedCount_vect_sorted_cum_sum > 0.25*smoothedCount_vect_sorted_cum_sum(end);
        threshold = min(smoothedCount_vect_sorted(mask1));
        [C1, h1] = contour(XGrid, YGrid, smoothedCounts, [-1000 threshold ], 'LineColor', 'b', 'Showtext', 'on'); hold on;
        mask2 = smoothedCount_vect_sorted_cum_sum > 0.5*smoothedCount_vect_sorted_cum_sum(end);
        threshold = min(smoothedCount_vect_sorted(mask2));
        [C2, h2] = contour(XGrid, YGrid, smoothedCounts, [-1000 threshold], 'LineColor', 'g', 'Showtext', 'on');  % 'k' for black lines
        mask3 = smoothedCount_vect_sorted_cum_sum > 0.75*smoothedCount_vect_sorted_cum_sum(end);
        threshold = min(smoothedCount_vect_sorted(mask3));
        [C3, h3] = contour(XGrid, YGrid, smoothedCounts, [-1000 threshold], 'LineColor', 'r', 'Showtext', 'on');
        for icase = 1:3
            if icase == 1
                C=C1;
            elseif icase == 2
                C=C2;
            elseif icase == 3
                C=C3;
            end
            contourLevels = C(1, :);
            polygonIndex = 1;
            ip = 1;
            while ip < size(C, 2)
                level = C(1, ip);
                numVertices = C(2, ip);
                xVertices = C(1, ip+1:ip+numVertices);
                yVertices = C(2, ip+1:ip+numVertices);
                contourData(polygonIndex).level = level;
                contourData(polygonIndex).polygon = [xVertices', yVertices'];
                polygonIndex = polygonIndex + 1;
                ip = ip + numVertices + 1;
                % Save the polygons into a .mat file
            end
            save(['30days_contours_smoothed/contour_polygons_index' num2str(icase) '_' num2str(i) '.mat'], 'contourData');
        end
    end
    toc
end

%
%%
figure
for i = istart : istart + 10
    hold off;
    for imod = 1:2
        if imod == 1
            str_mod = 'hycom';
        elseif imod == 2
            str_mod = 'glorys';
        end

        for icase = 1:3
            load([str_mod '/30days_contours_smoothed/contour_polygons_index' num2str(icase) '_' num2str(i) '.mat']);
            len_pol = [];
            count_pol = [];
            for j=1:length(contourData)
                len_pol = [len_pol; length(contourData(j).polygon)];
                count_pol = [count_pol; j];
            end
            max_pol = count_pol(len_pol == max(len_pol));
            max_pol = max_pol(1);
            if icase == 1
                c1_tab{i} = contourData;
                max_pol1(i) = max_pol;
            elseif icase == 2
                c2_tab{i} = contourData;
                max_pol2(i) = max_pol;
            elseif icase == 3
                c3_tab{i} = contourData;
                max_pol3(i) = max_pol;
            end
        end

        if imod == 1
            plot(c1_tab{i}(max_pol1(i)).polygon(:,1),c1_tab{i}(max_pol1(i)).polygon(:,2),'b--');hold on;
            plot(c2_tab{i}(max_pol2(i)).polygon(:,1),c2_tab{i}(max_pol2(i)).polygon(:,2),'g--');
            plot(c3_tab{i}(max_pol3(i)).polygon(:,1),c3_tab{i}(max_pol3(i)).polygon(:,2),'r--');
        else
            plot(c1_tab{i}(max_pol1(i)).polygon(:,1),c1_tab{i}(max_pol1(i)).polygon(:,2),'b:');hold on;
            plot(c2_tab{i}(max_pol2(i)).polygon(:,1),c2_tab{i}(max_pol2(i)).polygon(:,2),'g:');
            plot(c3_tab{i}(max_pol3(i)).polygon(:,1),c3_tab{i}(max_pol3(i)).polygon(:,2),'r:');
        end
    end
    axis([-160,-125,20,45])
    grid on;
    drawnow
end
%%
%%
for i = istart : iend
    i
    for imod = 2
        if imod == 1
            str_mod = 'hycom';
        elseif imod == 2
            str_mod = 'glorys';
        end

        for icase = 1:3
            load([str_mod '/30days_contours_smoothed/contour_polygons_index' num2str(icase) '_' num2str(i) '.mat']);
            surf = 0;
            for j=1:length(contourData)
                if contourData(j).level > 0
                    lat = contourData(j).polygon(:,2);
                    lon = contourData(j).polygon(:,1);
                    surf = surf + extract_area(lat,lon);
                end
            end
            if icase == 1
                gpgp_area1(i) = surf;
            elseif icase == 2
                gpgp_area2(i) = surf;
            elseif icase == 3
                gpgp_area3(i) = surf;
            end
        end
    end
end
%%
figure
plot(time(istart:iend),gpgp_area3(istart:iend));