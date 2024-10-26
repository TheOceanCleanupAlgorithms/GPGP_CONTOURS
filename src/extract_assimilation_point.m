icurrent = 0;
dx = 0.5;
gpgp_boundaries = [20 45 -160 -125];
gpgp_boundaries_p = [20-dx/2 45+dx/2 -160-dx/2 -125+dx/2];
iyear = 0;
icurrentyear = 0;
binEdgesX = linspace(gpgp_boundaries(3),gpgp_boundaries(4),(gpgp_boundaries(4)-gpgp_boundaries(3))/dx+1);
binEdgesY = linspace(gpgp_boundaries(1),gpgp_boundaries(2),(gpgp_boundaries(2)-gpgp_boundaries(1))/dx+1);
binEdgesX_p = linspace(gpgp_boundaries_p(3),gpgp_boundaries_p(4),(gpgp_boundaries_p(4)-gpgp_boundaries_p(3))/dx+1);
binEdgesY_p = linspace(gpgp_boundaries_p(1),gpgp_boundaries_p(2),(gpgp_boundaries_p(2)-gpgp_boundaries_p(1))/dx+1);
[YGrid, XGrid] = meshgrid(binEdgesY(1:end),binEdgesX(1:end));
%%
for i = 1:6051
    i
    load(['raw_heatmaps2/heatmap_' num2str(i) '.mat']);
    counts_tab{i} = counts;
end
%%
headers={'lat_id','lon_id','value','variance','time'};
for i = 1:120:365
    i
    lat_id = ones(1,size(XGrid,1))*(50-floor(size(XGrid,2)/2)); 	
    lon_id = 0:size(XGrid,1)-1;
    values = (counts_tab{i}(lon_id+1,floor(size(XGrid,2)/2)+1))';	
    variance = ones(1,size(XGrid,1))*0.1;
    time_mat = ones(1,size(XGrid,1))*(time(i)-time(1));
    if i == 1
        t = table(lat_id', lon_id', values', variance', time_mat', 'VariableNames', headers);
    else
        t1 = table(lat_id', lon_id', values', variance', time_mat', 'VariableNames', headers);
        t = [t;t1];
    end
end
writetable(t,'observations.csv')

