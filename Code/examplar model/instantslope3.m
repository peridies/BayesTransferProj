function slope = instantslope3(netX,splashX,Outlier)%instantenous slope3
%follwing 10
slope = zeros(size(netX,1),1);
for i = 1:size(netX,1)
    useTr = intersect((i:i+9),find(contains(Outlier,'no')));%take out outliers
    %         tmpNX = netX(useTr);
    %         tmpSX = splashX(useTr);
    if size(useTr,1) > 4
        fitdata = polyfit(splashX(useTr),netX(useTr),1);%linear fit, single individual
        slope(i) = fitdata(1);%instant Slope
    else
        slope(i) = nan;%instant Slope
    end
end