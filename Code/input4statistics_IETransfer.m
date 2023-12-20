%Read coin task data
%transfer phase 
%%house keeping
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_validation\'

cd(dir);
%% setting up some basics 
phase = 2; %1 learning 2 transfer
pickGood = 1; %0: all data, 1: pick based on R 2: use Mat isoutlier

T = readtable('IE Transfer Task.csv');
load('SubjDataIE.mat'); 
%%
nGrp = 2; %extrapolation %interpolation
subj = T.Participant;
type = T.LLType;
uSubj = unique(subj);
uType = unique(type);
nSubj = size(uSubj,1);
nType = size(uType,1);
nTrial = size(subj,1)/size(uSubj,1);
%% subject numbers, if pickGood == 1, remove bad participants
if phase == 1
    exSubj = unique(subj(find(contains(T.LLType,'Medium')))); %only ex experienced
    inSubj = unique(subj(find(contains(T.LLType,'Wide')))); %only in experienced
end
if pickGood == 1
badSubj = unique(subj(find(contains(T.Outlier,'bad'))));
exSubj = setdiff(exSubj,badSubj);
inSubj = setdiff(inSubj,badSubj);
nGdSubj = nSubj-size(badSubj,1);
end 
nExS = size(exSubj,1);
nInS = size(inSubj,1);
%% get weighting and LL
if phase == 1
    PrSD = T.PriorSD(1);
    LLSd = zeros(nType,1);
    Wop = zeros(nType,1);
    for ii = 1:nType
        if ii == 1
            LLSd(ii) = T.LikelihoodSD(find(contains(T.LLType,'Narrow'),1));
        elseif ii == 2
            LLSd(ii) = T.LikelihoodSD(find(contains(T.LLType,'Medium'),1));
        else
            LLSd(ii) = T.LikelihoodSD(find(contains(T.LLType,'Wide'),1));
        end
        Wop(ii) = PrSD^2/(PrSD^2+LLSd(ii)^2/5);
    end
else
end
%%
if pickGood == 1
    nSj = nGdSubj;
else
    nSj = nSubj;
end 
    plotData = zeros(nSj,nType*2+2);%1Wn 2Wm 3Ww 4 OIn 5 OIm 6 OIw 7grp 8subj 
    statData = zeros(nSj*nType,7);%1 sub 2 grp 3 LL 4 wts 5 oi 6 oic 7.NTr
for i = 1:nSj
    for k = 1:nType %1: narrow, 2: medium, 3: wide
        if i <= nExS
            plotData(i,7) = 1; %group/cond 1: ex 2: int
            plotData(i,8) = exSubj(i);
            statData(k+(i-1)*nType,1) =exSubj(i);
            statData(k+(i-1)*nType,2) =0;
            statData(k+(i-1)*nType,3) =LLSd(k);
            if k == 1 %narrow              
                useTmp = intersect(find(subj == exSubj(i)),find(contains(T.LLType,'Narrow')));
            elseif k == 2
                useTmp = intersect(find(subj == exSubj(i)),find(contains(T.LLType,'Medium')));
            else
                if phase == 1
                continue
                elseif phase == 2
                useTmp = intersect(find(subj == exSubj(i)),find(contains(T.LLType,'Wide'))); 
                end
            end 
        elseif i > nExS
            plotData(i,7) = 2;
            plotData(i,8) = inSubj(i-nExS);
            statData(k+(i-1)*nType,1) = inSubj(i-nExS);
            statData(k+(i-1)*nType,2) = 1;
            statData(k+(i-1)*nType,3) = LLSd(k);
            if k == 1 %narrow
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(T.LLType,'Narrow')));
            elseif k == 2
                if phase == 1
                continue
                elseif phase == 2
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(T.LLType,'Medium')));
                end
            else k == 3
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(T.LLType,'Wide')));
            end
        end
        if pickGood == 0
        useTr = useTmp;
        else
        useTr = intersect(useTmp,find(contains(T.Outlier,'no')));
        end
        statData(k+(i-1)*nType,8) = size(useTr,1);
        splashX = table2array(T(useTr,8)); %8 LLMean
        coinX = table2array(T(useTr,10));
        netX = table2array(T(useTr,16));
        fitdata = polyfit(splashX,netX,1);
        oitmp = table2array(T(useTr,19));
        oictmp = table2array(T(useTr,20));       
        Intercept(i,k) = fitdata(2);
        plotData(i,k) = fitdata(1);
        plotData(i,k+nType) = median(oitmp);
        statData(k+(i-1)*nType,4) = fitdata(1);
        statData(k+(i-1)*nType,5) = fitdata(2);
        statData(k+(i-1)*nType,6) = median(oitmp);
        statData(k+(i-1)*nType,7) = mean(oictmp);
    end
end
%% 2. summary of descriptive statistics
%statData %(40+35)*3 = 225
meanW = zeros(nGrp,nType); %2 X 3
meanOI = zeros(nGrp,nType); %2 X 3
sdW = zeros(nGrp,nType); %2 X 3
sdOI = zeros(nGrp,nType); %2 X 3
for aa = 1:nGrp
    if aa ==1
    meanW (aa,:) = mean(plotData(1:nExS,1:nType),1)
    meanOI (aa,:) = mean(plotData(1:nExS,nType+1:nType*2),1)
    sdW (aa,:) = std(plotData(1:nExS,1:nType),0,1)
    sdOI (aa,:) = std(plotData(1:nExS,nType+1:nType*2),0,1)
    elseif aa ==2
    meanW (aa,:) = mean(plotData(nExS+1:nSj,1:nType),1)
    meanOI (aa,:) = mean(plotData(nExS+1:nSj,nType+1:nType*2),1)
    sdW (aa,:) = std(plotData(nExS+1:nSj,1:nType),0,1)
    sdOI (aa,:) = std(plotData(nExS+1:nSj,nType+1:nType*2),0,1)
    end
end
%% 3.box plot 
sfacec = [0.9 0.59 0.33;0.6 0.4 0.22;0.3 0.2 0.11];

figure(1); %Weights
subplot(2,3,2);
hold on;
xge1 = repmat(1,nExS,1);
xge2 = repmat(2,nExS,1);
xge3 = repmat(3,nExS,1);


be1 = boxchart(xge1,plotData(1:nExS,1),'MarkerStyle','none') %col1: narrow col2 medium
be2 = boxchart(xge2,plotData(1:nExS,2),'MarkerStyle','none') %col1: narrow col2 medium
be3 = boxchart(xge3,plotData(1:nExS,3),'MarkerStyle','none') %col1: narrow

for jj =1:3
    xx = [jj-0.5 jj+0.5];
    yy = [Wop(jj) Wop(jj)]
    plot(xx,yy,'k','LineWidth',2);
    y1 = plotData(1:nExS,jj);
    x1=repmat(jj,size(y1,1),1);
    if jj== 1
        c=sfacec(1,:);
    elseif jj == 2
        c = sfacec(2,:);
    elseif jj == 3
        c = sfacec(3,:);
    end
    scatter(x1(:),y1,[],c,'filled','MarkerFaceAlpha',0.4',...
        'jitter','on','jitterAmount',0.15);
end
be1.BoxFaceColor = '#0072BD';
be2.BoxFaceColor = '#0072BD';
be3.BoxFaceColor = '#B3AF05';


ylim([0.4 1.1]);
xticks([1 2 3]);
xticklabels({'PwLn','PwLm','PwLw'});
ylabel('slope value');
if phase == 2
    title('Extra-Transfer');
elseif phase  == 1
    title('Extra-Learn');
end
axis square;
% box on;
grid on; 
ax =gca;
ax.FontSize = 16;

xgi1 = repmat(1,nInS,1);
xgi2 = repmat(2,nInS,1);
xgi3 = repmat(3,nInS,1);

subplot(2,3,1);
hold on; 
bi1 = boxchart(xgi1,plotData(nExS+1:nSj,1),'MarkerStyle','none') %col1: narrow 
bi2 = boxchart(xgi2,plotData(nExS+1:nSj,2),'MarkerStyle','none') %col2 medium
bi3 = boxchart(xgi3,plotData(nExS+1:nSj,3),'MarkerStyle','none') %col2 wide
for jj =1:3
    xx = [jj-0.5 jj+0.5];
    yy = [Wop(jj) Wop(jj)]
    plot(xx,yy,'k','LineWidth',2);    
    y1 = plotData(nExS+1:nSj,jj);
    x1=repmat(jj,size(y1,1),1);
    if jj== 1
        c=sfacec(1,:);
    elseif jj == 2
        c = sfacec(2,:);
    elseif jj == 3
        c = sfacec(3,:);
    end
    scatter(x1(:),y1,[],c,'filled','MarkerFaceAlpha',0.4',...
        'jitter','on','jitterAmount',0.15);
end
bi1.BoxFaceColor = '#0072BD';
bi2.BoxFaceColor = '#B3AF05';
bi3.BoxFaceColor = '#0072BD';


ylim([0.4 1.1]);
xticks([1 2 3]);
xticklabels({'Narrow','Medium','Wide'});
xticklabels({'PwLn','PwLm','PwLw'});
ylabel('slope value');
%
if phase == 2
    title('Inter-Transfer');
elseif phase  == 1
    title('Inter-Learn');
end

axis square;
grid on; 
ax =gca;
ax.FontSize = 16;
%% 4.
figure(2); %OIs
subplot(1,2,1);
hold on;
bo1 = boxchart(plotData(1:nExS,4:6))
yline(1);

ylim([0.4 1.1]);
ylabel('Optimality index');
title('Extrapolation');
box on;

subplot(1,2,2);
hold on; 
bo2 = boxchart(plotData(nExS+1:nSj,4:6))
yline(1);
bo2.BoxFaceColor = '#D95319';

ylim([0.4 1.1]);
ylabel('Optimality index');
title('Interpolation');
box on;
%% save data
clearvars ans coinX splashX dir fitdata i k T type ii jj xx yy 
DataTable = array2table(statData,...
    'VariableNames',{'Subjects','groups','Likelihoods','Weights','Intercept','OI','OIc','Tr'});
% if pickGood == 2
%     if phase == 1 %
%         LLuse = [LLSd(1) LLSd(2);LLSd(1) LLSd(3)]';%row1:e tow2:i
%         outlT = nan(nSj,4); %col1:outlier(0/1)LL1 col1:outlier(0/1)LL2 col3:Grp col4:subj
%         sbj = reshape(DataTable.Subjects,3,nSj);
%         outlT(:,4)  = sbj(1,:)';
%         gp = reshape(DataTable.groups,nSj,3);
%         outlT(:,3)  = gp(:,2);
%     elseif phase == 2
%         LLuse = repmat(LLSd,1,2);%row1:e tow2:i
%         outlT = nan(nSj,5); %col1:outlier(0/1)LL1 col1:outlier(0/1)LL2 col3:Grp col4:subj
%         sbj = reshape(DataTable.Subjects,3,nSj);
%         outlT(:,5)  = sbj(1,:)';
%         gp = reshape(DataTable.groups,nSj,3);
%         outlT(:,4)  = gp(:,2);
%     end
%     for gg =1:2 %1:Ex 2:In
%         for ll=1:size(LLuse,1) %LLuse
%             tmpi = intersect(find(DataTable.groups == gg-1),find(DataTable.Likelihoods == LLuse((gg-1)*size(LLuse,1)+ll)));
%             outlT(outlT(:,size(outlT,2)-1) == gg-1,ll) = isoutlier(DataTable.Weights(tmpi));
%         end
%     end
% end
% [ro,~]= find(outlT(:,1:3));
% sort(ro)
%% save data 
%         save('Transfer4Stats_Good.mat');