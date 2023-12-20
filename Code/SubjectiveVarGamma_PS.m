%Read Transfer task data
%starting with Lee's data
%internal variability: SD of their errors from the centroid of the dots
%%house keeping
%%use partial, everything the same except load and save _P 
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')

addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp1_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp1_validation\'

cd(dir);
%% setting up some basics
phase = 0; %0
pickGood = 1;
if phase == 0 && pickGood == 1
    T = readtable('SP Likelihood-Only Task.csv');
    load('Learning4Stats_GoodG_P.mat'); %serial then parallel
    load('SubjData.mat','nGdSubj','PSubj','SSubj','nSS','nPS','nType')
end
%% data table
subj = T.Participant;
type = unique(T.LLType);
nLType = size(type,1);

rmsq = zeros(nGdSubj,nLType+1); 
sLSd = zeros(nGdSubj,nLType+1); %subjective likelihood Sd
sLVar = zeros(nGdSubj,nLType+1); %subjective likelihood Var
sLGamma = zeros(nGdSubj,nLType+1); %%col1 narrow .060 col2 wide .15

sPrSd_p = zeros(nGdSubj,nLType+1);%with transformation 
sPrVar_p = zeros(nGdSubj,nLType+1);% p = prime

sPrGamma_p = zeros(nGdSubj,nLType+1);
sPrVtmp = zeros(nGdSubj,nLType+1);

sw_o = zeros(nGdSubj,nLType+1);
sw_op = zeros(nGdSubj,nLType+1);

sPrSd_r = zeros(nGdSubj,nLType+1); %raw without transformation
sPrVar_r = zeros(nGdSubj,nLType+1);
sPrGamma_r = zeros(nGdSubj,nLType+1);
%% get subjective likelihood sd, variance and gamma
nSj = nGdSubj;
for i = 1:nSj %row: subj
    if  i <= nSS
        tmpS = SSubj(i);
    else
        tmpS = PSubj(i-nSS);
    end
    rmsq(i,3) = tmpS;
    sLSd(i,3) = tmpS;
    sLVar(i,3) = tmpS;
    sLGamma(i,3) = tmpS;
    for k = 1:nLType %col1: narrow, col2: wide
        if i <= nSS
            %select trials
            if k == 1 %LL narrow
                useTmp = intersect(find(subj == SSubj(i)),find(contains(T.LLType,'Narrow')));
            elseif k == 2
                useTmp = intersect(find(subj == SSubj(i)),find(contains(T.LLType,'Wide')));
            end
        elseif i > nSS
            if k == 1 %narrow
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(T.LLType,'Narrow')));
            elseif k == 2
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(T.LLType,'Wide')));
            end
        end
        useTr = intersect(useTmp,find(contains(T.Outlier,'no')));

        splashX = table2array(T(useTr,8)); %8 LLMean
        netX = table2array(T(useTr,15)); %15
        delta = netX-splashX;
        sLSd(i,k) = sqrt(sum(delta.^2)/size(useTr,1));  
        sLVar(i,k) = sLSd(i,k)^2;
        sLGamma(i,k) = 1/sLVar(i,k);
    end
end
%%
sfacec = [0.3 0.2 0.11];

tmpSL = sLSd(:,1:2);
drawSL = tmpSL(:) ;
gele = 1:2;
xg = zeros(size(sLSd,1),nLType);
xg = repmat(gele,size(sLSd,1),1);
xgc = xg(:);

idx1 = isoutlier(sLSd(:,1),'quartiles');
outliers1 = sLSd(idx1,1);
idx2 = isoutlier(sLSd(:,2),'quartiles');
outliers2 = sLSd(idx2,2);
out = [idx1;idx2];

figure; hold on; 
bL = boxchart(xgc,drawSL);
bL.JitterOutliers = 'on';

ls = scatter(xgc(out == 0),drawSL(out == 0),[],sfacec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.1) 

for mm = 1:size(xg,1)
    plot(xg(mm,1:2),tmpSL(mm,1:2),"Color",[.55 .55 .55],"LineStyle","-")
end

axis square;
xtl = [1 2];
xtll = {'Ln','Lw'};
yt = [.01 .02 .03 .04 .05];
xlabel('Likelihood conditions')
set(gca(),'xtick', xtl, 'XTickLabel', xtll);% 'TickLabelInterpreter', 'latex');
set(gca(),'ytick', yt);
set(gca(),'FontSize', 16);
ylabel('Root Mean Sqaure Errors');
grid on; 

ylim([0.003 0.057])

% t1 = sLSd(:,1);
% t2 = sLSd(:,2);
% mean(sLSd)
% ttest(t1,t2,'Alpha',0.05)
%% 2. subjective prior sd and variance and gamma based on subjective likelihood and weight
for i = 1:nSj
    if  i <= nSS
        tmpS = SSubj(i);
    else
        tmpS = PSubj(i-nSS);
    end
    sPrSd_p(i,3) = tmpS;
    sPrVar_p(i,3) = tmpS;
    sPrGamma_p(i,3) = tmpS;
    sPrVtmp(i,3) = tmpS;
    sPrSd_r(i,3) = tmpS;
    sPrVar_r(i,3) = tmpS;
    sPrGamma_r(i,3) = tmpS;
    sw_o(i,3) = tmpS;
    sw_op(i,3) = tmpS;

    for k = 1:nLType %1: narrow, 2: wide
        if i <= nSS
            row = intersect(find(plttable.participant == SSubj(i)),find(plttable.LLabel == k));
        elseif i > nSS
            row = intersect(find(plttable.participant == PSubj(i-nSS)),find(plttable.LLabel == k));
        end

        sw = plttable.beta(row);
        sw_p = 1/(1+exp(-sw)); %sw'%transformation 
        if plttable.PLabel(row) == 1
            PrSD = .025; %pair Pr is narrow place
            p = 1;%pair Pr is narrow 
        else
            PrSD = .085;
            p = 2;%pair Pr is wide
        end

        sPrVtmp(i,p) = ((sLSd(i,k)^2/5)*sw_p)/(1-sw_p); %part1
        %transform back using transformed optimal weight
        sw_o(i,p) = PrSD^2/(PrSD^2+sLSd(i,k)^2/5); %objective prior
        sw_op(i,p) = 1/(1+exp(-sw_o(i,k)));
        sPrVar_p(i,p) = sPrVtmp(i,p)*sw_op(i,p)/(1-sw_op(i,p));
        sPrSd_p(i,p) = sqrt(sPrVar_p(i,p));
        sPrGamma_p(i,p) = 1/sPrVar_p(i,p);

        sPrVar_r(i,p) = ((sLSd(i,k)^2/5)*sw)/(1-sw); %subjective prior noise raw
        if sPrVar_r(i,p)>0
            sPrSd_r(i,p) = sqrt(sPrVar_r(i,p));
            sPrGamma_r(i,p) = 1/sPrVar_r(i,p);
        else
            sPrSd_r(i,p) = nan;
            sPrGamma_r(i,p) = nan;
        end
    end
end
    %%
figure(2); 
% subplot(1,2,1)
% hold on; 
% scatter(sPrVtmp(:),sPrVar_p(:),'filled');
% xlabel('logistically transformed prior variance')
% xlim([0 0.02])
% ylim([0 0.02])
% subplot(1,2,1)
% hold on;
% scatter(sPrVar_r(:),sPrVar_p(:),'filled');
% ylabel('logistically transformed prior variance')
% xlabel('un-transformed prior variance')
% xlim([0 0.02])
% ylim([0 0.02])
% axis square 
subplot(1,2,1)
hold on; 
xL = [-14 -4];
plot(xL,xL,'k');
s1 = scatter(log(sPrVar_r(1:nSS,1)),log(sPrVar_p(1:nSS,1)),...
    120,[.2 .2 .2],'filled');
s2 = scatter(log(sPrVar_r(1:nSS,2)),log(sPrVar_p(1:nSS,2)),...
    120,[.7 .7 .7],'filled');

% scatter(sPrGamma_r(nExS+1:end,1),sPrGamma_r(nExS+1:end,3),'filled');
xlabel('Subj prior variance-untransformed (log)')
ylabel('Subj prior variance-transformed(log)')
legend([s1 s2],{' narrow prior',' wide prior'})%'Location','northeastoutside'
title('serial-learning')
xlim([-14 -4])
ylim([-14 -4])
axis square 
set(gca,'LineWidth',1)
set(gca,'FontSize',16)


subplot(1,2,2)
hold on; 
plot(xL,xL,'k');
% scatter(sPrGamma_r(1:nExS,1),sPrGamma_r(1:nExS,2),'filled');
s3 = scatter(log(sPrVar_r(nSS+1:end,1)),log(sPrVar_p(nSS+1:end,1)),120,[.2 .2 .2],'filled');
s4= scatter(log(sPrVar_r(nSS+1:end,2)),log(sPrVar_p(nSS+1:end,2)),120,[.7 .7 .7],'filled');

xlabel('Subj prior variance-untransformed (log)')
ylabel('Subj prior variance-transformed(log)')
legend([s3 s4],{' narrow prior',' wide prior'})%'Location','northeastoutside'

% xlabel('Prior precision estimated by Likelihood(narrow)')
% ylabel('Subjective Prior precision (Lwide)')
title('parallel-learning')
xlim([-14 -4])
ylim([-14 -4])
axis square 
set(gca,'LineWidth',1)
set(gca,'FontSize',16)
%%
% scatter(sPrGamma_p(1:nExS,1),sPrGamma_p(1:nExS,2),'filled');
%xlabel('logistically transformed prior variance')
% xlim([0 0.02])
% ylim([0 0.02])
sPSdpn1 = sPrSd_p(1:nSS,1);%s:subjective P:prior%V: variance 
sPSdpw1 = sPrSd_p(1:nSS,2);%p: prime %serial group
sPSdpn2 = sPrSd_p(nSS+1:end,1);%s:subjective P:prior%V: variance 
sPSdpw2 = sPrSd_p(nSS+1:end,2);%p: prime %parallel group


itmp11 = find(sPrSd_r(1:nSS,1)>0);
itmp12 = find(sPrSd_r(1:nSS,2)>0);
itmp1 = intersect(itmp11,itmp12); %serial paritcipant to use

itmp21 = find(sPrSd_r(nSS+1:end,1)>0)+nSS;
itmp22 = find(sPrSd_r(nSS+1:end,1)>0)+nSS;
itmp2 = intersect(itmp21,itmp22); %parallel use

swtest(sPrSd_r(:))
sPSdrn1 = sPrSd_r(itmp1,1);%serial
sPSdrw1 = sPrSd_r(itmp1,2);
sPSdrn2 = sPrSd_r(itmp2,1);%2: parallel 
sPSdrw2 = sPrSd_r(itmp2,2);

% [h2r,p2r] = ttest(sPSdpn1,sPSdpw1)
ms_sPSdPn1 = nanmean(sPSdpn1)%subplot2(corrected), right(extra)
ms_sPSdPw1 = nanmean(sPSdpw1)
ers_sPSdPn1 = sqrt(std(sPSdpn1).^2/nSS);
ers_sPSdPw1 = sqrt(std(sPSdpw1).^2/nSS);

% % [h1r,p1r] = ttest(sPSdrn1,sPSdrm)
ms_sPSdRn1 = nanmean(sPSdrn1)%subplot1 right
ms_sPSdRw1 = nanmean(sPSdrw1)
ers_sPSdRn1 = sqrt(std(sPSdrn1).^2/nSS);%er3
ers_sPSdRw1 = sqrt(std(sPSdrw1).^2/nSS);%er4

% [h2l,p2l] = ttest(sPSdpn2,sPSdpw)
mp_sPSdPn2 = nanmean(sPSdpn2)%subplot2 left (inter)
mp_sPSdPw2 = nanmean(sPSdpw2)
erp_sPSdPn2 = sqrt(std(sPSdpn2).^2/nPS);
erp_sPSdPw2 = sqrt(std(sPSdpw2).^2/nPS);

% [h1l,p1l] = ttest(sPSdrn2,sPSdrw)
mp_sPSdRn2 = nanmean(sPSdrn2)%subplot1(uncorrected) left(inter)
mp_sPSdRw2 = nanmean(sPSdrw2)
erp_sPSdRn2 = sqrt(std(sPSdrn2).^2/nPS);%er1
erp_sPSdRw2 = sqrt(std(sPSdrw2).^2/nPS);%er2

xt = [1-0.18,1+0.18,2-0.18,2+0.18]; %xtick
xtg = [1 2]; %xtick

figure(3); 
subplot(1,2,1);
hold on;
b1 = bar(xt(1),ms_sPSdRn1,0.3,'FaceColor',[.2 .2 .2])%serial uncorrected
er1 = errorbar(xt(1),ms_sPSdRn1,ers_sPSdRn1,ers_sPSdRn1);%narrow
er1.Color = [0 0 0];                            
er1.LineStyle = 'none'; 
er1.LineWidth = 1; 

b2 = bar(xt(2),ms_sPSdRw1,0.3,'FaceColor',[.75 .75 .75])%serial 
er2 = errorbar(xt(2),ms_sPSdRw1,ers_sPSdRw1,ers_sPSdRw1);%wide
er2.Color = [0 0 0];                            
er2.LineStyle = 'none'; 
er2.LineWidth = 1; 

b3 = bar(xt(3),mp_sPSdRn2,0.3,'FaceColor',[.2 .2 .2])%parallel %uncorrrected 
er3 = errorbar(xt(3),mp_sPSdRn2,erp_sPSdRn2,erp_sPSdRn2);%narrow 
er3.Color = [0 0 0];                            
er3.LineStyle = 'none'; 
er3.LineWidth = 1; 

b4 = bar(xt(4),mp_sPSdRw2,0.3,'FaceColor',[.75 .75 .75])%parallel 
er4 = errorbar(xt(4),mp_sPSdRw2,erp_sPSdRw2,erp_sPSdRw2);%wide
er4.Color = [0 0 0];                            
er4.LineStyle = 'none'; 
er4.LineWidth = 1; 

xline(1.5,'-')
yline(.025,'--','LineWidth',.8)
yline(.085,'--','LineWidth',.8)

ylim([0 0.09])
axis square
title('Un-transformed')
ylabel('estimated  \sigma_{Pr}')
x1 = {'serial','parallel'};
set(gca(),'xtick', xtg, 'XTickLabel', x1);% 'TickLabelInterpreter', 'latex');
set(gca,'LineWidth',1)
set(gca,'FontSize',16)


subplot(1,2,2);
hold on;
b5 = bar(xt(1),ms_sPSdPn1,0.3,'FaceColor',[.2 .2 .2])
er5 = errorbar(xt(1),ms_sPSdPn1,ers_sPSdPn1,ers_sPSdPn1);
er5.Color = [0 0 0];                            
er5.LineStyle = 'none';er5.LineWidth = 1; 
 
b6 = bar(xt(2),ms_sPSdPw1,0.3,'FaceColor',[.75 .75 .75])
er6 = errorbar(xt(2),ms_sPSdPw1,ers_sPSdPw1,ers_sPSdPw1);
er6.Color = [0 0 0];                            
er6.LineStyle = 'none';er6.LineWidth = 1; 

b7 = bar(xt(3),mp_sPSdPn2,0.3,'FaceColor',[.2 .2 .2])
er7 = errorbar(xt(3),mp_sPSdPn2,erp_sPSdPn2,erp_sPSdPn2);
er7.Color = [0 0 0];                            
er7.LineStyle = 'none';er7.LineWidth = 1; 

b8 = bar(xt(4),mp_sPSdPw2,0.3,'FaceColor',[.75 .75 .75])%parallel 
er8 = errorbar(xt(4),mp_sPSdPw2,erp_sPSdPw2,erp_sPSdPw2);%wide
er8.Color = [0 0 0];                            
er8.LineStyle = 'none';er8.LineWidth = 1; 

xline(1.5,'-')
yline(.025,'--','LineWidth',.8)
yline(.085,'--','LineWidth',.8)

axis square
title('Transformed')
ylim([0 0.09])
% xlim([0.5 2.5])
axis square
x2 = x1;
% text(1-.06, 0.12,'***','FontSize',12);text(2-.06, 0.06,'***','FontSize',12);
% line([xt(1),xt(2)],[.118,.118],'color','k')
% line([xt(3),xt(4)],[.058,.058],'color','k')
legend([b1 b2],{'narrow prior','wide prior'})%'Location','northeastoutside'

set(gca(),'xtick', xtg, 'XTickLabel', x2);% 'TickLabelInterpreter', 'latex');
set(gca,'LineWidth',1)
set(gca,'FontSize',18)
%% get mean precision (average) 
msPrGamma_p = sum(sPrGamma_p,2);%/2
sPrGamma_r(sPrGamma_r == 0) = NaN;
msPrGamma_r = nanmean(sPrGamma_r,2);
figure;scatter(sPrGamma_r(:),sPrGamma_p(:));%
%displayed that there is an inflation of prior precision (around 10 fold)
%using the adjusted one (Gamma_prime) 

% save('subjective_variances1_P.mat','sLSd','sLVar',...
%     'sLGamma','sPrSd_p','sPrVar_p','sPrGamma_p',...
%     'msPrGamma_p','sPrVtmp','sw_o','sw_op','sPrSd_r',...
% 'sPrVar_r','sPrGamma_r','msPrGamma_r');