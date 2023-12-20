%Read Transfer task data
%starting with Lee's data 
%internal variability: SD of their errors from the centroid of the dots 
%%house keeping
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_validation\'

cd(dir);
%% setting up some basics 
phase = 0; %0 
pickGood = 1; 
if phase == 0 && pickGood == 1
T = readtable('IE Likelihood-Only Task.csv');
load('Learning4Stats_Good_P.mat','exSubj','inSubj','nExS',...
    'nInS','badSubj','nSubj','nGdSubj','nType',...
    'nCond','DataTable','LLSd','Wop','PrSD');
end 
%% data table 
subj = T.Participant;
type = T.LLType;

rmsq = zeros(nGdSubj,nType); 
sLSd = zeros(nGdSubj,nType); %subjective likelihood Sd
sLVar = zeros(nGdSubj,nType); %subjective likelihood Var
sLGamma = zeros(nGdSubj,nType); %

sPrSd_p = zeros(nGdSubj,nType);
sPrVar_p = zeros(nGdSubj,nType);
sPrGamma_p = zeros(nGdSubj,nType);
msPrGamma_p = zeros(nGdSubj,1);
sPrVtmp = zeros(nGdSubj,nType);
sw_o = zeros(nGdSubj,nType);
sw_op = zeros(nGdSubj,nType);

sPrSd_r = zeros(nGdSubj,nType); %raw without transformation
sPrVar_r = zeros(nGdSubj,nType);
sPrGamma_r = zeros(nGdSubj,nType);
msPrGamma_r = zeros(nGdSubj,1);
%% get subjective likelihood sd and variance 
nSj = nGdSubj; 
for i = 1:nSj
    for k = 1:nType %1: narrow, 2: medium, 3: wide
        if i <= nExS
            %select trials
            if k == 1 %narrow              
                useTmp = intersect(find(subj == exSubj(i)),find(contains(T.LLType,'Narrow')));
            elseif k == 2
                useTmp = intersect(find(subj == exSubj(i)),find(contains(T.LLType,'Medium')));
            elseif k == 3
                useTmp = intersect(find(subj == exSubj(i)),find(contains(T.LLType,'Wide')));
            end 
        elseif i > nExS
            if k == 1 %narrow
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(T.LLType,'Narrow')));
            elseif k == 2
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(T.LLType,'Medium')));
            elseif k == 3
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(T.LLType,'Wide')));
            end
        end
        useTr = intersect(useTmp,find(contains(T.Outlier,'no')));

        splashX = table2array(T(useTr,7)); %7 LLMean %8 LSD
        netX = table2array(T(useTr,14));
        delta = netX-splashX;

%         rmsq(i,k) = rms(delta);
sLSd(i,k) = sqrt(sum(delta.^2)/size(useTr,1)); %rms
%         sLSd(i,k) = std(delta);
%         sLSd(i,k) = rmsq(i,k)
        sLVar(i,k) = sLSd(i,k)^2;
        sLGamma(i,k) = 1/sLVar(i,k);
    end
end
%% draw subjective likelihood 
sfacec = [0.3 0.2 0.11];
drawSL = sLSd(:);
gele = 1:3;
xg = zeros(size(sLSd,1),size(sLSd,2));
xg = repmat(gele,size(sLSd,1),1);
xgc = xg(:);
figure(1); hold on;
bL = boxchart(xgc,drawSL);
bL.JitterOutliers = 'on';

idx1 = isoutlier(sLSd(:,1),'quartiles');
outliers1 = sLSd(idx1,1);
idx2 = isoutlier(sLSd(:,2),'quartiles');
outliers2 = sLSd(idx2,2);
idx3 = isoutlier(sLSd(:,3),'quartiles');
outliers3 = sLSd(idx2,3);
out = [idx1;idx2;idx3];

ls = scatter(xgc(out == 0),drawSL(out == 0),[],sfacec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.16) 
axis square;
xlim([0 4]);
ylim([0 .06])
xtl = [1 2 3];
xtll = {'Ln','Lm','Lw'};
yt = [.01 .02 .03 .04 .05];
xlabel('Likelihood conditions')
set(gca(),'xtick', xtl, 'XTickLabel', xtll);% 'TickLabelInterpreter', 'latex');
set(gca(),'ytick', yt);
set(gca(),'FontSize', 16);

ylabel('Root Mean Sqaure Errors');
grid on;
%%
tmpg = (reshape(DataTable.groups,3,size(DataTable,1)/3))';
t = table(tmpg(:,1),sLSd(:,1),sLSd(:,2),sLSd(:,3),...
'VariableNames',{'group','L1','L2','L3'});
LL = table([1 2 3]','VariableNames',{'Likelihood'});
rm = fitrm(t,'L1-L3~group','WithinDesign',LL);
ranovatbl = ranova(rm)

%% 2. subjective prior sd and variance and gamma  
% for k = 1:3
% Wop_p(k) = 1/(1+exp(-Wop(k)))
% end
for i = 1:nSj
    for k = 1:nType %1: narrow, 2: medium, 3: wide
        if i <= nExS
            row = intersect(find(DataTable.Subjects == exSubj(i)),find(DataTable.Likelihoods == LLSd(k)));
        elseif i > nExS
            row = intersect(find(DataTable.Subjects == inSubj(i-nExS)),find(DataTable.Likelihoods == LLSd(k)));
        end
        
        sw = DataTable.Weights(row);
        sw_p = 1/(1+exp(-sw)); %sw'
        if sw == 0 %learning phase 0 means no this likelihood condition
            continue
        else
            sPrVtmp(i,k) = ((sLSd(i,k)^2/5)*sw_p)/(1-sw_p); %part1
            %transform back using transformed optimal weight
            sw_o(i,k) = PrSD^2/(PrSD^2+sLSd(i,k)^2/5);
            sw_op(i,k) = 1/(1+exp(-sw_o(i,k)));
            sPrVar_p(i,k) = sPrVtmp(i,k)*sw_op(i,k)/(1-sw_op(i,k));
            sPrSd_p(i,k) = sqrt(sPrVar_p(i,k));
            sPrGamma_p(i,k) = 1/sPrVar_p(i,k);

            sPrVar_r(i,k) = ((sLSd(i,k)^2/5)*sw)/(1-sw); %part1
            if sPrVar_r(i,k)>0
            sPrSd_r(i,k) = sqrt(sPrVar_r(i,k));
            sPrGamma_r(i,k) = 1/sPrVar_r(i,k);
            else
            sPrSd_r(i,k) = nan;  
            sPrGamma_r(i,k) = nan;
            end
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
% subplot(2,2,4)
% hold on; 
% scatter(sPrGamma_r(1:nExS,1),sPrGamma_r(1:nExS,2),'filled');
% scatter(sPrGamma_p(1:nExS,1),sPrGamma_p(1:nExS,2),'filled');

% scatter(sPrGamma_r(nExS+1:end,1),sPrGamma_r(nExS+1:end,3),'filled');
xlabel('Subjective Prior precision (Lnarrow)')
ylabel('Subjective Prior precision (Lmedium)')
title('extrapolation-learning')
axis square 

subplot(2,2,2)
hold on; 
% scatter(sPrGamma_r(1:nExS,1),sPrGamma_r(1:nExS,2),'filled');
scatter(sPrGamma_r(nExS+1:end,1),sPrGamma_r(nExS+1:end,3),'filled');
scatter(sPrGamma_p(nExS+1:end,1),sPrGamma_p(nExS+1:end,3),'filled');

% xlabel('Prior precision estimated by Likelihood(narrow)')
ylabel('Subjective Prior precision (Lwide)')
title('Interpolation-learning')
axis square 

% scatter(sPrGamma_p(1:nExS,1),sPrGamma_p(1:nExS,2),'filled');
%xlabel('logistically transformed prior variance')
% xlim([0 0.02])
% ylim([0 0.02])
sPSdpn1 = sPrSd_p(1:nExS,1);%s:subjective P:prior%V: variance 
sPSdpm = sPrSd_p(1:nExS,2);%p: prime 
sPSdpn2 = sPrSd_p(nExS+1:end,1);%s:subjective P:prior%V: variance 
sPSdpw = sPrSd_p(nExS+1:end,3);%p: prime 
itmp1 = find(sPrSd_r(1:nExS,1)>0);
itmp2 = find(sPrSd_r(nExS+1:end,1)>0)+nExS;

sPSdrn1 = sPrSd_r(itmp1,1);
sPSdrm = sPrSd_r(itmp1,2);
sPSdrn2 = sPrSd_r(itmp2,1);
sPSdrw = sPrSd_r(itmp2,3);

%[h2r,p2r] = ttest(sPSdpn1,sPSdpm)
mex_sPSdPn1 = mean(sPSdpn1,"omitnan")%subplot2(corrected), right(extra)
mex_sPSdPm = mean(sPSdpm,"omitnan")
erex_sPSdPn1 = sqrt(nanstd(sPSdpn1).^2/nExS);
erex_sPSdPm = sqrt(nanstd(sPSdpm).^2/nExS);

% [h1r,p1r] = ttest(sPSdrn1,sPSdrm)
mex_sPSdRn1 = mean(sPSdrn1,"omitnan")%subplot1 right
mex_sPSdRm = mean(sPSdrm,"omitnan")
erex_sPSdRn1 = sqrt(nanstd(sPSdrn1).^2/nExS);%er3
erex_sPSdRm = sqrt(nanstd(sPSdrm).^2/nExS);%er4

% [h2l,p2l] = ttest(sPSdpn2,sPSdpw)
min_sPSdPn2 = mean(sPSdpn2,"omitnan")%subplot2 left (inter)
min_sPSdPw = mean(sPSdpw,"omitnan")
erin_sPSdPn2 = sqrt(nanstd(sPSdpn2).^2/nInS);
erin_sPSdPw = sqrt(nanstd(sPSdpw).^2/nInS);

% [h1l,p1l] = ttest(sPSdrn2,sPSdrw)
min_sPSdRn2 = mean(sPSdrn2,"omitnan")%subplot1(uncorrected) left(inter)
min_sPSdRw = mean(sPSdrw,"omitnan")
erin_sPSdRn2 = sqrt(nanstd(sPSdrn2).^2/nInS);%er1
erin_sPSdRw = sqrt(nanstd(sPSdrw).^2/nInS);%er2

xt = [1-0.18,1+0.18,2-0.18,2+0.18]; %xtick
xtg = [1 2]; %xtick

figure; 
subplot(1,2,1);
hold on;
b1 = bar(xt(1),min_sPSdRn2,0.3,'FaceColor',[.15 .15 .15])
er1 = errorbar(xt(1),min_sPSdRn2,erin_sPSdRn2,erin_sPSdRn2);
er1.Color = [0 0 0];                            
er1.LineStyle = 'none'; 
er1.LineWidth = 1; 

b2 = bar(xt(2),min_sPSdRw,0.3,'FaceColor',[.85 .85 .85])
er2 = errorbar(xt(2),min_sPSdRw,erin_sPSdRw,erin_sPSdRw);
er2.Color = [0 0 0];                            
er2.LineStyle = 'none'; 
er2.LineWidth = 1; 

b3 = bar(xt(3),mex_sPSdRn1,0.3,'FaceColor',[.15 .15 .15])
er3 = errorbar(xt(3),mex_sPSdRn1,erex_sPSdRn1,erex_sPSdRn1);
er3.Color = [0 0 0];                            
er3.LineStyle = 'none'; 
er3.LineWidth = 1; 

b4 = bar(xt(4),mex_sPSdRm,0.3,'FaceColor',[.5 .5 .5])
er4 = errorbar(xt(4),mex_sPSdRm,erex_sPSdRm,erex_sPSdRm);
er4.Color = [0 0 0];                            
er4.LineStyle = 'none'; 
er4.LineWidth = 1; 
xline(1.5,'-')
yline(.025,'--','LineWidth',.8)
legend([b1 b4 b2],{'Lnarrow','Lmedian','Lwide'})%'Location','northeastoutside'



ylim([0 0.12])
axis square
title('Un-transformed')
ylabel('estimated  \sigma_{Pr}')
x1 = {'Interpolation','Extrapolation'};

set(gca(),'xtick', xtg, 'XTickLabel', x1);% 'TickLabelInterpreter', 'latex');
set(gca,'LineWidth',1)
set(gca,'FontSize',16)


subplot(1,2,2);
hold on;
b5 = bar(xt(1),min_sPSdPn2,0.3,'FaceColor',[.2 .2 .2])
er5 = errorbar(xt(1),min_sPSdPn2,erin_sPSdPn2,erin_sPSdPn2);
er5.Color = [0 0 0];                            
er5.LineStyle = 'none';er5.LineWidth = 1; 
 
b6 = bar(xt(2),min_sPSdPw,0.3,'FaceColor',[.5 .5 .5])
er6 = errorbar(xt(2),min_sPSdPw,erin_sPSdPw,erin_sPSdPw);
er6.Color = [0 0 0];                            
er6.LineStyle = 'none';er6.LineWidth = 1; 

b7 = bar(xt(3),mex_sPSdPn1,0.3,'FaceColor',[.2 .2 .2])
er7 = errorbar(xt(3),mex_sPSdPn1,erex_sPSdPn1,erex_sPSdPn1);
er7.Color = [0 0 0];                            
er7.LineStyle = 'none';er7.LineWidth = 1; 

b8 = bar(xt(4),mex_sPSdPm,0.3,'FaceColor',[.75 .75 .75])
er8 = errorbar(xt(4),mex_sPSdPm,erex_sPSdPm,erex_sPSdPm);
er8.Color = [0 0 0];                            
er8.LineStyle = 'none';er8.LineWidth = 1; 

xline(1.5,'-')
yline(.025,'--','LineWidth',.8)
axis square
title('Transformed')
ylim([0 0.12])
% xlim([0.5 2.5])
axis square
x2 = x1;
% text(1-.06, 0.12,'***','FontSize',12);text(2-.06, 0.06,'***','FontSize',12);
% line([xt(1),xt(2)],[.118,.118],'color','k')
% line([xt(3),xt(4)],[.058,.058],'color','k')
legend([b5 b8 b6],{'Ln','Lm','Lw'})%'Location','northeastoutside'

set(gca(),'xtick', xtg, 'XTickLabel', x2);% 'TickLabelInterpreter', 'latex');
set(gca,'LineWidth',1)
set(gca,'FontSize',16)
%% get mean precision (average) 
msPrGamma_p = mean(sPrGamma_p,2);%/2
sPrGamma_r(sPrGamma_r == 0) = NaN;
msPrGamma_r = nanmean(sPrGamma_r,2);
rmPrGamma = zeros(nSj,2);%reconstruct matrix 75X2
rmPrGamma(1:nExS,:) = sPrGamma_r(1:nExS,1:2);
rmPrGamma(nExS+1:end,1) = sPrGamma_r(nExS+1:end,1);
rmPrGamma(nExS+1:end,2) = sPrGamma_r(nExS+1:end,3);
v1 = [min(sPrGamma_r(:)),min(sPrGamma_r(:))];%line with slope = 1
v2 = [max(sPrGamma_r(:)),max(sPrGamma_r(:))];
pt = [rmPrGamma(:,1),rmPrGamma(:,2)];
distance_2D = point_to_line_distance(pt, v1, v2);
% pt(isnan(pt))=0;
figure;
hold on; scatter(pt(:,1),pt(:,2));
hL=plot([min(sPrGamma_r(:)) max(sPrGamma_r(:))],[min(sPrGamma_r(:)) max(sPrGamma_r(:))]); 
axis square
figure;scatter(msPrGamma_r(:),msPrGamma_p(:),'filled');%
axis square;
xlabel('untransformed mean prior precision')
ylabel('transformed mean prior precision')
set(gca,'FontSize',16)

%displayed that there is an inflation of prior precision (around 10 fold)
%using the adjusted one (Gamma_prime) 
%%
% save('subjective_variances1_P.mat','sLSd','sLVar',...
%     'sLGamma','sPrSd_p','sPrVar_p','sPrGamma_p',...
%     'msPrGamma_p','sPrVtmp','sw_o','sw_op','sPrSd_r',...
% 'sPrVar_r','sPrGamma_r','msPrGamma_r');