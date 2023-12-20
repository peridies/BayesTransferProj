%% thus we compute expected optimal weighting using Gamma_r(raw)
%Kiryakova 2020 equation 
%using subjective likelihood 
%%house keeping
%can choose whether use _P.mat or not 
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')

addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp1_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp1_validation\'

cd(dir);
%%
load('subjective_variances1_P.mat');
load('SubjData.mat','SSubj','PSubj','nSS',...
    'nPS');

load('OldNewTable_Good_P.mat','oldtable','newtable');
%%
sSWsbs = zeros(nSS,2+6+1);% predicted/measuredN/OldLearn/OTrans/Ttype/pr/ll/condition/subj
for i = 1:nSS
    sSWsbs(i,2) = newtable.beta(newtable.participant == SSubj(i));%real data new trans
    prtmp = newtable.PrLabel(newtable.participant == SSubj(i));
    chrow =  intersect(find(oldtable.participant == SSubj(i)),find(oldtable.PrLabel == prtmp))
    sSWsbs(i,3) = oldtable.beta(chrow(1));%real olddata learning%same subj 
    sSWsbs(i,4) = oldtable.beta(chrow(2));%real olddata transfer 
    sSWsbs(i,9) = SSubj(i);%subj 
end

pSWsbs = zeros(nPS,2+6+1);
for j = 1:nPS
    pSWsbs(j,2) = newtable.beta(newtable.participant == PSubj(j));%real data: trans new
    prtmp = newtable.PrLabel(newtable.participant == PSubj(j));
    chrow =  intersect(find(oldtable.participant == PSubj(j)),find(oldtable.PrLabel == prtmp))
    pSWsbs(j,3) = oldtable.beta(chrow(1));%real data
    pSWsbs(j,4) = oldtable.beta(chrow(2));%real data
    pSWsbs(j,9) = PSubj(j);%subj 
end
%%
%picking up the Pr and likelihood from the newtaable 
spWtr = zeros(nSS,1);
for  ii = 1:nSS
    sSWsbs(ii,5) = newtable.TLabel(newtable.participant == SSubj(ii));
    pr = newtable.PrLabel(newtable.participant == SSubj(ii));
    sSWsbs(ii,6) = pr; %col 3 pr
    ll = newtable.LLabel(newtable.participant == SSubj(ii));
    sSWsbs(ii,7) = ll; %col 4 ll
    lG = sLGamma(sLGamma(:,3) == SSubj(ii),ll);
    pG = sPrGamma_r(sPrGamma_r(:,3) == SSubj(ii),pr);
    if pG == nan
        spWtr(ii) = nan;
    else
        spWtr(ii) = (lG*5)/(lG*5+pG);
    end
end
sSWsbs(:,1) = spWtr; %col1 predicted by subjective prior and likelihood
sSWsbs(:,8) = repmat(1,size(sSWsbs,1),1); %col5 group serial  == 1 (rtSC)
[h1,p1] = ttest(sSWsbs(:,1),sSWsbs(:,2))
%%
ppWtr = zeros(nPS,1);
for  jj = 1:nPS
    pSWsbs(jj,5) = newtable.TLabel(newtable.participant == PSubj(jj));
    pr = newtable.PrLabel(newtable.participant == PSubj(jj));
    pSWsbs(jj,6) = pr; %col 3 pr
    ll = newtable.LLabel(newtable.participant == PSubj(jj));
    pSWsbs(jj,7) = ll; %col 4 ll
    lG = sLGamma(sLGamma(:,3) == PSubj(jj),ll);
    pG = sPrGamma_r(sPrGamma_r(:,3) == PSubj(jj),pr);
    if pG == nan
        ppWtr(jj) = nan;
    else
        ppWtr(jj) = (lG*5)/(lG*5+pG);
    end
end
pSWsbs(:,1) = ppWtr;
pSWsbs(:,8) = repmat(2,size(pSWsbs,1),1); %col5 %group parallel = 2 

stattPN = [sSWsbs;pSWsbs];%col1 predict col2 measured 
rnantPN = stattPN(sum(isnan(stattPN),2)==0,:);%remove nan 
% N: 44+47 = 91 %  

%statistics mixed_anova
% compPMData = rnantPN(:,1:2); %1prediction 2measure
% compNOData = [rnantPN(:,2) rnantPN(:,4)]; %2measure %4transfer
% between_factors1 = rnantPN(:,6:8);
% tbl1 = simple_mixed_anova(compPMData,between_factors1,{'Prediction'},{'Pr_Types','LL_Types','Cog_load'})
% tbl2 = simple_mixed_anova(compNOData,between_factors1,{'Experience'},{'Pr_Types','LL_Types','Cog_load'})
% 
% mean(rnantPN(rnantPN(:,7) == 1,2)) %narrow LL
% mean(rnantPN(rnantPN(:,7) == 2,2))
% 
% mean(rnantPN(rnantPN(:,6) == 1,2)) %narrow Pr
% mean(rnantPN(rnantPN(:,6) == 2,2))
%%
stattS = rnantPN(rnantPN(:,8)==1,1:7);
stattP = rnantPN(rnantPN(:,8)==2,1:7);
between_factorS = stattS(:,6:7); %pr and LL
between_factorP = stattP(:,6:7);

tblS = simple_mixed_anova(stattS(:,1:2),between_factorS,{'Prediction'},{'Pr_Types','LL_Types'});
tblP = simple_mixed_anova(stattP(:,1:2),between_factorP,{'Prediction'},{'Pr_Types','LL_Types'});

mean(stattS(stattS(:,6) == 1,1))
mean(stattS(stattS(:,6) == 2,1))
mean(stattS(stattS(:,6) == 1,2)) %M
mean(stattS(stattS(:,6) == 2,2))

mean(stattP(stattP(:,6) == 1,1))
mean(stattP(stattP(:,6) == 2,1))
mean(stattP(stattP(:,6) == 1,2)) %M
mean(stattP(stattP(:,6) == 2,2))
%% transfer score 
tSC = zeros(size(rnantPN,1),1);%transfer score 
pChg = (rnantPN(:,1)-rnantPN(:,4));
mChg = (rnantPN(:,2)-rnantPN(:,4));
tSC = ((rnantPN(:,2)-rnantPN(:,4))./(rnantPN(:,1)-rnantPN(:,4)));
tSCtable = [tSC,rnantPN(:,5:9),pChg,mChg];%5:9 =>  2:6 7:expected 8:measured
idx1 = isoutlier(tSC,'median'); %remove outliers %new_val
% idx1 = isoutlier(tSC,'movmedian',16); %remove outliers %new_val
% idx1 = isoutlier(tSC,'grubbs'); %remove outliers %new_val
rtSC = tSCtable(idx1 == 0,:);
% rtSC = tSCtable;
%%
swtest(rtSC(:,1))
med_rtSC = median(rtSC(:,1))

[psrank0, hsrank0] = signrank(rtSC(:,1),0,'tail','left')
[psrank1, ~] = signrank(rtSC(:,1),0,'tail','right')

iqr_rtsc = iqr(rtSC(:,1))
[p,h,stats] = signtest(rtSC(:,1),0,'Method','approximate','tail','right')
psrank0l = signrank(rtSC(:,1),0,'tail','right')

psrank1r = signrank(rtSC(:,1),1,'tail','right')
psrank1b = signrank(rtSC(:,1),1,'tail','both')
psrank1l = signrank(rtSC(:,1),1,'tail','left')
ranksum(rtSC(1:44,1),rtSC(45:end,1),'tail','right')

lrtSC = log(rtSC(:,1));
p = ranksum(rtSC(rtSC(:,5) == 1,1),rtSC(rtSC(:,5) == 2,1))
medS_ts = median(rtSC(rtSC(:,5) == 1,1))
medP_ts = median(rtSC(rtSC(:,5) == 2,1))
iqrS_ts = iqr(rtSC(rtSC(:,5) == 1,1))
iqrP_ts = iqr(rtSC(rtSC(:,5) == 2,1))
%%
bfData = [rtSC(:,1),rtSC(:,3:6)];
%1Wpl 2WpL 3WPl 4WPL 5OIpl 6OIpL 7OIPl 8PL
%9Group 10subj
varTypes = ["double","logical","logical","logical","double"];
varNames = ["rtSC","Pr","LL","CogLoad","participant"];
bftable = array2table(bfData,'VariableNames',varNames);

lt00 = fitlme(bftable,'rtSC~1+(1|participant)')
ltsP = fitlme(bftable,'rtSC~Pr+(1|participant)')
ltsL = fitlme(bftable,'rtSC~LL+(1|participant)')
ltsC = fitlme(bftable,'rtSC~CogLoad+(1|participant)')

lts0 = fitlme(bftable,'rtSC~Pr+LL+(1|participant)')
ltsm = fitlme(bftable,'rtSC~Pr*LL+(1|participant)')
lts3 = fitlme(bftable,'rtSC~CogLoad+Pr+LL+(1|participant)');
ltsi = fitlme(bftable,'rtSC~CogLoad*Pr*LL+(1|participant)');
% LLts03 = compare(lts0,lts3,'CheckNesting',true,'NSim',100) %p = .91
% LLts3i = compare(lts0,ltsi,'CheckNesting',true,'NSim',100) %p = .68
% LLtsm3 = compare(lts3,ltsm,'CheckNesting',true,'NSim',100) %p = .33
% LLtsmi = compare(ltsm,ltsi,'CheckNesting',true,'NSim',100) %p = .53
% LLts00 = compare(lt00,lts0,'CheckNesting',true,'NSim',100) %p = .73
LLts0P = compare(lt00,ltsP,'CheckNesting',true,'NSim',100) %val p = .73
LLts0L = compare(lt00,ltsL,'CheckNesting',true,'NSim',100) %val p = .73
LLts0C = compare(lt00,ltsC,'CheckNesting',true,'NSim',100) %val p = .73
LLts0I = compare(lt00,ltsi,'CheckNesting',true,'NSim',100) %val p = .73
LLtsPI = compare(ltsP,ltsi,'CheckNesting',true,'NSim',100) %val p = .73
%%
% writetable(bftable,'ps_TransScoreTable.csv','Delimiter',',','QuoteStrings',true)
% modelFull = fitlme(bftable,'rtSC~Pr*LL*CogLoad');
% modelM = fitlme(bftable,'rtSC~Pr+LL+CogLoad');

% Keep main of freq and ori:freq interaction.
%The evidence for the main effect is the ratio of the Bayes Factors.
%bfMain = bfFull/bfRestricted

%[bf10,pValue] =%bf.ttest2(rtSC(1:44,1),rtSC(45:end,1),'tail','right')%Trang
%% parametric 
newTS = (bftable.rtSC-mean(bftable.rtSC))/std(bftable.rtSC);
[hks,pks,~,~] = kstest(newTS)
[~,tbl_bm,stats_bm] = anovan(rtSC(:,1),{rtSC(:,3) rtSC(:,4), rtSC(:,5)},"Model","interaction", ...
    "Varnames",["PrType","LLType","Cog_load"]);

bftable.Pr = categorical(bftable.Pr);
bftable.LL = categorical(bftable.LL);
bftable.CogLoad1 = categorical(bftable.CogLoad);
modelFull = fitlme(bftable,'rtSC~Pr*LL*CogLoad1')
modelFull.anova
% [bfWithMain] = bf.anova(sleep,'extra~group+ID+group:ID','treatAsRandom',{'ID','group:ID'},'scale',1);
bfFull = bf.anova(bftable,'rtSC~Pr*LL');
% The Bayes Factor shows that the Full model as comapred to the Null model 
% (i.e. intercept only model)
bfFull
%% adapted from rouderfigure5.m
bfLLFixed  =   bf.anova(bftable,'rtSC~LL+Pr:LL');  % Keep main of freq and ori:freq interaction.
bfPrFixed  =   bf.anova(bftable,'rtSC~Pr+Pr:LL');  % Keep main of freq and ori:freq interaction.
% The evidence for the main effect is the ratio of the Bayes Factors. 
bfMain1 = bfFull/bfLLFixed %evidence of Pr main effect%validation
bfMain2 = bfFull/bfPrFixed %evidence of LL main effect
%%
m_rtSC = mean(rtSC(:,1))
se_rtSC = sqrt(std(rtSC(:,1))^2/size(rtSC,1))
[hrtSC, prtSC]= ttest(rtSC(:,1),0,"Tail","right") %if significantly larger than 0
[hrtSCm1, prtSCm1] = ttest(rtSC(:,1),1,"Tail","left") %if significantly different from 1

[results13,~,~,gnames13] = multcompare(stats_bm,"Dimension",[1 3]);%Pr Cog
[results23,~,~,gnames23] = multcompare(stats_bm,"Dimension",[2 3]);%LL Cog

sPrnId = intersect(find(rtSC(:,5) == 1),find(rtSC(:,3)== 1));
sPrwId = intersect(find(rtSC(:,5) == 1),find(rtSC(:,3)== 2));
pPrnId = intersect(find(rtSC(:,5) == 2),find(rtSC(:,3)== 1));
pPrwId = intersect(find(rtSC(:,5) == 2),find(rtSC(:,3)== 2));

snrtsc = rtSC(sPrnId,1);swrtsc = rtSC(sPrwId,1);
pnrtsc = rtSC(pPrnId,1);pwrtsc = rtSC(pPrwId,1);
t1 = [snrtsc;swrtsc]%all serial 
t2 = [pnrtsc;pwrtsc]%all parallel
se_t1 = sqrt(std(t1)^2/size(t1,1))
se_t2 = sqrt(std(t2)^2/size(t2,1))
[~,p]= ttest2(t1,t2)
[bf10tPS,ptPS] = bf.ttest2(t1,t2,'tail','both')
mean(t1)
mean(t2)
mean(snrtsc)
mean(swrtsc)
mean(pnrtsc)
mean(pwrtsc)
se_pn = sqrt(std(pnrtsc)^2/size(pnrtsc,1))
se_pw = sqrt(std(pwrtsc)^2/size(pwrtsc,1))

[~,ptsc_PrS] = ttest2(snrtsc,swrtsc)
[~,ptsc_PrP]= ttest2(pnrtsc,pwrtsc)

bfData = [rtSC(:,1),rtSC(:,3:5)-1,rtSC(:,6)];
%%
sfacec = [0.9 0.59 0.33;0.3 0.2 0.11];

figure; %Weights
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m, n, p, [0.0 0.01], [0.1 0.01], [0.2 0.2]);
if ~make_it_tight,  clear subplot;  end

Titles = ["Serial","Parallel"];
sidx = find(rnantPN(:,6) == 1);
% Titles = ["Group1 (TransferNew)","Group2 (TransferNew)"];
for kk  = 1:2 %group serial versus parallel
    if kk  == 1%1 %p
        uidx = find(rnantPN(:,8) == 1);
    elseif kk == 2% 2 %S
        uidx = find(rnantPN(:,8) == 2);
    end
    use = rnantPN(uidx,:);

    subplot(1,2,kk)
    hold on;
    bsp = boxchart(use(:,5)*2-0.32,use(:,1));%predicted
    bsp.BoxFaceColor =[.75 .75 .75];% '#EFEF21';%'#A9A9AD';%bNcolor(2);
    bsm = boxchart(use(:,5)*2+0.32,use(:,2));
    bsm.BoxFaceColor ='#B3AF05';%'#EFEF21';%'#A9A9AD'; %bNcolor(1);
    xlabel('trial types');
    hold on;

    for jj =1:max(rnantPN(:,5)) %4 types
        for mm = 1:2 %1: predicted 2: measured
            if mm == 1
                y1p = use((use(:,5) == jj),1); %col1: predicted
                x1p=repmat(jj*2-0.32,size(y1p,1),1); %col2: measured
                scatter(x1p(:),y1p(:),[],[.5 .5 .5],...
                    'filled','MarkerFaceAlpha',0.7','jitter','on','jitterAmount',0);
            else
                y1m = use((use(:,5) == jj),2);
                x1m=repmat(jj*2+0.32,size(y1p,1),1);
                if mod(jj,2) == 1
                    scatter(x1m(:),y1m(:),[],sfacec(1,:),...
                        'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0);
                else mod(jj,2) == 0
                    scatter(x1m(:),y1m(:),[],sfacec(2,:),...
                        'filled','MarkerFaceAlpha',0.4','jitter','on','jitterAmount',0);
                end
            end
        end
        pltx = [x1p x1m];
        plty = [y1p y1m];
        for nn = 1:size(pltx,1)
            plot(pltx(nn,:),plty(nn,:),"Color",[.65 .65 .65],"LineStyle","-")
        end
    end

% ylim([0 5.5]);
ylim([0 1.18]);
xticks([2 4 6 8]);
xticklabels({'PnLn','PnLw','PwLn','PwLw'});

if kk == 2 
    legend([bsp bsm],{'Predicted','Measured'})%'Location','northeastoutside'
    set(gca,'yticklabels',[]);
else
    ylabel('likelihood weight');
end 

title(Titles(kk));
box on;
axis square;
grid on; 
ax = gca;
ax.FontSize = 16;
end
%% transfer score plots 
remove = 1;
xposi = [2 4 6 8];
figure;
% t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
Titles = ["Serial","Parallel"];
if remove == 1
ptSC = rtSC;
else
ptSC = tSCtable;
end 
sjN = nan(2,4);
for kl  = 1:2 %group serial versus parallel
    if kl  == 1%serial
        scidx = find(ptSC(:,5) == 1);%score index
    elseif kl == 2%parallel
        scidx = find(ptSC(:,5) == 2);
    end
    use2 = ptSC(scidx,:);

    subplot(1,2,kl)
    hold on;
%     bsc = boxchart(use2(:,2)*2,use2(:,1));%predicted
%     bsc.BoxFaceColor =[.5 .5 .5];% '#EFEF21';%'#A9A9AD';%bNcolor(2);
    xlabel('trial types')
    hold on;

    for jk =1:max(rtSC(:,2)) %4 types
        y1p = use2((use2(:,2) == jk),1); %col1: predicted
        x1p =  repmat(jk*2,size(y1p,1),1); %col2: measured
        ypi(jk) = violin(y1p,'x',[xposi(jk)],'facecolor',[.6 .6 .6],...
            'edgecolor','k','plotlegend',0);
        scatter(x1p(:),y1p(:),22,[0.4 0.2 0.11],...
            'filled','MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0);
        sjN(kl,jk) = size(y1p,1);
    end

    asy_p = 1.7 

    yline(0,'--','LineWidth',1.5)
    yline(1,'--','LineWidth',1.5)

    xlim([1 9]);
%     ylim([-10 10.5]);
    ylim([-2.1 3.2]);
    xticks([2 4 6 8]);
    xticklabels({'PnLn','PnLw','PwLn','PwLw'});
    if kl == 1
        ylabel('transfer score (a.u.)');
    else
        set(gca,'yticklabels',[]);
    end

    title(Titles(kl));
    box on;
    axis square;
    grid on;
    ax = gca;
    ax.FontSize = 18;
end
%%
figure; 
xpos = [1 2];
for ij =1:2
    yline(0,'--','LineWidth',2.5)
    yline(1,'--','LineWidth',2.5)
    y1v = bftable.rtSC(bftable.CogLoad == ij) %-1
    x1p = repmat(xpos(ij),size(y1v,1),1)        
    vp(ij) = violin(y1v,'x',[xpos(ij)],'facecolor',[.6 .6 .6],...
            'edgecolor','k','plotlegend',0);
    erp = errorbar(xpos(ij),median(y1v),...
        quantile(y1v,0.25)-median(y1v),quantile(y1v,0.75)-median(y1v));
%     erp = errorbar(xpos(ij),mean(y1v),...
%         std(y1v)/sqrt(size(y1v,1)));
    erp.Color = [0.1 0.1 0.1];
    erp.LineStyle = 'none';
    erp.LineWidth = 4.5;
    s(ij) = scatter(x1p(:),y1v,48,[0.4 0.2 0.11],...
        'filled','MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15);
end

xlim([0.5 2.5]); %0.3 5
ylim([-1.5 2.5]);
yticks([0 1]);

ylabel('transfer score (a.u.)');
box off; 
grid on;
grid minor
ax = gca;
ax.FontSize = 45; %48
ax.LineWidth = 1;
set(gca,'xticklabel',{[]})