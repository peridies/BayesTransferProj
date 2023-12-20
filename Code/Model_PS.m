%% compare real data with modelled data 
% 1. linear model (1)mixed linear (2)likelihood-only
% 2. exampler
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
    'nPS')%,

load('OldNewtable_Good_P.mat','oldtable','newtable','qtable');
Roldtable = oldtable; 
Rnewtable = newtable;
Rqtable = qtable;
clear oldtable newtable qtable 
load('OldNewtable_Exemplar_P.mat','oldtable','newtable','qtable');
Moldtable = oldtable; 
Mnewtable = newtable;
Mqtable = qtable;
clear oldtable newtable qtable 
%%
sSWsbs = zeros(nSS,16);
% 1predB/2predL/3predDF/4predE/5measuredNew/6OldLearn/7OTrans
% /8TLabel/9Pr/10LL/11Grp/12d1/13d2/14d3/15d4/16subj
for i = 1:nSS
    sSWsbs(i,5) = Rnewtable.beta(Rnewtable.participant == SSubj(i));%real data new trans
    prtmp = Rnewtable.PrLabel(Rnewtable.participant == SSubj(i));
    chrow =  intersect(find(Roldtable.participant == SSubj(i)),find(Roldtable.PrLabel == prtmp))
    sSWsbs(i,6) = Roldtable.beta(chrow(1));%real olddata learning%same subj
    sSWsbs(i,7) = Roldtable.beta(chrow(2));%real olddata transfer
    sSWsbs(i,16) = SSubj(i);%subj
    sSWsbs(i,4) = Mnewtable.beta(Mnewtable.participant == SSubj(i));
end

pSWsbs = zeros(nPS,16);
for j = 1:nPS
    pSWsbs(j,5) = Rnewtable.beta(Rnewtable.participant == PSubj(j));%real data: trans new
    prtmp = Rnewtable.PrLabel(Rnewtable.participant == PSubj(j));
    chrow =  intersect(find(Roldtable.participant == PSubj(j)),find(Roldtable.PrLabel == prtmp))
    pSWsbs(j,6) = Roldtable.beta(chrow(1));%real data learn
    pSWsbs(j,7) = Roldtable.beta(chrow(2));%real data trans
    pSWsbs(j,16) = PSubj(j);%subj
    pSWsbs(j,4) = Mnewtable.beta(Mnewtable.participant == PSubj(j));%real data: trans new
end
%% picking up the Pr and likelihood from the newtaable
all = [1 2];
spWtr = zeros(nSS,1);
spWSlM = zeros(nSS,1); %predicted by simple slope model
spWSDM = zeros(nSS,1); %predicted by direct fit model
for  ii = 1:nSS
    thisSb = SSubj(ii);
    sSWsbs(ii,8) = Rnewtable.TLabel(Rnewtable.participant == thisSb);
    prn = Rnewtable.PrLabel(Rnewtable.participant == thisSb);%new pr
    sSWsbs(ii,9) = prn; %col 3 pr
    lln = Rnewtable.LLabel(Rnewtable.participant == thisSb);%new ll
    sSWsbs(ii,10) = lln; %col 4 ll
    lGn = sLGamma(sLGamma(:,3) == thisSb,lln);
    pGn = sPrGamma_r(sPrGamma_r(:,3) == thisSb,prn);
    ol1 = setdiff(all,lln);
    op1 = prn; %one old pair has to be newPr/oldLL
    ol2 = lln; %the other is orthogonal to 1 (oldPr/newLL)
    op2 = setdiff(all,prn);%the other is orthogonal to 1
    lgo1 = sLGamma(sLGamma(:,3) == thisSb,ol1);
    pgo1 = sPrGamma_r(sPrGamma_r(:,3) == thisSb,op1);
    lgo2 = sLGamma(sLGamma(:,3) == thisSb,ol2);
    pgo2 = sPrGamma_r(sPrGamma_r(:,3) == thisSb,op2);
    ob1 = sSWsbs(sSWsbs(:,16) == thisSb,6);
    ob2 = Rqtable.beta(Rqtable.participant == thisSb);
    lgo = [lgo1;lgo2];
    pgo = [pgo1;pgo2];
    ob = [ob1;ob2];
    if isnan(pGn)
        spWtr(ii) = nan;
    else
        spWtr(ii) = (lGn*5)/(lGn*5+pGn);
    end

    spWSlM(ii) = ob1*lGn/lgo1;
    if spWSlM(ii)>1
        spWSlM(ii)=1;
    end

    if isnan(pgo)
        spWSDM(ii) = nan;
    else
        spWSDM(ii) = directfit(lgo,pgo,lGn,pGn,ob); %predicted by direct fit model
        if spWSDM(ii)>1
            spWSDM(ii)=1;
        elseif spWSDM(ii)<0
            spWSDM(ii)=0;
        end
    end
end
sSWsbs(:,1) = spWtr; %col1 predicted by subjective prior and likelihood
sSWsbs(:,2) = spWSlM;
sSWsbs(:,3) = spWSDM;
sSWsbs(:,11) = repmat(1,size(sSWsbs,1),1); %col5 group serial  == 1 (rtSC)
[h1,p1] = ttest(sSWsbs(:,1),sSWsbs(:,2))
%%
ppWtr = zeros(nPS,1);
ppWSlM = zeros(nPS,1);
ppWSDM = zeros(nPS,1); %predicted by direct fit model
for  jj = 1:nPS
    thisSb = PSubj(jj);
    pSWsbs(jj,8) = Rnewtable.TLabel(Rnewtable.participant == thisSb);
    prn = Rnewtable.PrLabel(Rnewtable.participant == thisSb);%new pr
    pSWsbs(jj,9) = prn; %col 3 pr
    lln = Rnewtable.LLabel(Rnewtable.participant == thisSb);%new ll
    pSWsbs(jj,10) = lln; %col 4 ll
    lGn = sLGamma(sLGamma(:,3) == thisSb,lln);
    pGn = sPrGamma_r(sPrGamma_r(:,3) == thisSb,prn);
    ol1 = setdiff(all,lln);
    op1 = prn; %one old pair has to be newPr/oldLL
    ol2 = lln; %the other is orthogonal to 1 (oldPr/newLL)
    op2 = setdiff(all,prn);%the other is orthogonal to 1
    lgo1 = sLGamma(sLGamma(:,3) == thisSb,ol1);
    pgo1 = sPrGamma_r(sPrGamma_r(:,3) == thisSb,op1);
    lgo2 = sLGamma(sLGamma(:,3) == thisSb,ol2);
    pgo2 = sPrGamma_r(sPrGamma_r(:,3) == thisSb,op2);
    ob1 = pSWsbs(pSWsbs(:,16) == thisSb,6);
    ob2 = Rqtable.beta(Rqtable.participant == thisSb);
    lgo = [lgo1;lgo2];
    pgo = [pgo1;pgo2];
    ob = [ob1;ob2];
    if isnan(pGn)
        ppWtr(jj) = nan;
    else
        ppWtr(jj) = (lGn*5)/(lGn*5+pGn);
    end

    ppWSlM(jj) = ob1*lGn/lgo1;
    if ppWSlM(jj)>1
        ppWSlM(jj)=1;
    end

    if isnan(pgo)
        ppWSDM(jj) = nan;
    else
        ppWSDM(jj) = directfit(lgo,pgo,lGn,pGn,ob); %predicted by direct fit model
        if ppWSDM(jj)>1
            ppWSDM(jj)=1;
        elseif ppWSDM(jj)<0
            ppWSDM(jj)=0;
        end
    end
end
pSWsbs(:,1) = ppWtr;
pSWsbs(:,2) = ppWSlM;
pSWsbs(:,3) = ppWSDM;
pSWsbs(:,11) = repmat(2,size(pSWsbs,1),1); %col5 %group parallel = 2

stattPN = [sSWsbs;pSWsbs];%col1 predict col2 measured
rnantPN = stattPN(sum(isnan(stattPN),2)==0,:);%remove nan N: 45+47 = 92

rnantPN(:,12) = abs(rnantPN(:,5)-rnantPN(:,1)); %abs error data-bayes prediction
rnantPN(:,13)= abs(rnantPN(:,5)-rnantPN(:,2)); %abs data-heuristic1:LL only 
rnantPN(:,14)= abs(rnantPN(:,5)-rnantPN(:,3)); %abs data-heuristic2:direct fit
rnantPN(:,15)= abs(rnantPN(:,5)-rnantPN(:,4)); %abs data-heuristic2:direct fit


swtest(rnantPN(:,12))
swtest(rnantPN(:,13))
signrank(rnantPN(:,12),rnantPN(:,13),'tail','left')
signrank(rnantPN(:,12),rnantPN(:,14),'tail','left')
signrank(rnantPN(:,12),rnantPN(:,15),'tail','left')
%% compare models
nSample = size(rnantPN,1) %col 4 real measured data 
LLb = nSample*(log(sum((rnantPN(:,5)-rnantPN(:,1)).^2)/nSample)) %Bayes LL
LLl = nSample*(log(sum((rnantPN(:,5)-rnantPN(:,2)).^2)/nSample))%heuristic LL
LLd = nSample*(log(sum((rnantPN(:,5)-rnantPN(:,3)).^2)/nSample)) %heuristic direct
LLm = nSample*(log(sum((rnantPN(:,5)-rnantPN(:,4)).^2)/nSample)) %exemplar

BICb = LLb
BICl = LLl
BICd = 3*log(nSample) + LLd
BICm = LLm
%%
sfacec = [0.9 0.59 0.33;0.3 0.2 0.11];
c_1 = [175 233 221]/255; %old [0.1 0.3 0.5] old_1 [0.35 0.1 0.45]
c_2 = [0,0.2,1]; %old [0.4 0.1 0.5]
c_3 = [0.5,0.5,0.5];
c_4 = [0 204 255]/255;
c_5 = ([141 95 211]/255)*0.7;

figure; %Weights
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m, n, p, [0.0 0.01], [0.1 0.01], [0.2 0.2]);
if ~make_it_tight,  clear subplot;  end

Titles = ["serial","parallel"];
xpm = [1.1 1.9 2.7 3.5];
for kk  = 1:2 %group serial versus parallel
    if kk  == 1%1 %serial
        uidx = find(rnantPN(:,10) == 1);
    elseif kk == 2% 2 %parallel
        uidx = find(rnantPN(:,10) == 2);
    end
    use = rnantPN(uidx,:);

    subplot(1,2,kk)
    hold on;
        y1v = use(:,12);
        y2v = use(:,13);
        y3v = use(:,14);
        y4v = use(:,15);
        v1 = violin(y1v,'x',[xpm(1)],'facecolor',c_3,...
            'edgecolor',[.3 .3 .3],'plotlegend',0);%predicted
        erv1 = errorbar(xpm(1),median(y1v),...
        quantile(y1v,0.25)-median(y1v),quantile(y1v,0.75)-median(y1v));
        v2 = violin(y2v,'x',[xpm(2)],'facecolor',c_4,...
            'edgecolor',[.5 .5 .5],'plotlegend',0);%predicted
        erv2 = errorbar(xpm(2),median(y2v),...
        quantile(y2v,0.25)-median(y2v),quantile(y2v,0.75)-median(y2v));
        v3 = violin(y3v,'x',[xpm(3)],'facecolor',c_1,...
            'edgecolor',[.5 .5 .5],'plotlegend',0);%predicted
        erv3 = errorbar(xpm(3),median(y3v),...
        quantile(y3v,0.25)-median(y3v),quantile(y3v,0.75)-median(y3v));
        v4 = violin(y4v,'x',[xpm(4)],'facecolor',c_5,...
            'edgecolor',[.5 .5 .5],'plotlegend',0);%predicted
        erv4 = errorbar(xpm(4),median(y4v),...
        quantile(y4v,0.25)-median(y4v),quantile(y4v,0.75)-median(y4v));
        
        erv1.Color = [0.1 0.1 0.1];
        erv1.LineStyle = 'none';
        erv1.LineWidth = .7;
        erv2.Color = [0.1 0.1 0.1];
        erv2.LineStyle = 'none';
        erv2.LineWidth = .7;
        erv3.Color = [0.1 0.1 0.1];
        erv3.LineStyle = 'none';
        erv3.LineWidth = .7;
        erv4.Color = [0.1 0.1 0.1];
        erv4.LineStyle = 'none';
        erv4.LineWidth = .7;

        x1v=repmat(xpm(1),size(y1v,1),1);
        s1 = scatter(x1v(:),y1v(:),[],sfacec(2,:),...
            'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0.04);
        x2v=repmat(xpm(2),size(y2v,1),1);
        s2 = scatter(x2v(:),y2v(:),[],sfacec(2,:),...
            'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0.04);
        x3v=repmat(xpm(3),size(y3v,1),1);
        s3 = scatter(x3v(:),y3v(:),[],sfacec(2,:),...
            'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0.04);
        x4v=repmat(xpm(4),size(y4v,1),1);
        s4 = scatter(x4v(:),y4v(:),[],sfacec(2,:),...
            'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0.04);
        xticks([1.1 1.9 2.7 3.5]+0.2);
        xticklabels({'Bayes','linear','likelihood','exemplar'});
        
        hold on;

    xlim([0.5 4.1]);
    ylim([-0.5 1.5]);
    yline(0)

    if kk == 2
        set(gca,'yticklabels',[]);
    else
        ylabel('|data - model prediction|');
    end

    title(Titles(kk));
    box on;
    axis square;
    grid on;
    grid minor;
    ax = gca;
    ax.FontSize = 24;
end
%% plot BIC results
BIC = [BICb,BICd,BICl BICm];
dBIC = BIC-(min(BIC))
figure; 
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m, n, p, [0.0 0.01], [0.15 0.01], [0.2 0.2]);
if ~make_it_tight,  clear subplot;  end


subplot(1,2,1);hold on; 
bB = bar(1,BICb,0.3,'FaceColor',c_3)%'EdgeColor',[1 1 1])%extrapolation
bd = bar(2,BICd,0.3,'FaceColor',c_4)%'EdgeColor',[1 1 1])%interpolation 
bl = bar(3,BICl,0.3,'FaceColor',c_1)%,'EdgeColor',[1 1 1])
bm = bar(4,BICm,0.3,'FaceColor',c_5,'FaceAlpha',0.5)%,'EdgeColor',[1 1 1])

xlim([0.5 4.5]);
% ylim([-450 205])
% ylim([-500 228])
ylim([-300 137])
xticks([1 2 3 4]+0.2);
xticklabels({'Bayes','linear','likelihood','exemplar'});
ax = gca;
ax.FontSize = 18;
% axis square;
grid on; 
grid minor;