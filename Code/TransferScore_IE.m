%% thus we compute expected optimal weighting using Gamma_r(raw)
%Kiryakova 2020 equation 
%using subjective likelihood and priors
%compare with measured SW to create transfer score 
%house keeping
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_validation\'

cd(dir);
%%
load('subjective_variances1_P.mat'); %subjective variances using  root mean square instead of SD
load('Learning4Stats_Good_P.mat','exSubj','inSubj','nExS',...
    'nInS');

load('Transfer4Stats_Good.mat','DataTable'); %already MAD outlier removed 

exSWsbs = zeros(nExS,2+4);%col1 predicted col2 measured ex: Lw %col3 to compare use Ln(the ll that both had)
for i = 1:nExS
    exSWsbs(i,2) = DataTable.Weights(i*3);%real data Lw (trans)
    exSWsbs(i,3) = DataTable.Weights(i*3-2);%real data Ln
    exSWsbs(i,4) = DataTable.Weights(i*3-1);%real data Lm
    exSWsbs(i,5) = 2;
    exSWsbs(i,6) = DataTable.Subjects(i*3); 
end

inSWsbs = zeros(nInS,2+4);
for j = 1:nInS
    inSWsbs(j,2) = DataTable.Weights((j+nExS)*3-1);%real data Lm (trans)
    inSWsbs(j,3) = DataTable.Weights((j+nExS)*3-2);%real data Ln
    inSWsbs(j,4) = DataTable.Weights((j+nExS)*3);%real data Lw
    inSWsbs(j,5) = 1;
    inSWsbs(j,6) = DataTable.Subjects(j*3); 
end
%% compute expected SW
idxe = zeros(nExS,1);
idxi = zeros(nInS,1);
%extra
expWtr = zeros(nExS,1); %expected W of extrapolation 
for  ii = 1:nExS
expWtr(ii) = (sLGamma(ii,3)*5)/(sLGamma(ii,3)*5+msPrGamma_r(ii));
end
exSWsbs(:,1) = expWtr; %col1 predicted slope
[nani1,~]= find(isnan(exSWsbs(:,1)));%exclude nan
% idxe = isoutlier(exSWsbs(:,2),'median'); %remove outliers
idxe(nani1) = 1;
texSWsbs = exSWsbs(idxe == 0,:); %g = good:not nan  
% iswO1 = isoutlier(texSWsbs(:,2),'median');%exclude performance outliners 
% gexSWsbs = texSWsbs(iswO1 == 0,:); %g = good:median 
gexSWsbs = texSWsbs; %g = good:median 

%inter 
inpWtr = zeros(nInS,1);
for  jj = 1:nInS
inpWtr(jj) = (sLGamma(jj+nExS,2)*5)/(sLGamma(jj+nExS,2)*5+msPrGamma_r(jj+nExS));
end
inSWsbs(:,1) = inpWtr;
[nani2,~]= find(isnan(inSWsbs(:,1)));
idxi(nani2) = 1;
tinSWsbs = inSWsbs(idxi == 0,:);%not nan

ginSWsbs = tinSWsbs; %g = good:median 

[h1,p1] = ttest(gexSWsbs(:,1),gexSWsbs(:,2));

tSWtable =  [ginSWsbs;gexSWsbs]; %i.e. non-nan 
isw = isoutlier(tSWtable(:,2),'median');%exclude performance outliners 
SWtable = tSWtable; %g = nan removed 
%% transfer score  
tSC = zeros(size(SWtable,1),2);%transfer score Col1 predicted
tmpu1 = SWtable(:,2)-SWtable(:,3); %measured changes
tmpl1 = SWtable(:,1)-SWtable(:,3); %predicted changes 
tSC(:,1) = tmpu1./tmpl1; %scores measure/predicted Ln

tmpu2 = SWtable(:,2)-SWtable(:,4); %measured changes
tmpl2 = SWtable(:,1)-SWtable(:,4); %predicted changes 
tSC(:,2) = tmpu2./tmpl2; %scores predicted from Lw(in) Lm(ex)

v1 = [min(tSC(:)),min(tSC(:))];%line with slope = 1
v2 = [max(tSC(:)),max(tSC(:))];
pt = [tSC];
distance_2D = point_to_line_distance(pt, v1, v2);
% itsDall = isoutlier(distance_2D,'median');%discrenpancies outliners using 3MAD
itsDall = isoutlier(distance_2D,'movmedian',16);%discrenpancies outliners using 3MAD
itsM = isoutlier(distance_2D,'mean');%exclude outlinersusing 3 std
itsG = isoutlier(distance_2D,'grubbs');%exclude outlinersusing Grub
tSCtable = [mean(tSC,2),SWtable(:,5:6)];%1-31 and then end tSCtable(itsd == 0,:);%exclude MAD outliers 
itsD1 = isoutlier(distance_2D(tSCtable(:,2)== 1),'median');%exclude outliners using 3MAD
itsD2 = isoutlier(distance_2D(tSCtable(:,2)== 2),'median');%exclude outliners using 3MAD
itsD = [itsD1;itsD2];
its12 = isoutlier(tSCtable(:,1),'median');%exclude outliners %not in use

cTS = [tSC,its12,distance_2D,itsDall,itsD,itsM];%col5:itsDall
tp_strict = sum(cTS(:,[5 6]),2);
tp_strict(tp_strict==1) = 0;
cTS(:,8) = tp_strict;
tmp_rtSC = tSCtable(itsDall == 0,:);%exclude discrepancies outliner
% tmp_rtSC = tSCtable(itsD == 0,:);%exclude discrepancies outliner
% tmp_rtSC = tSCtable(cTS(:,8) == 0,:);%exclude discrepancies outliner %tSCtable already took the mean 
%cTS col8 is the most liberal discrepancies exclusion criterion
its1 = isoutlier(tmp_rtSC(:,1),'median');%exclude value outliners 
% its1 = isoutlier(tmp_rtSC(:,1),'movmedian',30); %remove outliers %new_val
% rtSC = tmp_rtSC(its1 == 0,:);%exclude discrepancies outliner
rtSC = tmp_rtSC;%
%%
[H, pValue, ~] = swtest(rtSC(:,1))%discovery: normal, validation: normal
in_rtSC= (rtSC(rtSC(:,2)==1,1));
ex_rtSC= (rtSC(rtSC(:,2)==2,1));
[h_swrt,p_swrt,swrtstats]= swtest(in_rtSC(:,1))
[h_swrt,p_swrt,swrtstats]= swtest(ex_rtSC(:,1))

Med1ts = median(rtSC(rtSC(:,2)==1,1))
Med2ts = median(rtSC(rtSC(:,2)==2,1))
iqr1ts = iqr(rtSC(rtSC(:,2)==1,1))
iqr2ts = iqr(rtSC(rtSC(:,2)==2,1))
p_rs = ranksum(rtSC(rtSC(:,2)==1,1),rtSC(rtSC(:,2)==2,1))
p_rs1 = ranksum(rtSC(rtSC(:,2)==1,1),rtSC(rtSC(:,2)==2,1),'tail','right')

signrank(in_rtSC(:,1),1,'Tail','left')
signrank(ex_rtSC(:,1),1,'Tail','left')
signrank(in_rtSC(:,1),0,'Tail','right')
signrank(ex_rtSC(:,1),0,'Tail','right')

[p,h,stats] = signtest(rtSC(:,1),0,'Tail','right')
psrank1 = signrank(rtSC(:,1),1,'tail','right')
[~,ptie] = ttest2(in_rtSC,ex_rtSC,'tail','right')
p = ranksum(in_rtSC,ex_rtSC,'tail','right')

newTS = (rtSC(:,1)-mean(rtSC(:,1)))/std(rtSC(:,1));
[hks,pks,~,~] = kstest(newTS)
% rtSC = tSCtable;
varTypes = ["double","logical","double"];
varNames = ["rtSC","CogLoad","participant"];
rtSCTable = array2table(rtSC,'VariableNames',varNames);

% writetable(bftable,'ps_TransScoreTable.csv','Delimiter',',','QuoteStrings',true)
% [~,~,stats] = anovan(rtSC(:,1),{rtSC(:,2)},"Model","interaction", ...
%     "Varnames",["In_Ex"]);

in = size(find(rtSC(:,2)==1),1)
mean(rtSC(1:in,1))%
mean(rtSC(in+1:end,1))
SE_rtSCi = std(rtSC(1:in,1))/sqrt(size(rtSC(1:in,1),1))
SE_rtSCe = std(rtSC(in+1:end,1))/sqrt(size(rtSC(in+1:end,1),1))

in = length(find(rtSC(:,2)== 1));
[hsc,psc,cc,statsc]= ttest2(rtSC(1:in,1),rtSC(in+1:end,1))
[ht2,pt2,ci,stats] = ttest2(rtSC(1:in,1),rtSC(in+1:end,1),"Tail","right")

[~,ptsi1] = ttest(rtSC(1:in,1),1,'Tail','left')
[~,ptse1]=ttest(rtSC(in+1:end,1),1)
[~,ptse1]=ttest(rtSC(in+1:end,1),1,'tail','left')
[~,ptsi0]=ttest(rtSC(1:in,1),0,'tail','right')
[~,ptse0]=ttest(rtSC(in+1:end,1),0,'tail','right')

[hscm1,pscm1] = ttest(rtSC(:,1),1,'tail','left')
[hscm0,pscm0] = ttest(rtSC(:,1),0,'tail','right')

[bf10,pValue] = bf.ttest2(rtSC(1:in,1),rtSC(in+1:end,1),'tail','right')%
[bf10t,pValt] = bf.ttest2(rtSC(1:in,1),rtSC(in+1:end,1))%

[bf10_all,pValue_all] = bf.ttest(rtSC(:,1),1,'tail','left')%
[bf10_one,pValue_one] = bf.ttest(rtSC(1:in,1),1,'tail','left')

[hallsc,pallsc] = ttest(rtSC(:,1))
[hallsc1,pallsc1] = ttest(rtSC(:,1)-1)

psrank0R = signrank(rtSC(:,1),0,'tail','right')
psrank1B = signrank(rtSC(1:in,1),1,'tail','both')
psrank1B = signrank(rtSC(in+1:end,1),1,'tail','both')

psrank0L = signrank(rtSC(:,1),1,'tail','left')

[p_0,h_0,stats_0] = signtest(rtSC(:,1),0,'tail','right')
psrank1 = signrank(rtSC(:,1)-1,0)
%%
MSC = mean(rtSC(:,1))
medianTS = median(rtSC(:,1))
msrtsc = median(rtSC(1:in,1))
mprtsc = median(rtSC(in+1:end,1))
seSCs = sqrt((std(rtSC(1:in,1)))^2/in)
seSCp = sqrt((std(rtSC(in+1:end,1)))^2/(size(rtSC,1)-in))
%%prepare for the figure 5
aa = 1;
gnExS =size(gexSWsbs,1);
gnInS =size(ginSWsbs,1);
%%test normality of predicted SW
newgEb = (gexSWsbs(:,1)-mean(gexSWsbs(:,1)))/std(gexSWsbs(:,1));
[heks,peks,~,~] = kstest(newgEb) %discovery: Gaussian
newgIb = (ginSWsbs(:,1)-mean(ginSWsbs(:,1)))/std(ginSWsbs(:,1));
[hiks,piks,~,~] = kstest(newgIb) %discovery: not Gaussian
[psre,hsre] = signrank(gexSWsbs(:,1),gexSWsbs(:,2))%1:predicted%2:measured
[psri,hsri] = signrank(ginSWsbs(:,1),ginSWsbs(:,2))

xin = repmat(aa,size(ginSWsbs,1),2);
xin(:,1) = xin(:,1)-.18;
xin(:,2) = xin(:,2)+.18;
xex = repmat(aa*2,size(gexSWsbs,1),2);
xex(:,1) = xex(:,1)-.18;
xex(:,2) = xex(:,2)+.18;
mexpW = median(gexSWsbs(:,1));
mexrW = median(gexSWsbs(:,2));
err_expW = sqrt(std(gexSWsbs(:,1)).^2/size(gexSWsbs,1));
err_exrW = sqrt(std(gexSWsbs(:,2)).^2/size(gexSWsbs,1));
minpW = median(ginSWsbs(:,1));
minrW = median(ginSWsbs(:,2))
err_inpW = sqrt(std(ginSWsbs(:,1)).^2/size(ginSWsbs,1));
err_inrW = sqrt(std(ginSWsbs(:,2)).^2/size(ginSWsbs,1));
%%
figure; 
hold on;
bep = bar(2-0.18,mexpW,0.3,'FaceColor',[.5 .5 .5],'EdgeColor',[1 1 1])%extrapolation
bem = bar(2+0.18,mexrW,0.3,'FaceColor','#B3AF05','EdgeColor',[1 1 1])
bip = bar(1-0.18,minpW,0.3,'FaceColor',[.5 .5 .5],'EdgeColor',[1 1 1])%interpolation 
bim = bar(1+0.18,minrW,0.3,'FaceColor','#B3AF05','EdgeColor',[1 1 1])


erip = errorbar(1-0.18,minpW,err_inpW,err_inpW);
erip.Color = [.2 .2 .2];                            
erip.LineStyle = 'none';
erip.LineWidth = 1.5;
erir = errorbar(1+0.18,minrW,err_inrW,err_inrW);
erir.Color = [.2 .2 .2];                            
erir.LineStyle = 'none'; 
erir.LineWidth = 1.5;

erep = errorbar(2-0.18,mexpW,err_expW,err_expW);
erep.Color = [.2 .2 .2];                            
erep.LineStyle = 'none'; 
erep.LineWidth = 1.5;
erer = errorbar(2+0.18,mexrW,err_exrW,err_exrW);
erer.Color = [.2 .2 .2];                            
erer.LineStyle = 'none'; 
erer.LineWidth = 1.5;

for mm = 1:gnInS
    plot(xin(mm,:),ginSWsbs(mm,1:2),"Color",[.55 .55 .55],"LineStyle","-")
end
% asy_in = max(max(ginSWsbs(:,1:2)))+.02;
% text(1-.01,asy_in,'*','FontSize',20);

for nn = 1:gnExS
    plot(xex(nn,:),gexSWsbs(nn,1:2),"Color",[.55 .55 .55],"LineStyle","-")
end 
asy_ex = max(max(gexSWsbs(:,1:2)))+.06;
text(2-.16,asy_ex,'***','FontSize',48);

% line([xin(1,1),xin(1,2)],[asy_in-.02,asy_in-.02],'color','k')
line([xex(1,1),xex(1,2)],[asy_ex-.02,asy_ex-.02],'color','k','LineWidth',1)

ylim([0 max(max(gexSWsbs(:,1:2)))+.2])
xlim([0.5 2.5])
yticks([0 .2 .4 .6 .8 1]);
ylabel('slope')
xticks([1 2]);
xticklabels({'inter - PwLm','extra - PwLw'});
legend([bip bim],'predicted','measured','location','northwest')
axis square
grid on;
box on; 
set(gca,'LineWidth',1)
set(gca,'FontSize',20)
%% transfer score plots
figure(4);
t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
txts = ["perfect transfer","no transfer"];
sfacec = [0.1 0.3 0.9;0.6 0.4 0.22;0.4 0.2 0.11];

hold on;
for kl =1:max(rtSC(:,2))
    y1p = rtSC(rtSC(:,2)== kl,1);%score
    x1p=repmat(kl,size(y1p,1),1); %col2: measured
    vp(jj) = violin(y1p,'x',[kl],'facecolor',[.5 .5 .5],...
        'edgecolor','k','plotlegend',0,'plotmean',1,'plotmedian',0);    
    vs(jj) = scatter(x1p(:),y1p(:),40,sfacec(3,:),...
        'filled','MarkerFaceAlpha',0.9,'jitter','on','jitterAmount',0.15);
    erp = errorbar(kl,median(y1p),...
        quantile(y1p,0.25)-median(y1p),quantile(y1p,0.75)-median(y1p));
%     erp.Color = [0.35 0.35 0.35];
%     erp = errorbar(kl,mean(y1p),...
%         std(y1p)/(sqrt(size(y1p,1))));
    erp.Color = [0.1 0.1 0.1];
    erp.LineStyle = 'none';
    erp.LineWidth = 2;
    %     vs(jj) = scatter(vx1(:),vy1(:),12,...
    %                     sfacec(3,:),'filled','MarkerFaceAlpha',0.8,...
    %                     'jitter','on','jitterAmount',0.1);
    %     us(jj) = scatter(ux1(:),uy1(:),12,...
    %                     sfacec(1,:),'filled','MarkerFaceAlpha',0.9,...
    %                     'jitter','on','jitterAmount',0);
end

yline(0,'--','LineWidth',2)
yline(1,'--','LineWidth',2)
% asy = max(max(rtSC(:,1)))+.1;
% text(1.5,asy,'*','FontSize',40);
% line([1.1,1.97],[asy-.2,asy-.2],'color','k','LineWidth',2)

% txt11 = text(2.4,1.17,txts(1),'FontSize',18);
% txt21 = text(2.4,-.1,txts(2),'FontSize',18);

xlim([0.5 5.2]);
ylim([-1.5 2.5]);
yticks([0 1]);
xticks([1 2]);
% xticklabels({'inter-PwLm','extra-PwLw'});
set(gca,'xticklabel',{[]})
ylabel('transfer score (a.u.)');

box off;
grid minor;
ax = gca;
ax.FontSize = 32;
ax.LineWidth = 1;