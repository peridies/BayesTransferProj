%Read Bayes Transfer task data
%starting with Lee's data (which has 1 prior,two/three LL in learning/transfer)
%%house keeping
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_validation\'

cd(dir);
%% 1 read csv %set up some basics
%load learning
load('Learning4Stats_Good_P.mat','exSubj','inSubj','nExS',...
    'nInS','badSubj','nSubj','nGdSubj','nType',...
    'nC ond','DataTable','LLSd','Wop','PrSD');
LStatTable = DataTable; %groups 0 = extra; 1 = inter
LStatTable.LLabel = zeros(size(LStatTable,1),1);
LStatTable.LLabel(LStatTable.Likelihoods == 0.0240) = 1;
LStatTable.LLabel(LStatTable.Likelihoods == 0.0600) = 2;
LStatTable.LLabel(LStatTable.Likelihoods == 0.1500) = 3;
LStatTable.Gamma = zeros(size(LStatTable,1),1);
LStatTable.Gamma(LStatTable.Likelihoods == 0.0240) = 1/(0.0240^2);
LStatTable.Gamma(LStatTable.Likelihoods == 0.0600) = 1/(0.0600^2);
LStatTable.Gamma(LStatTable.Likelihoods == 0.1500) = 1/(0.1500^2);
LStatTable.Phase = zeros(size(LStatTable,1),1);
LearnTable = LStatTable;
LStatTable(~LStatTable.Weights,:) = [];
clear DataTable;
[InEL,~] = ismember(LearnTable.Subjects,exSubj);
[InIL,~] = ismember(LearnTable.Subjects,inSubj);
LearnEx = LearnTable(find(InEL == 1),:);
LearnIn = LearnTable(find(InIL == 1),:);
%load transfer
load("Transfer4stats_Good.mat",'DataTable');
TransTable = DataTable;
TransTable.LLabel(TransTable.Likelihoods == 0.0240) = 1;
TransTable.LLabel(TransTable.Likelihoods == 0.0600) = 2;
TransTable.LLabel(TransTable.Likelihoods == 0.1500) = 3;
% TransTable.Gamma = zeros(225,1);
TransTable.Gamma(TransTable.Likelihoods == 0.0240) = 1/(0.0240^2);
TransTable.Gamma(TransTable.Likelihoods == 0.0600) = 1/(0.0600^2);
TransTable.Gamma(TransTable.Likelihoods == 0.1500) = 1/(0.1500^2);
TransTable.Phase = ones(size(TransTable,1),1);
TStatTable = TransTable;
newExI = intersect(find(TStatTable.groups == 0), find(TStatTable.LLabel == 3))%new:extraPwLw
newInI = intersect(find(TStatTable.groups == 1), find(TStatTable.LLabel == 2))%new:interPwLm
newI = [newExI;newInI];%index of trasnfer-new trials
TStatTable(newI,:) = [];%remove transfer-new from TStatTable
LTOldtable = [LStatTable;TStatTable]; %merge %later for stats
clear DataTable;
[InET,~] = ismember(TransTable.Subjects,exSubj);
[InIT,~] = ismember(TransTable.Subjects,inSubj);
TransEx = TransTable(find(InET == 1),:);
TransIn = TransTable(find(InIT == 1),:);
nSj = nGdSubj
%% statistics
swtest(LStatTable.Weights(:))
newLbeta = (LStatTable.Weights-mean(LStatTable.Weights))/std(LStatTable.Weights);
[hks,pks,~,~] = kstest(newLbeta)%discovery validation  h=1;  

%Learn
l0C = fitlme(LStatTable,'Weights~LLabel+(1|Subjects)')
lgC= fitlme(LStatTable,'Weights~LLabel+groups+(1|Subjects)')
lgiC= fitlme(LStatTable,'Weights~LLabel*groups+(1|Subjects)')

% l0 = fitlme(LStatTable,'Weights~Likelihoods+(1|Subjects)')
% lg = fitlme(LStatTable,'Weights~Likelihoods+groups+(1|Subjects)')
% lgi = fitlme(LStatTable,'Weights~Likelihoods*groups+(1|Subjects)')
% 
% lN0 = fitlme(LStatTable,'Weights~Gamma+(1|Subjects)')
% lNg = fitlme(LStatTable,'Weights~Gamma+groups+(1|Subjects)')
% lNgi = fitlme(LStatTable,'Weights~Gamma*groups+(1|Subjects)')

LL0gC = compare(l0C,lgiC,'CheckNesting',true,'NSim',100) 
% LL0gC = compare(l0C,lgC,'CheckNesting',true,'NSim',100) 
LLgiC = compare(lgC,lgiC,'CheckNesting',true,'NSim',100)%discovery<0.001
% LL0g = compare(l0,lg,'CheckNesting',true,'NSim',100) 
% LLgi = compare(lg,lgi,'CheckNesting',true,'NSim',100) 
% LLogN = compare(lN0,lNg,'CheckNesting',true,'NSim',100)
% LLgiN = compare(lNg,lNgi,'CheckNesting',true,'NSim',100) 
%%
LSEx = LStatTable(LStatTable.groups == 0,:);
LSIn = LStatTable(LStatTable.groups == 1,:);
% le = fitlme(LSEx,'Weights~Likelihoods+(1|Subjects)')
% li = fitlme(LSIn,'Weights~Likelihoods+(1|Subjects)')
leC = fitlme(LSEx,'Weights~LLabel+(1|Subjects)')
liC = fitlme(LSIn,'Weights~LLabel+(1|Subjects)')
mLSEx1 = median(LSEx.Weights(LSEx.LLabel == 1)) 
mLSEx2 = median(LSEx.Weights(LSEx.LLabel == 2)) 
mLSIn1 = median(LSIn.Weights(LSIn.LLabel == 1)) 
mLSIn3 = median(LSIn.Weights(LSIn.LLabel == 3))
iqrLSEx1 = iqr(LSEx.Weights(LSEx.LLabel == 1))
iqrLSEx2 = iqr(LSEx.Weights(LSEx.LLabel == 2))
iqrLSIn1 = iqr(LSIn.Weights(LSIn.LLabel == 1)) 
iqrLSIn3 = iqr(LSIn.Weights(LSIn.LLabel == 3))

[pLE,~] = signrank(LSEx.Weights(LSEx.LLabel == 1),LSEx.Weights(LSEx.LLabel == 2))
[pLI,~] = signrank(LSIn.Weights(LSIn.LLabel == 1),LSIn.Weights(LSIn.LLabel == 3))
%%
LToldtable = [LStatTable;TStatTable];
LToEx = LToldtable(LToldtable.groups == 0,:);
LToIn = LToldtable(LToldtable.groups == 1,:);
%LToldtable
lo0 = fitlme(LToldtable,'Weights~LLabel+(1|Subjects)')
% lop= fitlme(LToldtable,'Weights~LLabel+Phase+(1|Subjects)')
log = fitlme(LToldtable,'Weights~LLabel+groups+(1|Subjects)')
% lopi= fitlme(LToldtable,'Weights~LLabel*Phase+(1|Subjects)')
logi = fitlme(LToldtable,'Weights~LLabel*groups+(1|Subjects)')
loai = fitlme(LToldtable,'Weights~LLabel*groups*Phase+(1|Subjects)')%not better than gi
% LLop = compare(lo0,lop,'CheckNesting',true,'NSim',100) %.43
% LLog = compare(lo0,log,'CheckNesting',true,'NSim',100) %0.019
% LLopi = compare(lo0,lopi,'CheckNesting',true,'NSim',100) %.33
LLog = compare(lo0,log,'CheckNesting',true,'NSim',100) %0.0099
LLggi = compare(log,logi,'CheckNesting',true,'NSim',100) %disco0.11
% LLoai = compare(lo0,loai,'CheckNesting',true,'NSim',100) %0.0099
LLagi = compare(logi,loai,'CheckNesting',true,'NSim',100) %discovery%0.65 %validation.35
LLgai = compare(log,loai,'CheckNesting',true,'NSim',100) %discovery%.38
% LLggi = compare(log,logi,'CheckNesting',true) %discovery0.096 %validation 3.25e-06
%discovery L+transfer old trials: best 

mLToEx1 = median(LToEx.Weights(LToEx.LLabel == 1)) 
mLToEx2 = median(LToEx.Weights(LToEx.LLabel == 2)) 
mLToIn1 = median(LToIn.Weights(LToIn.LLabel == 1)) 
mLToIn3 = median(LToIn.Weights(LToIn.LLabel == 3))
iLToEx1 = iqr(LToEx.Weights(LToEx.LLabel == 1))
iLToEx2 = iqr(LToEx.Weights(LToEx.LLabel == 2))
iLToIn1 = iqr(LToIn.Weights(LToIn.LLabel == 1)) 
iLToIn3 = iqr(LToIn.Weights(LToIn.LLabel == 3))
[pltoE,~] = signrank(LToEx.Weights(LToEx.LLabel == 1),LToEx.Weights(LToEx.LLabel == 2))
[pltoI,~] = signrank(LToIn.Weights(LToIn.LLabel == 1),LToIn.Weights(LToIn.LLabel == 3))

[pLT,~] =signrank(LStatTable.Weights,TStatTable.Weights)
median(LStatTable.Weights)
median(TStatTable.Weights)
iqr(LStatTable.Weights)
iqr(TStatTable.Weights)

writetable(LToldtable,'ie_oldtable.csv','Delimiter',',','QuoteStrings',true)
%% posthoc for the interaction
TSEx = TStatTable(TStatTable.groups == 0,:);
TSIn = TStatTable(TStatTable.groups == 1,:);
te = fitlme(TSEx,'Weights~Likelihoods+(1|Subjects)')
ti = fitlme(TSIn,'Weights~Likelihoods+(1|Subjects)')
letC = fitlme(TSEx,'Weights~LLabel+(1|Subjects)')
litC = fitlme(LSIn,'Weights~LLabel+(1|Subjects)')
mTSEx1 = median(TSEx.Weights(TSEx.LLabel == 1)) 
mTSEx2 = median(TSEx.Weights(TSEx.LLabel == 2)) 
mTSIn1 = median(TSIn.Weights(TSIn.LLabel == 1)) 
mTSIn3 = median(TSIn.Weights(TSIn.LLabel == 3))
iqrTSEx1 = iqr(TSEx.Weights(TSEx.LLabel == 1))
iqrTSEx2 = iqr(TSEx.Weights(TSEx.LLabel == 2))
iqrTSIn1 = iqr(TSIn.Weights(TSIn.LLabel == 1)) 
iqrTSIn3 = iqr(TSIn.Weights(TSIn.LLabel == 3))

[pTE,~] = signrank(TSEx.Weights(TSEx.LLabel == 1),TSEx.Weights(TSEx.LLabel == 2))
[pTI,~] = signrank(TSIn.Weights(TSIn.LLabel == 1),TSIn.Weights(TSIn.LLabel == 3))
%%
lt0 = fitlme(TransTable,'Weights~LLabel+(1|Subjects)')
ltg = fitlme(TransTable,'Weights~LLabel+groups+(1|Subjects)')
lti = fitlme(TransTable,'Weights~LLabel*groups+(1|Subjects)')
LLtg = compare(lt0,ltg,'CheckNesting',true,'NSim',100) %discovery%.1%validationLR15.55 p<.01
LLti = compare(ltg,lti,'CheckNesting',true,'NSim',100) %validation.67
LLt0i = compare(lt0,lti,'CheckNesting',true,'NSim',100) %discovery%.089%validation
LLtgF = compare(lt0,ltg,'CheckNesting',true) %discovery%.1%validationLR15.55 p<.01

median(TransTable.Weights(TransTable.groups == 0)) %ex
median(TransTable.Weights(TransTable.groups == 1)) %in
iqr(TransTable.Weights(TransTable.groups == 0)) %ex
iqr(TransTable.Weights(TransTable.groups == 1)) %in
[pTIE,~] = ranksum(TransTable.Weights(TransTable.groups == 1),TransTable.Weights(TransTable.groups == 0))
%%
TAllEx = TransTable(TransTable.groups == 0,:);
TAllIn = TransTable(TransTable.groups == 1,:);
% lte0 = fitlme(TAllEx,'Weights~Likelihoods+(1|Subjects)')
% lti0 = fitlme(TAllIn,'Weights~Likelihoods+(1|Subjects)')
lteC = fitlme(TAllEx,'Weights~LLabel+(1|Subjects)')
ltiC = fitlme(TAllIn,'Weights~LLabel+(1|Subjects)')
mTAEx1 = median(TAllEx.Weights(TAllEx.LLabel == 1)) 
mTAEx2 = median(TAllEx.Weights(TAllEx.LLabel == 2)) 
mTAEx3 = median(TAllEx.Weights(TAllEx.LLabel == 3)) 
mTAIn1 = median(TAllIn.Weights(TAllIn.LLabel == 1))
mTAIn2 = median(TAllIn.Weights(TAllIn.LLabel == 2))
mTAIn3 = median(TAllIn.Weights(TAllIn.LLabel == 3))
iqrTAEx1 = iqr(TAllEx.Weights(TAllEx.LLabel == 1))
iqrTAEx2 = iqr(TAllEx.Weights(TAllEx.LLabel == 2))
iqrTAEx3 = iqr(TAllEx.Weights(TAllEx.LLabel == 3))
iqrTAIn1 = iqr(TAllIn.Weights(TAllIn.LLabel == 1)) 
iqrTAIn2 = iqr(TAllIn.Weights(TAllIn.LLabel == 2))
iqrTAIn3 = iqr(TAllIn.Weights(TAllIn.LLabel == 3))

[pTE,~] = signrank(TSEx.Weights(TSEx.LLabel == 1),TSEx.Weights(TSEx.LLabel == 2))
[pTI,~] = signrank(TSIn.Weights(TSIn.LLabel == 1),TSIn.Weights(TSIn.LLabel == 3))
%% compare to optimal
pvalul = zeros(4,1);
medianBl = zeros(4,1); %intern -> w -> extran -> medium 
for mm =1:2 %mm1: inter
    for nn = 1:2
        if mm == 1 && nn == 2%mm1 inter, 13 mm2:extra 12
            pickb = LSIn.Weights(LSIn.LLabel == nn+1);
            testb = pickb-Wop(nn+1);
        elseif mm == 1 && nn == 1
            pickb = LSIn.Weights(LSIn.LLabel == nn);
            testb = pickb-Wop(nn);
        else
           pickb = LSEx.Weights(LSEx.LLabel == nn);
           testb = pickb-Wop(nn);
        end
    medianBl((mm-1)*2+nn) = median(pickb);
    pvalul((mm-1)*2+nn) = signrank(testb,0,'tail','right')
    end
end
%% 3.box/violin plot of learn/trans-old/trans-new
figure; %Weights
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m, n, p, [0.01 0.01], [0.15 0.1], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end
myhexvalues = ['#F42272';'#226CE0';'#7030A0']; %PwLn Lm PwLw
myrgb = hex2rgb(myhexvalues)
% Titles = ["interpolation","extrapolation"];
xt = ["learning","transfer"];
xpos = [2 3 4;
    6 7 8]; %234 learning %678 transfer

sfacec = [0.1 0.3 0.9;0.6 0.4 0.22;0.4 0.2 0.11];
xmin = 0.9; xmax = 9.1;

for ij  = 1:2%size(Titles,2)
    subplot(1,2,ij)
    for ph = 1:size(xt,2)
        if ij == 1 && ph == 1 %inter, learning
            use = LearnIn;
        elseif ij == 1 && ph == 2%inter, transfer
            use = TransIn;
        elseif ij == 2 && ph == 1 %extra, learning
            use = LearnEx;
        elseif ij == 2 && ph == 2 %%extra, transfer
            use = TransEx;
        end
        hold on;
        for jj = 1:nType %scatter plot
            if ph == 1
            xx = [xmin xmax];
            yy = [Wop(jj) Wop(jj)];
            hl(jj) = plot(xx,yy,'--','Color',[myrgb(jj,:) .5],'LineWidth',3);
            end
            y1 = use.Weights(use.LLabel == jj);
            idxol = isoutlier(y1,'median');%exclude outliners 
            vy1 = y1(idxol == 0);
            uy1 = y1(idxol == 1);%outliers
            vx1=  repmat(xpos(ph,jj),size(vy1,1),1);
            ux1=  repmat(xpos(ph,jj),size(uy1,1),1);
            if mean(y1) ~= 0
                vs(jj) = scatter(vx1(:),vy1(:),18,...
                    sfacec(3,:),'filled','MarkerFaceAlpha',0.9,...
                    'jitter','on','jitterAmount',0.3);
                us(jj) = scatter(ux1(:),uy1(:),18,...
                    sfacec(1,:),'filled','MarkerFaceAlpha',0.9,...
                    'jitter','on','jitterAmount',0);%outliers
                %             b(jj) = boxchart(x1,y1,'BoxWidth',.75,'LineWidth',1.5);
                %             b(jj).BoxFaceColor = metroclr(jj,:);
                vp(jj) = violin(vy1,'x',[xpos(ph,jj)],'facecolor',myrgb(jj,:),...
                        'edgecolor','k','plotlegend',0,'plotmean',0,'plotmedian',1);
                erp = errorbar(xpos(ph,jj),median(vy1),...
                    quantile(vy1,0.25)-median(vy1),quantile(vy1,0.75)-median(vy1));
                erp.Color = [0.35 0.35 0.35];
                erp.LineStyle = 'none';
                erp.LineWidth = 1;
            end
        end
    end

    box off;
    grid on;
    ax =gca;
    ax.FontSize = 24;
    ax.LineWidth = 1;

    xlim([xmin xmax]);
    ylim([0 1.23]);

    if ij == 1 %interpolation
        xticks([2 4 6 7 8]);
        %xticklabels({'PwLn','PwLw','PwLn','PwLm','PwLw'});
        set(gca,'Xticklabel',[])
        xtickangle(45)
        yticks([0 .2 .4 .6 .8 1]);
        ylabel('slope');
%         txt1 = text(xpos(1,2)-1,1.18,xt(1),'FontSize',24);
    elseif ij  == 2
        xticks([2 3 6 7 8]);
        %         xticklabels({'PwLn','PwLm','PwLn','PwLm','PwLw'});
        set(gca,'Xticklabel',[])
        xtickangle(45)
        set(gca,'Yticklabel',[])
 %       txt1 = text(xpos(1,2)-1.3,1.18,xt(1),'FontSize',24);
        %         oh=legend([vp(1) vp(3) hl(3)],'learning/old condition',...
        %             'new condition','Optimal sw',...
        %             'Location','southeast');
    end
%    txt2 = text(xpos(2,2)-0.9,1.18,xt(2),'FontSize',24);
%     tl = title(Titles(ij));
%     tl.FontSize = 28;
    xline(5)
end