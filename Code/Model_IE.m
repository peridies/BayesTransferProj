%% compare models in the exp 2
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp2_validation\'

cd(dir);
%%
load('subjective_variances1_P.mat'); %subjective variances using  root mean square instead of SD
load('SubjDataIE.mat','exSubj','inSubj','nExS',...
    'nInS');
load('Transfer4Stats_Good.mat','DataTable'); %
RealData = DataTable;
clear DataTable
load('ExemplarTransfer4Stats_P.mat');
ExemTData = DataTable
clear DataTable
%% esxWsbs %extraloation sensory weight side-by-side col1 predicted col2 real data from transfer
exSWsbs = zeros(nExS,13);%col1 H0 predictedB %col2 H1 predicted LL
%col3 H2 predicted direct_fit %col4 H3 predicted exemplar
%col5 measured transfer ex: Lw; in:Lm  
%col6 to compare use Ln(the ll that both had)
%col7 Lm/Lw
%col8 grp %col9 subj 
%col10diff1 %col11%diff2%col12 differ3%col13 differ4
for i = 1:nExS
    exSWsbs(i,4) = ExemTData.Weights(i*3);%exem data Lw (trans)
    exSWsbs(i,5) = RealData.Weights(i*3);%real data Lw (trans)
    exSWsbs(i,6) = RealData.Weights(i*3-2);%real data Ln
    exSWsbs(i,7) = RealData.Weights(i*3-1);%real data Lm
    exSWsbs(i,8) = 2;%grp
    exSWsbs(i,9) = RealData.Subjects(i*3);
end

inSWsbs = zeros(nInS,13);
for j = 1:nInS
    inSWsbs(j,4) = ExemTData.Weights((j+nExS)*3-1);%exem data Lm (trans)
    inSWsbs(j,5) = RealData.Weights((j+nExS)*3-1);%real data Lm (trans)
    inSWsbs(j,6) = RealData.Weights((j+nExS)*3-2);%real data Ln
    inSWsbs(j,7) = RealData.Weights((j+nExS)*3);%real data Lw
    inSWsbs(j,8) = 1;
    inSWsbs(j,9) = RealData.Subjects(j*3);
end
%% compute expected SW
idxe = zeros(nExS,1);
idxi = zeros(nInS,1);
%extra
expWtr = zeros(nExS,1); %expected W%Bayes
exWSlM = zeros(nExS,1); %expected W%L
exWSDM = zeros(nExS,1); %predicted %direct fit model
for  ii = 1:nExS
    lgn = sLGamma(ii,3);
    pgn = msPrGamma_r(ii);
    expWtr(ii) = (lgn*5)/(lgn*5+pgn);
    lG = sLGamma(ii,1:2);%old trial LLgamma
    bo = exSWsbs(ii,5:6);
    exWSlM(ii) = betafit(lG,lgn,bo);
    lgo = lG';
    pgo = repmat(msPrGamma_r(ii),2,1);
    ob = bo';
    exWSDM(ii) = directfit(lgo,pgo,lgn,pgn,ob);
end
exSWsbs(:,1) = expWtr; %col1 predicted slope Bayes
exSWsbs(:,2) = exWSlM; %col2 predicted slope heuristicLL
exSWsbs(:,3) = exWSDM; %col2 predicted slope heuristic

[nani1,~]= find(isnan(exSWsbs(:,1)));%exclude nan
idxe(nani1) = 1;
texSWsbs = exSWsbs(idxe == 0,:); %g = good:not nan
iswO1 = isoutlier(texSWsbs(:,2),'median');%exclude performance outliners
gexSWsbs = texSWsbs; %g = good:median
%% interpolation 
inpWtr = zeros(nInS,1);
inWSlM = zeros(nInS,1);
inWSDM = zeros(nInS,1); %pred%direct fit model

for  jj = 1:nInS
    lgn   = sLGamma(jj+nExS,2);
    pgn = msPrGamma_r(jj+nExS);
    inpWtr(jj) = (lgn*5)/(lgn*5+pgn);
    lG = [sLGamma(jj+nExS,1),sLGamma(jj+nExS,3)];%old trial LLgamma
    bo = exSWsbs(jj,5:6);
    inWSlM(jj) = betafit(lG,lgn,bo);
    lgo = lG';
    pgo = repmat(msPrGamma_r(jj+nExS),2,1);
    ob = bo';
    inWSDM(jj) = directfit(lgo,pgo,lgn,pgn,ob);
end
inSWsbs(:,1) = inpWtr;%col1 predicted Bayes
inSWsbs(:,2) = inWSlM;%col2 predicted heuristic
inSWsbs(:,3) = inWSDM; %col2 predicted slope heuristic
[nani2,~]= find(isnan(inSWsbs(:,1)));
idxi(nani2) = 1;
tinSWsbs = inSWsbs(idxi == 0,:);%not nan
iswO2 = isoutlier(tinSWsbs(:,2),'median');%exclude performance outliners 
ginSWsbs = tinSWsbs; %g = good:median 

% [h1,p1] = ttest(gexSWsbs(:,1),gexSWsbs(:,2));
%%
tSWtable =  [ginSWsbs;gexSWsbs];
isw = isoutlier(tSWtable(:,2),'median');%exclude performance outliners 
SWtable = tSWtable(isw == 0,:); %g = good:median 
SWtable(:,10) = SWtable(:,5)-SWtable(:,1); %abs error data-bayes prediction
SWtable(:,11)= SWtable(:,5)-SWtable(:,2); %abs data-heuristic modelLL
SWtable(:,12)= SWtable(:,5)-SWtable(:,3); %abs data-heuristic model-direct
SWtable(:,13)= SWtable(:,5)-SWtable(:,4); %abs data-exemplar model
%%
removeO = 1;
if removeO
    iswb = isoutlier(SWtable(:,10),'median');%exclude outliners Bayes
    iswl = isoutlier(SWtable(:,11),'median');%exclude outliners LL model
    iboth = [iswb, iswl];
    [row1,col1] = find(iboth)
    SWtable(row1,:) = [];
end
swtest(SWtable(:,10))
swtest(SWtable(:,11))
[hbm,pbm]= ttest(SWtable(:,10),SWtable(:,11))
%%
nSample = size(SWtable,1);
LLb = nSample*log(sum((SWtable(:,5)-SWtable(:,1)).^2)/nSample) %Bayes LL
LLl = nSample*log(sum((SWtable(:,5)-SWtable(:,2)).^2)/nSample) %heuristic LL
LLd = nSample*log(sum((SWtable(:,5)-SWtable(:,3)).^2)/nSample) %heuristic Dire
LLm = nSample*log(sum((SWtable(:,5)-SWtable(:,4)).^2)/nSample) %exemplar

BICb = LLb
BICl = 2*log(nSample) + LLl
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

Titles = ["interpolation","extrapolation"];
xpm = [1.1 1.9 2.7 3.5];
for kk  = 1:2 %group inter vs extra
    if kk  == 1%1 %interpolation
        uidx = find(SWtable(:,8) == 1);
    elseif kk == 2% 2 %extrapolation
        uidx = find(SWtable(:,8) == 2);
    end
    use = SWtable(uidx,:);

    subplot(1,2,kk)
    hold on;
    y1v = use(:,10);
    y2v = use(:,12);%direct
    y3v = use(:,11);%linear model
    y4v = use(:,13);%linear model
    v1 = violin(y1v,'x',[xpm(1)],'facecolor',c_3,...
        'edgecolor',c_3,'plotlegend',0);%predicted
    erv1 = errorbar(xpm(1),median(y1v),...
        quantile(y1v,0.25)-median(y1v),quantile(y1v,0.75)-median(y1v));
    v2 = violin(y2v,'x',[xpm(2)],'facecolor',c_4,...
        'edgecolor',c_3,'plotlegend',0);%predicted
    erv2 = errorbar(xpm(2),median(y2v),...
        quantile(y2v,0.25)-median(y2v),quantile(y2v,0.75)-median(y2v));
    v3 = violin(y3v,'x',[xpm(3)],'facecolor',c_1,...
        'edgecolor',c_3,'plotlegend',0);%predicted
    erv3 = errorbar(xpm(3),median(y3v),...
        quantile(y3v,0.25)-median(y3v),quantile(y3v,0.75)-median(y3v));
    v4 = violin(y4v,'x',[xpm(4)],'facecolor',c_5,...
        'edgecolor',c_3,'plotlegend',0);%predicted
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
    xticks([1.1+.2 1.9+.2 2.7+.2 3.5+.2]);
    xticklabels({'Bayes','linear','likelihood','exemplar'});

    x1v=repmat(xpm(1),size(y1v,1),1);
    s1 = scatter(x1v(:),y1v(:),[],sfacec(2,:),...
        'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0.02);
    x2v=repmat(xpm(2),size(y2v,1),1);
    s2 = scatter(x2v(:),y2v(:),[],sfacec(2,:),...
        'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0.02);
    x3v=repmat(xpm(3),size(y3v,1),1);
    s3 = scatter(x3v(:),y3v(:),[],sfacec(2,:),...
        'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0.02);
    x4v=repmat(xpm(4),size(y4v,1),1);
    s4 = scatter(x4v(:),y4v(:),[],sfacec(2,:),...
        'filled','MarkerFaceAlpha',0.5','jitter','on','jitterAmount',0.02);

    hold on;

    xlim([0.5 4.1]);
%     ylim([-0.2 0.7]);
    ylim([-0.5 0.5]);

    yline(0)

    if kk == 2
        %         legend([bsv2 bsm],{'Predicted','Measured'})%'Location','northeastoutside'
        set(gca,'yticklabels',[]);
    else
        ylabel('data - model prediction');
    end

    title(Titles(kk));
    box on;
    axis square;
%     grid on;
    ax = gca;
    ax.FontSize = 22;
end
%% plot BIC results
BIC = [BICb,BICd,BICl BICm];
dBIC = BIC-(min(BIC))
figure; 
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m, n, p, [0.0 0.01], [0.15 0.01], [0.2 0.2]);
if ~make_it_tight,  clear subplot;  end
subplot(1,2,1);hold on; 
bB = bar(1,BICb,0.3,'FaceColor',c_3)%,'EdgeColor',[1 1 1])%extrapolation
bd = bar(2,BICd,0.3,'FaceColor',c_4)%,'EdgeColor',[1 1 1])%interpolation 
bl = bar(3,BICl,0.3,'FaceColor',c_1)%,'EdgeColor',[1 1 1])
bm = bar(4,BICm,0.3,'FaceColor',c_5,'FaceAlpha',0.5)%,'EdgeColor',[1 1 1])

xlim([0.5 4.5]);
% ylim([-250 114])
% ylim([-290 132])
ylim([-600 273])
% ylim([-660 301])
xticks([1 2 3 4]);
xticklabels({'Bayes','linear','likelihood','exemplar'});
ax = gca;
ax.FontSize = 18;
% axis square;
grid on; 
grid minor;