%modelling learning phase data using an exemplar model 
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\MatlabAnalyses')
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Online Data\ExtraInterLee'
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\data_YY_IE22\csv4matlab\'

cd(dir);
%% setting up some basics
phase = 1; %0: LL only 1 learning 2 transfer
pickAllGood = 1; % 0: all trs 1: all good trs 2:new PL combination&&good
%% read csv
T = readtable('IE Learning Task.csv');
load('subjective_variances1_P.mat','sLSd')
load('Learning4Stats_Good.mat','exSubj','inSubj','nExS','nInS','LLSd','Wop');

Subj = [inSubj;exSubj];
[LabelG,~] = ismember(T.Participant,Subj);
IndG= find(LabelG == 1);
gT = T(IndG,:);
nGdSubj = nInS + nExS;
nSubj = nGdSubj;

gT.LLabel = zeros(size(gT,1),1);
gT.LLabel(find(contains(gT.LLType,'Narrow'))) = 1;
gT.LLabel(find(contains(gT.LLType,'Medium'))) = 2;
gT.LLabel(find(contains(gT.LLType,'Wide'))) = 3;

gT.GLabel = zeros(size(gT,1),1); %1: interpolation 2: extrapolation 
GiInd = find(ismember(gT.Participant,inSubj));
GeInd = find(ismember(gT.Participant,exSubj));
gT.GLabel(GiInd) = 1;
gT.GLabel(GeInd) = 2;

clear T IndG 
%%
nGrp = size(unique(gT.GLabel),1); %extrapolation %interpolation
uSubj = unique(gT.Participant);
nLType = size(unique(gT.LLType),1);
nTrial = size(gT,1)/size(uSubj,1);
TakeTy = nLType-1; %becuse if phase == 1
trpt = nTrial/TakeTy; %trials per type
%%
gT.instS = nan(size(gT,1),1);
gT.instS2 = zeros(size(gT,1),1);
%
mu = zeros(size(gT.NetX,1),1);
for ii = 1:size(gT,1)
if  contains(gT.Outlier(ii),'no')
gT.instS(ii) = instantslope(gT.NetX(ii),gT.LLMean(ii),mu(ii)); %too volatile 
end
end
%%
for jj = 1:nSubj %inst: latest 10 trs
    pT = unique(gT.Participant,'stable');%75 
    tT = gT(gT.Participant == pT(jj),:);%200 trs
    indPt = find(gT.Participant == pT(jj));%200 indeces
    for kk = 1:2
        tmpinstS = instantslope3(tT.NetX((kk-1)*100+1:kk*100),tT.LLMean((kk-1)*100+1:kk*100),...
            tT.Outlier((kk-1)*100+1:kk*100));
            gT.instS2(indPt((kk-1)*100+1:kk*100)) = tmpinstS ;
    end
end
%%
N = 5;%fixed N
gT.eN = zeros(size(gT,1),1); %modelled net position%exemplar model
gT.sLSD = zeros(size(gT,1),1);
for jk = 1:nSubj %instS2: latest 10 trs
    pT1 = unique(gT.Participant,'stable');%75
    indPt1 = find(gT.Participant == pT1(jk));%200 indeces
    for ik = 1:size(indPt1)
        gT.sLSD(indPt1(ik)) = sLSd(jk,gT.LLabel(indPt1(ik)));
    end
end
%%
for ji = 1:nSubj %instS3
    pT = unique(gT.Participant,'stable');%75
    indPt = find(gT.Participant == pT(ji));%200 indeces
    tT = gT(gT.Participant == pT(ji),:);%200 trs
    tT1 = sortrows(tT,{'TrialNum'},{'ascend'});%sort forexemplar fitting
    tmpeN = simpexemplar(tT1.LLMean,tT1.CoinX,tT1.sLSD,N);
    tT1.eN = tmpeN;
    tT2= tT1(tT.TrialNum,:);%reverse the trial order
    gT.eN(indPt) = tT2.eN;
end
%%
gTI = gT(gT.GLabel == 1,:);
gTE = gT(gT.GLabel == 2,:);
%%
subj = gT.Participant;
nType=3;
nSj = nSubj;
statData = zeros(nSj*3,7);%1 sub 2 grp 3 LL 4 wts 5 oi 6 oic 7.NTr
for i = 1:nSj
    for k = 1:nType %1: narrow, 2: medium, 3: wide
        if i <= nExS
            statData(k+(i-1)*nType,1) =exSubj(i);
            statData(k+(i-1)*nType,2) =0;
            statData(k+(i-1)*nType,3) =LLSd(k);
            if k == 1 %narrow              
                useTmp = intersect(find(subj == exSubj(i)),find(contains(gT.LLType,'Narrow')));
            elseif k == 2
                useTmp = intersect(find(subj == exSubj(i)),find(contains(gT.LLType,'Medium')));
            else
                if phase == 1
                continue
                elseif phase == 2
                useTmp = intersect(find(subj == exSubj(i)),find(contains(gT.LLType,'Wide'))); 
                end
            end 
        elseif i > nExS
            statData(k+(i-1)*nType,1) = inSubj(i-nExS);
            statData(k+(i-1)*nType,2) = 1;
            statData(k+(i-1)*nType,3) = LLSd(k);
            if k == 1 %narrow
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(gT.LLType,'Narrow')));
            elseif k == 2
                if phase == 1
                continue
                elseif phase == 2
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(gT.LLType,'Medium')));
                end
            else k == 3
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(gT.LLType,'Wide')));
            end
        end

        useTmpP = useTmp(26:end); 
        nanI = isnan(gT.eN(useTmpP));
        useTmpP(nanI == 1) = []; %remove nan trials
        useTr = intersect(useTmpP,find(contains(gT.Outlier,'no')));

        statData(k+(i-1)*nType,8) = size(useTr,1);
        splashX = table2array(gT(useTr,8)); %8 LLMean
        coinX = table2array(gT(useTr,10));
        eN = table2array(gT(useTr,26));%modelled net positions 
        fitdata = polyfit(splashX,eN,1);
        oitmp = table2array(gT(useTr,19));
        oictmp = table2array(gT(useTr,20));       
        statData(k+(i-1)*nType,4) = fitdata(1);
        statData(k+(i-1)*nType,5) = fitdata(2);
        statData(k+(i-1)*nType,6) = median(oitmp);
        statData(k+(i-1)*nType,7) = mean(oictmp);
    end
end
%%
clearvars -except statData
DataTable = array2table(statData,...
    'VariableNames',{'Subjects','groups','Likelihoods','Weights','Intercept','OI','OIc','Tr'});
% save('Exempler4Stats_P.mat');