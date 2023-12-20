%Read Transfer task data for iSW
%parallel/serial data (which has 2 prior 2 LL in learning/transfer)
%%house keeping
%%output ..iSW...mat
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\MatlabAnalyses')
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Online Data\ExtraInterLee'
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\data_YY_IE22\csv4matlab\'

cd(dir);
%% setting up some basics
pickAllGood = 1; % 0: all trs 1: all good trs 2:new PL combination&&good
%% read csv
T = readtable('IE Transfer Task.csv');
load('subjective_variances1_P.mat','sLSd')
load('Learning4Stats_Good.mat','exSubj','inSubj','nExS','nInS','LLSd','Wop');

Subj = [inSubj;exSubj];
[LabelG,~] = ismember(T.Participant,Subj);
IndG= find(LabelG == 1);
gtT = T(IndG,:);
nGdSubj = nInS + nExS;
nSubj = nGdSubj;

gtT.LLabel = zeros(size(gtT,1),1);
gtT.LLabel(find(contains(gtT.LLType,'Narrow'))) = 1;
gtT.LLabel(find(contains(gtT.LLType,'Medium'))) = 2;
gtT.LLabel(find(contains(gtT.LLType,'Wide'))) = 3;

gtT.GLabel = zeros(size(gtT,1),1); %1: interpolation 2: extrapolation 
GiInd = find(ismember(gtT.Participant,inSubj));
GeInd = find(ismember(gtT.Participant,exSubj));
gtT.GLabel(GiInd) = 1;
gtT.GLabel(GeInd) = 2;

clear T IndG
load('LearnGT.mat','gT')
glT = gT;
clear gT
%%
nGrp = size(unique(gtT.GLabel),1); %extrapolation %interpolation
uSubj = unique(gtT.Participant);
nLType = size(unique(gtT.LLType),1);
nTrial = size(gtT,1)/size(uSubj,1);
TakeTy = nLType; %becuse if phase == 1
trpt = nTrial/TakeTy; %trials per type
%%
gtT.eN = zeros(size(gtT,1),1); %modelled net position%exemplar model
gtT.sLSD = zeros(size(gtT,1),1);
for jk = 1:nSubj %instS2: latest 10 trs
    pT1 = unique(gtT.Participant,'stable');%75
    indPt1 = find(gtT.Participant == pT1(jk));%200 indeces
    for ik = 1:size(indPt1)
        gtT.sLSD(indPt1(ik)) = sLSd(jk,gtT.LLabel(indPt1(ik)));
    end
end
%%
N = 5;%fixed N=5 N=20
for ji = 1:nSubj %instS2: latest 10 trs
    pT = unique(gtT.Participant,'stable');%75
    indPt = find(gtT.Participant == pT(ji));%
    tT = gtT(gtT.Participant == pT(ji),:);%
    tT1 = sortrows(tT,{'TrialNumOrig'},{'ascend'});%sort forexemplar fitting
    CoinE = glT.CoinX(glT.Participant == pT(ji));
    LLE = glT.LLMean(glT.Participant == pT(ji));
    tmpeN = simpexemplarT(LLE,CoinE,tT1.LLMean,tT1.sLSD,N);
    tT1.eN = tmpeN;
    tT2= tT1(tT.TrialNumOrig,:);%reverse the trial order
    gtT.eN(indPt) = tT2.eN;
end
%%
gtTI = gtT(gtT.GLabel == 1,:);
gtTE = gtT(gtT.GLabel == 2,:);
%%
phase =2;
subj = gtT.Participant;
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
                useTmp = intersect(find(subj == exSubj(i)),find(contains(gtT.LLType,'Narrow')));
            elseif k == 2
                useTmp = intersect(find(subj == exSubj(i)),find(contains(gtT.LLType,'Medium')));
            else
                if phase == 1
                continue
                elseif phase == 2
                useTmp = intersect(find(subj == exSubj(i)),find(contains(gtT.LLType,'Wide'))); 
                end
            end 
        elseif i > nExS
            statData(k+(i-1)*nType,1) = inSubj(i-nExS);
            statData(k+(i-1)*nType,2) = 1;
            statData(k+(i-1)*nType,3) = LLSd(k);
            if k == 1 %narrow
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(gtT.LLType,'Narrow')));
            elseif k == 2
                if phase == 1
                continue
                elseif phase == 2
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(gtT.LLType,'Medium')));
                end
            else k == 3
                useTmp = intersect(find(subj == inSubj(i-nExS)),find(contains(gtT.LLType,'Wide')));
            end
        end

        useTmpP = useTmp(1:end); %the only change => remove first 50 trials
        nanI = isnan(gtT.eN(useTmpP));
        useTmpP(nanI == 1) = []; %remove nan trials
        useTr = intersect(useTmpP,find(contains(gtT.Outlier,'no')));

        statData(k+(i-1)*nType,8) = size(useTr,1);
        splashX = table2array(gtT(useTr,8)); %8 LLMean
        coinX = table2array(gtT(useTr,10));
        eN = table2array(gtT(useTr,24));
        fitdata = polyfit(splashX,eN,1);
        oitmp = table2array(gtT(useTr,19));
        oictmp = table2array(gtT(useTr,20));       
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
% save('ExemplerTransfer4Stats.mat');