%simulate exemplar model data, transfer phase 
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\MatlabAnalyses')
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Online Data\PSTrang\'
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\data_SC_PS22\csv4matlab\'

cd(dir);
%% setting up some basics
pickAllGood = 1; % 0: all trs 1: all good trs 2:find Mat isoutlier
T = readtable('SP Transfer Task.csv');
load('subjective_variances1_P.mat','sLSd')

T.TLabel = zeros(size(T,1),1);
T.TLabel(find(contains(T.TType,'pl'))) = 1;
T.TLabel(find(contains(T.TType,'pL'))) = 2;
T.TLabel(find(contains(T.TType,'Pl'))) = 3;
T.TLabel(find(contains(T.TType,'PL'))) = 4;
load("SubjData.mat");

if pickAllGood == 1
    nSj = nGdSubj;
else
    nSj = nSubj;
end

Group = T.Group;
block = T.Block;
type = T.TType;
uGroup = unique(Group);
uType = flip(unique(type)); %pl pL Pl PL
nType = size(uType,1);
nTrial = size(type,1)/nSubj;

Subj = [SSubj;PSubj];
[LabelG,~] = ismember(T.Participant,Subj);
IndG= find(LabelG == 1);
gtT = T(IndG,:);

gtT.PLabel = zeros(size(gtT,1),1);
gtT.PLabel(find(gtT.PriorSD == 0.0250)) = 1;
gtT.PLabel(find(gtT.PriorSD == 0.0850)) = 2;

gtT.LLabel = zeros(size(gtT,1),1);
gtT.LLabel(find(contains(gtT.LLType,'Narrow'))) = 1;
gtT.LLabel(find(contains(gtT.LLType,'Wide'))) = 2;

gtT.GLabel = zeros(size(gtT,1),1); %1: interpolation 2: extrapolation 
GsInd = find(ismember(gtT.Participant,SSubj));
GpInd = find(ismember(gtT.Participant,PSubj));
gtT.GLabel(GsInd) = 1;
gtT.GLabel(GpInd) = 2;

load('LearnGT.mat','gT')
glT = gT;
clear gT
%%
N = 5;%fixed N=5/20
gtT.eN = zeros(size(gtT,1),1); %modelled net position%exemplar model
gtT.sLSD = zeros(size(gtT,1),1);
for jk = 1:nSj 
    pT1 = unique(gtT.Participant,'stable');%75
    indPt1 = find(gtT.Participant == pT1(jk));%200 indeces
    for ik = 1:size(indPt1)
        gtT.sLSD(indPt1(ik)) = sLSd(jk,gtT.LLabel(indPt1(ik)));
    end
end

for ji = 1:nSj %instS2: latest 10 trs
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
%% %%
subj = gtT.Participant;

sz = [nSj*2 11];
varTypes = ["double","double","double","double",...
    "string","double","string","double","string",...
    "string","double"];%"double","double"
varNames = ["beta","OI","OIc","TLabel",...
    "TType","PLabel","PType","LLabel","LType",...
    "group","participant",];%"GLabel","Tr"
plttable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
Intercept = zeros(nSj,nType); 

for i = 1:nSj
    for k = 1:nType %1pl 2pL 3Pl 4PL%not related to Var ABCD1234 
        if i <= nSS
            if k == 1 %pl
                useTmp = intersect(find(subj == SSubj(i)),find(contains(gtT.TType,'pl')));
                tmptp = "pl"; tmppr = "p"; tmpl = "l";
                tmppnum = 1; tmplnum = 1;
            elseif k == 2 %pL
                useTmp = intersect(find(subj == SSubj(i)),find(contains(gtT.TType,'pL')));
                tmptp = "pL"; tmppr = "p";tmpl = "L";
                tmppnum = 1; tmplnum = 2;
            elseif k == 3 %Pl
                useTmp = intersect(find(subj == SSubj(i)),find(contains(gtT.TType,'Pl')));
                tmptp = "Pl"; tmppr = "P"; tmpl = "l";                
                tmppnum = 2; tmplnum = 1;
            elseif k == 4 %PL
                useTmp = intersect(find(subj == SSubj(i)),find(contains(gtT.TType,'PL')));
                tmptp = "PL";tmppr = "P";tmpl = "L";
                tmppnum = 2; tmplnum = 2;
            end

        elseif i > nSS
            if k == 1 %narrow/narrow
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(gtT.TType,'pl')));
                tmptp = "pl"; tmppr = "p"; tmpl = "l";
                tmppnum = 1; tmplnum = 1;
            elseif k == 2 %n/W
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(gtT.TType,'pL')));
                tmptp = "pL"; tmppr = "p";tmpl = "L";
                tmppnum = 1; tmplnum = 2;
            elseif k == 3 %W/n
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(gtT.TType,'Pl')));
                tmptp = "Pl"; tmppr = "P"; tmpl = "l";                
                tmppnum = 2; tmplnum = 1;
            elseif k == 4 %W/W
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(gtT.TType,'PL')));
                tmptp = "PL";tmppr = "P";tmpl = "L";
                tmppnum = 2; tmplnum = 2;
            end
        end

        if isempty(useTmp)
            continue
        end

        if i <= nSS %%
            if any(plttable.participant == SSubj(i),'all') %1row already got data from the currrent subj
                plttable.group(i*2) = "serial"; %group/cond serial parallel
                plttable.participant(i*2) = SSubj(i);
                plttable.TType(i*2) = tmptp;
                plttable.TLabel(i*2) = k;
                plttable.PType(i*2) = tmppr; 
                plttable.PLabel(i*2) = tmppnum;                
                plttable.LType(i*2) = tmpl;
                plttable.LLabel(i*2) = tmplnum ;
                tmpmk = 2;
            else %first of this subj
                plttable.group((i-1)*2+1) = "serial"; %group/cond serial parallel
                plttable.participant((i-1)*2+1) = SSubj(i);
                plttable.TType((i-1)*2+1) = tmptp;
                plttable.TLabel((i-1)*2+1) = k;                
                plttable.PType((i-1)*2+1) = tmppr; 
                plttable.PLabel((i-1)*2+1) = tmppnum;                
                plttable.LType((i-1)*2+1) = tmpl;
                plttable.LLabel((i-1)*2+1) = tmplnum ;
                tmpmk = 1;
            end
        elseif i > nSS
            if any(plttable.participant == PSubj(i-nSS),'all') %1row having data of the currrent subj
                plttable.group(i*2) = "parallel"; %group/cond serial parallel
                plttable.participant(i*2) = PSubj(i-nSS);
                plttable.TType(i*2) = tmptp;
                plttable.TLabel(i*2) = k;                
                plttable.PType(i*2) = tmppr; 
                plttable.PLabel(i*2) = tmppnum;                
                plttable.LType(i*2) = tmpl;
                plttable.LLabel(i*2) = tmplnum ;
                tmpmk = 2;
            else %first time
                plttable.group((i-1)*2+1) = "parallel"; %group/cond serial parallel                
                plttable.participant((i-1)*2+1) = PSubj(i-nSS);
                plttable.TType((i-1)*2+1) = tmptp;
                plttable.TLabel((i-1)*2+1) = k;                
                plttable.PType((i-1)*2+1) = tmppr; 
                plttable.PLabel((i-1)*2+1) = tmppnum;                
                plttable.LType((i-1)*2+1) = tmpl;
                plttable.LLabel((i-1)*2+1) = tmplnum ;
                tmpmk = 1;
            end
        end

        useTmpP1 = intersect(useTmp,find(contains(gtT.Outlier,'no')));%take out outliers
        useTr = setdiff(useTmpP1,find(isnan(gtT.eN)));

        splashX = table2array(gtT(useTr,10)); % slash pos = LLMean(10)
        coinX = table2array(gtT(useTr,12));% coin pos ()
        netX = table2array(gtT(useTr,18));% coin pos ()

        eN = table2array(gtT(useTr,28));
        fitdata = polyfit(splashX,eN,1);%fit for modelled slope 

        oitmp = table2array(gtT(useTr,21));%OI (21)
        oictmp = table2array(gtT(useTr,22));%OIc (22)
        Intercept(i,k) = fitdata(2);

        if tmpmk == 1
            plttable.beta((i-1)*2+1) = fitdata(1); %group/cond serial parallel
            plttable.OI((i-1)*2+1) = median(oitmp);
            plttable.OIc((i-1)*2+1) = median(oictmp);
        elseif tmpmk == 2
            plttable.beta(i*2) = fitdata(1); %group/cond serial parallel
            plttable.OI(i*2) = median(oitmp);
            plttable.OIc(i*2) = median(oictmp);
        end
    end
end

clearvars -except plttable

%save('ExemplerTransfer4Stats.mat');