%simulate exemplar model data, learning phase 
clear all;
close all;

addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\MatlabAnalyses')
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Online Data\PSTrang\'
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\data_SC_PS22\csv4matlab\'

cd(dir);
%% setting up some basics
phase = 1; %0: LL only 1 learning 2 transfer
pickAllGood = 1; % 0: all trs 1: all good trs 2:find Mat isoutlier
%read csv, get table ready for instant slope
if phase == 1
    T = readtable('SP Learning Task.csv');
elseif phase == 2
    T = readtable('SP Transfer Task.csv');
end
T.TLabel = zeros(size(T,1),1);
T.TLabel(find(contains(T.TType,'pl'))) = 1;
T.TLabel(find(contains(T.TType,'pL'))) = 2;
T.TLabel(find(contains(T.TType,'Pl'))) = 3;
T.TLabel(find(contains(T.TType,'PL'))) = 4;
load("SubjData.mat");
load('subjective_variances1_P.mat','sLSd') %will be used in the exemplar model

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
gT = T(IndG,:);

gT.PLabel = zeros(size(gT,1),1);
gT.PLabel(find(gT.PriorSD == 0.0250)) = 1;
gT.PLabel(find(gT.PriorSD == 0.0850)) = 2;

gT.LLabel = zeros(size(gT,1),1);
gT.LLabel(find(contains(gT.LLType,'Narrow'))) = 1;
gT.LLabel(find(contains(gT.LLType,'Wide'))) = 2;

gT.GLabel = zeros(size(gT,1),1); %1: interpolation 2: extrapolation 
GsInd = find(ismember(gT.Participant,SSubj));
GpInd = find(ismember(gT.Participant,PSubj));
gT.GLabel(GsInd) = 1;
gT.GLabel(GpInd) = 2;
%%
N = 5;%fixed N = 5; N = 20
gT.eN = zeros(size(gT,1),1); %modelled net position%exemplar model
gT.sLSD = zeros(size(gT,1),1);
for jk = 1:nSj 
    pT1 = unique(gT.Participant,'stable');%75
    indPt1 = find(gT.Participant == pT1(jk));%200 indeces
    for ik = 1:size(indPt1)
        gT.sLSD(indPt1(ik)) = sLSd(jk,gT.LLabel(indPt1(ik)));
    end
end
%%
for ji = 1:nSj %
    pT = unique(gT.Participant,'stable');%
    indPt = find(gT.Participant == pT(ji));%200 indeces
    pl = unique(gT.PLabel(indPt),'stable');%there are 2
    for pp = 1:2
    indPp = indPt(((pp-1)*200+1):pp*200); 
    tT = gT(indPp,:);%200 trs
    tmpeN = simpexemplar(tT.LLMean,tT.CoinX,tT.sLSD,N);
    tT.eN = tmpeN;
    gT.eN(indPp) = tT.eN;
    end
end
%% compute modelled slope values 
subj = gT.Participant;
sz = [nSj*2 11];
varTypes = ["double","double","double","double",...
    "string","double","string","double","string",...
    "string","double"];%"double","double"
varNames = ["beta","OI","OIc","TLabel",...
    "TType","PLabel","PType","LLabel","LType",...
    "group","participant",];%"GLabel","Tr"
plttable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
Intercept = zeros(nSj,nType); 
% if phase == 1
for i = 1:nSj
    for k = 1:nType %1pl 2pL 3Pl 4PL%not related to Var ABCD1234 
        if i <= nSS
            if k == 1 %pl
                useTmp = intersect(find(subj == SSubj(i)),find(contains(gT.TType,'pl')));
                tmptp = "pl"; tmppr = "p"; tmpl = "l";
                tmppnum = 1; tmplnum = 1;
            elseif k == 2 %pL
                useTmp = intersect(find(subj == SSubj(i)),find(contains(gT.TType,'pL')));
                tmptp = "pL"; tmppr = "p";tmpl = "L";
                tmppnum = 1; tmplnum = 2;
            elseif k == 3 %Pl
                useTmp = intersect(find(subj == SSubj(i)),find(contains(gT.TType,'Pl')));
                tmptp = "Pl"; tmppr = "P"; tmpl = "l";                
                tmppnum = 2; tmplnum = 1;
            elseif k == 4 %PL
                useTmp = intersect(find(subj == SSubj(i)),find(contains(gT.TType,'PL')));
                tmptp = "PL";tmppr = "P";tmpl = "L";
                tmppnum = 2; tmplnum = 2;
            end

        elseif i > nSS
            if k == 1 %narrow/narrow
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(gT.TType,'pl')));
                tmptp = "pl"; tmppr = "p"; tmpl = "l";
                tmppnum = 1; tmplnum = 1;
            elseif k == 2 %n/W
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(gT.TType,'pL')));
                tmptp = "pL"; tmppr = "p";tmpl = "L";
                tmppnum = 1; tmplnum = 2;
            elseif k == 3 %W/n
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(gT.TType,'Pl')));
                tmptp = "Pl"; tmppr = "P"; tmpl = "l";                
                tmppnum = 2; tmplnum = 1;
            elseif k == 4 %W/W
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(gT.TType,'PL')));
                tmptp = "PL";tmppr = "P";tmpl = "L";
                tmppnum = 2; tmplnum = 2;
            end
        end

        if isempty(useTmp)
            continue
        end
        useTmpP = useTmp(51:end); 

        if i <= nSS %%
            if any(plttable.participant == SSubj(i),'all') %1row already got data from the currrent subj
                plttable.group(i*2) = "serial"; %group/cond serial parallel
%                 plttable.GLabel(i*2) = 1; %group/cond serial parallel
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
%                 plttable.GLabel((i-1)*2+1) = 1; %group/cond serial parallel
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
%                 plttable.GLabel(i*2) = 2; %group/cond serial parallel
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
%                 plttable.GLabel((i-1)*2+1) = 2; %group/cond serial parallel
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

        useTmpP1 = intersect(useTmpP,find(contains(gT.Outlier,'no')));%take out outliers
        useTr = setdiff(useTmpP1,find(isnan(gT.eN)));

        splashX = table2array(gT(useTr,10)); % slash pos = LLMean(10)
        coinX = table2array(gT(useTr,12));% coin pos ()
        netX = table2array(gT(useTr,18));

        eN = table2array(gT(useTr,28));% modelled net pos ()
        fitdata = polyfit(splashX,eN,1);% 

        oitmp = table2array(gT(useTr,21));%OI (21)
        oictmp = table2array(gT(useTr,22));%OIc (22)
        Intercept(i,k) = fitdata(2);

        if tmpmk == 1
            plttable.beta((i-1)*2+1) = fitdata(1); %modelled beta
            plttable.OI((i-1)*2+1) = median(oitmp);
            plttable.OIc((i-1)*2+1) = median(oictmp);
        elseif tmpmk == 2
            plttable.beta(i*2) = fitdata(1); 
            plttable.OI(i*2) = median(oitmp);
            plttable.OIc(i*2) = median(oictmp);
        end
    end
end

clearvars -except plttable
% save('Exempler4Stats_P.mat');