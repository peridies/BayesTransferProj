%get modelled slopes=> create oldtable newtable qtable   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SophieLin 
clear all;
close all;
addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')
addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\MatlabAnalyses')
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Online Data\PSTrang\'
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\data_SC_PS22\csv4matlab\'

cd(dir);
%% 1 read csv
    load("Exempler4Stats_P.mat");
    plttable.pNum = zeros(size(plttable,1),1);
    plttable.lNum = zeros(size(plttable,1),1); 
    plttable.pNum(plttable.PLabel == 1) = .025; %narrow
    plttable.pNum(plttable.PLabel == 2) = .085; %wide
    plttable.lNum(plttable.LLabel == 1) = .06; %narrow
    plttable.lNum(plttable.LLabel == 2) = .15; %wide
    LearnTable = plttable;     
    LearnSerial = plttable(contains(plttable.group,'serial'),:);
    LearnParallel= plttable(contains(plttable.group,'parallel'),:);
    clearvars plttable 

    load("ExemplerTransfer4Stats.mat");
    plttable.pNum = zeros(size(plttable,1),1);
    plttable.lNum = zeros(size(plttable,1),1); 
    plttable.pNum(plttable.PLabel == 1) = .025; %narrow
    plttable.pNum(plttable.PLabel == 2) = .085; %wide
    plttable.lNum(plttable.LLabel == 1) = .06; %narrow
    plttable.lNum(plttable.LLabel == 2) = .15; %wide 
    TransTable = plttable;
    TransSerial = plttable(contains(plttable.group,'serial'),:);
    TransParallel= plttable(contains(plttable.group,'parallel'),:);
    clearvars plttable 
    load('SubjData.mat')
    nSj = nGdSubj;

%% 2. new table: only novelly introduced conditions in the transfer 
%old table: learnt trials, compared learning and trasnfer phases 
szn = [nSj*1 8+6];
szo = [nSj*2 10+6];

varTypesn = ["double","double","double","double","string",...
    "double","string","double","string",...
    "double","string","double","double","double"];
varNamesn = ["beta","OI","OIc","TLabel","TType",...
    "PrLabel","PrType","LLabel","LType",...
    "GLabel","group","pNum","lNum","participant"];
varTypeso = ["double","double","double","double","string",...
    "double","string","double","string",...
    "double","string","double","string",...
    "double","double","double"];
varNameso = ["beta","OI","OIc","TLabel","TType",...
    "PrLabel","PrType","LLabel","LType",...
    "GLabel","group","PhLabel","phase",...
    "pNum","lNum","participant"];
newtable = table('Size',szn,'VariableTypes',varTypesn,'VariableNames',varNamesn);%new trs in transfer 
oldtable = table('Size',szo,'VariableTypes',varTypeso,'VariableNames',varNameso);%experienced 

szq = [nSj*1 10+6];
varTypesq = varTypeso; 
varNamesq = varNameso;
qtable = table('Size',szq,'VariableTypes',varTypesq,...
    'VariableNames',varNamesq);%unique


for ii =1:nSj
    if ii <= nSS
        isubj = SSubj(ii);
    elseif  ii > nSS
        isubj = PSubj(ii-nSS);
    end
    tmpLtlabel = LearnTable.TLabel(LearnTable.participant == isubj);
    tmpTtlabel = TransTable.TLabel(TransTable.participant == isubj);
    tmpLtable = LearnTable(LearnTable.participant == isubj,:);
    tmpTtable = TransTable(TransTable.participant == isubj,:);
    newType = setdiff(tmpTtlabel,tmpLtlabel);%identify subj specific new trials
    oldType = intersect(tmpTtlabel,tmpLtlabel);
    qType = setdiff(tmpLtlabel,oldType);
    newtable(ii,1:9) = tmpTtable(tmpTtable.TLabel == newType,1:9);
    newtable(ii,11) = tmpTtable(tmpTtable.TLabel == newType,10);
    newtable(ii,12:13) = tmpTtable(tmpTtable.TLabel == newType,12:13);
    newtable.participant(ii) = isubj;

    oldtable(ii*2-1,1:9) = tmpLtable(tmpLtable.TLabel == oldType,1:9);
    oldtable(ii*2-1,11) = tmpLtable(tmpLtable.TLabel == oldType,10);
    oldtable(ii*2-1,14:15) = tmpLtable(tmpLtable.TLabel == oldType,12:13);
    oldtable(ii*2,1:9) = tmpTtable(tmpTtable.TLabel == oldType,1:9);
    oldtable(ii*2,11) = tmpTtable(tmpTtable.TLabel == oldType,10);
    oldtable(ii*2,14:15) = tmpTtable(tmpTtable.TLabel == oldType,12:13);

    oldtable.phase(ii*2-1) = "Learn";
    oldtable.PhLabel(ii*2-1) = 1;
    oldtable.phase(ii*2) = "Trans";
    oldtable.PhLabel(ii*2) = 2;
    oldtable.participant(ii*2-1:ii*2) = repmat(isubj,2,1);

    qtable(ii,1:9) = tmpLtable(tmpLtable.TLabel == qType,1:9);
    qtable(ii,11) = tmpLtable(tmpLtable.TLabel == qType,10);
    qtable(ii,14:15) = tmpLtable(tmpLtable.TLabel == qType,12:13);
    qtable.phase(ii) = "Learn";
    qtable.PhLabel(ii) = 1;
    qtable.participant(ii) = isubj;
end
oldtable.GLabel(find(contains(oldtable.group,'serial'))) = 1;
oldtable.GLabel(find(contains(oldtable.group,'parallel'))) = 2;
newtable.GLabel(find(contains(newtable.group,'serial'))) = 1;
newtable.GLabel(find(contains(newtable.group,'parallel'))) = 2;

qtable.GLabel(find(contains(qtable.group,'serial'))) = 1;
qtable.GLabel(find(contains(qtable.group,'parallel'))) = 2;

clearvars -except oldtable newtable qtable 
% save('OldNewTable_Exemplar_P.mat','oldtable','newtable','qtable');