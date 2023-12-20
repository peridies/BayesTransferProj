%parallel/serial data (which has 2 prior 2 LL in learning/transfer)
%separate old and new trials in transfer phase for statstics
%%house keeping 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SophieLin 
clear all;
close all;
addpath('C:\Users\CSLIN\Documents\MATLAB\FileExchange')

addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp1_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp1_validation\'

cd(dir);
%% 1 read csv
    load("Learning4stats_GoodG_P.mat");
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

    load("Transfer4stats_GoodG.mat");
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

swtest(LearnTable.beta(:))
newLbeta = (LearnTable.beta-mean(LearnTable.beta))/std(LearnTable.beta);
[hks,pks,~,~] = kstest(newLbeta)
%%
%eg lme = fitlme(shift,'QCDev ~ Shift + (1|Operator)',...
%'FitMethod','REML','DummyVarCoding','effects')
% lme0 = fitlme(LearnTable,'beta~pNum+lNum+(1|participant)')%keep only the best model
% lme0G = fitlme(LearnTable,'beta~pNum+lNum+group+(1|participant)')
% lme0I = fitlme(LearnTable,'beta~pNum*lNum+group+(1|participant)')
% lmeIAll = fitlme(LearnTable,'beta~pNum*lNum*group+(1|participant)')

lce0 = fitlme(LearnTable,'beta~PType+LType+(1|participant)')%keep only the best model
% lce0G = fitlme(LearnTable,'beta~PType+LType+group+(1|participant)')
lcePI = fitlme(LearnTable,'beta~PType*LType+(1|participant)')
lceIAll = fitlme(LearnTable,'beta~PType*LType*group+(1|participant)')

% LL_0G = compare(lme0,lme0G,'CheckNesting',true,'NSim',100)
% LL_0I = compare(lme0,lme0I,'CheckNesting',true,'NSim',100)
% LL_0A = compare(lme0,lmeIAll,'CheckNesting',true,'NSim',100)
LL_0P = compare(lce0,lcePI,'CheckNesting',true,'NSim',100)
LL_0I = compare(lce0,lceIAll,'CheckNesting',true,'NSim',100)
LL_PA = compare(lcePI,lceIAll,'CheckNesting',true,'NSim',100)

betaLS = LearnTable.beta(find(contains(LearnTable.group,'serial')));
betaLP = LearnTable.beta(find(contains(LearnTable.group,'parallel')));
[pL,hL] = ranksum(betaLS,betaLP)
median(betaLS)
median(betaLP)
iqr(betaLS)
iqr(betaLP)
%%
LTableS = LearnTable(find(contains(LearnTable.group,'serial')),:);
% lls0 = fitlme(LTableS,'beta~pNum+lNum+(1|participant)')%keep only the best model
% llscc = fitlme(LTableS,'beta~pNum*lNum+(1|participant)')%keep only the best model
lls0c = fitlme(LTableS,'beta~PType+LType+(1|participant)')%keep only the best model
llsc3 = fitlme(LTableS,'beta~PType*LType+(1|participant)')%keep only the best model
LL_s0c = compare(lls0c,llsc3,'CheckNesting',true)

LTableP = LearnTable(find(contains(LearnTable.group,'parallel')),:);
% llp0 = fitlme(LTableP,'beta~pNum+lNum+(1|participant)')%keep only the best model
% llpcc = fitlme(LTableP,'beta~pNum*lNum+(1|participant)')%keep only the best model
llp0c = fitlme(LTableP,'beta~PType+LType+(1|participant)')%keep only the best model
llpc3 = fitlme(LTableP,'beta~PType*LType+(1|participant)')%keep only the best model

LL_p0c = compare(llp0c,llpc3,'CheckNesting',true,'NSim',100)
%%
median(LTableS.beta(LTableS.PLabel == 1))
median(LTableS.beta(LTableS.PLabel == 2))
median(LTableS.beta(LTableS.LLabel == 1)) %
median(LTableS.beta(LTableS.LLabel == 2))

median(LTableS.beta(LTableS.TLabel == 1))%PnLn
median(LTableS.beta(LTableS.TLabel == 2))%PnLw
median(LTableS.beta(LTableS.TLabel == 3))%PwLn
median(LTableS.beta(LTableS.TLabel == 4))%PwLw
iqr(LTableS.beta(LTableS.TLabel == 1))%PnLn
iqr(LTableS.beta(LTableS.TLabel == 2))%PnLw
iqr(LTableS.beta(LTableS.TLabel == 3))%PwLn
iqr(LTableS.beta(LTableS.TLabel == 4))%PwLw

median(LTableP.beta(LTableP.PLabel == 1))
median(LTableP.beta(LTableP.PLabel == 2))
median(LTableP.beta(LTableP.LLabel == 1)) %M
median(LTableP.beta(LTableP.LLabel == 2))

median(LTableP.beta(LTableP.TLabel == 1))%PnLn
median(LTableP.beta(LTableP.TLabel == 2))%PnLw
median(LTableP.beta(LTableP.TLabel == 3))%PwLn
median(LTableP.beta(LTableP.TLabel == 4))%PwLw
iqr(LTableP.beta(LTableP.TLabel == 1))%PnLn
iqr(LTableP.beta(LTableP.TLabel == 2))%PnLw
iqr(LTableP.beta(LTableP.TLabel == 3))%PwLn
iqr(LTableP.beta(LTableP.TLabel == 4))%PwLw
%%
newOI = LearnTable.OI-1;
[poi,hoi] = signrank(newOI,0,'tail','left');

pvalul = zeros(4,1);
medianBl = zeros(4,1);
for mm =1:4
    pickb = LearnTable.beta(LearnTable.TLabel == mm);
    medianBl(mm) = median(pickb);
    testb = pickb-Wop(mm);
    pvalul(mm) = signrank(testb,0,'tail','right')
end
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
oldtable = table('Size',szo,'VariableTypes',varTypeso,'VariableNames',varNameso);%experienced learn/transfer-old 

szq = [nSj*1 10+6];
varTypesq = varTypeso; 
varNamesq = varNameso;
qtable = table('Size',szq,'VariableTypes',varTypesq,...
    'VariableNames',varNamesq);%unique combinations only presented in the learning phase


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

% save('OldNewTable_Good_P.mat','oldtable','newtable','qtable');
% writetable(oldtable,'ps_oldtable.csv','Delimiter',',','QuoteStrings',true)
% writetable(newtable,'ps_newtable.csv','Delimiter',',','QuoteStrings',true)
%%
ottable = oldtable(oldtable.PhLabel == 2,:);
swtest(oldtable.beta(:))
% lmo0 = fitlme(oldtable,'beta~pNum+lNum+(1|participant)')
lmoc0 = fitlme(oldtable,'beta~PrType+LType+(1|participant)')
lmocI = fitlme(oldtable,'beta~PrType*LType+(1|participant)')%winning
lmpp = fitlme(oldtable,'beta~PrType*LType+phase+(1|participant)')

lmop = fitlme(oldtable,'beta~PrType*LType*phase+(1|participant)')
lmcG = fitlme(oldtable,'beta~PrType+LType+group+(1|participant)')
lmiG = fitlme(oldtable,'beta~PrType*LType*group+(1|participant)')
lmiA = fitlme(oldtable,'beta~PrType*LType*phase*group+(1|participant)')

mobeta1 = median(oldtable.beta(oldtable.TLabel == 1))
mobeta2 = median(oldtable.beta(oldtable.TLabel == 2))
mobeta3 = median(oldtable.beta(oldtable.TLabel == 3))
mobeta4 = median(oldtable.beta(oldtable.TLabel == 4))
qobeta1 = iqr(oldtable.beta(oldtable.TLabel == 1))
qobeta2 = iqr(oldtable.beta(oldtable.TLabel == 2))
qobeta3 = iqr(oldtable.beta(oldtable.TLabel == 3))
qobeta4 = iqr(oldtable.beta(oldtable.TLabel == 4))

%lmop = fitlme(oldtable,'beta~pNum+lNum+phase+(1|participant)')
%lmog = fitlme(oldtable,'beta~pNum+lNum+group+(1|participant)');
%lmopg = fitlme(oldtable,'beta~pNum+lNum+phase+group+(1|participant)');
%lmoI = fitlme(oldtable,'beta~pNum*lNum*phase+group+(1|participant)');
%lmig = fitlme(oldtable,'beta~pNum*lNum*group+(1|participant)');
%lmiA = fitlme(oldtable,'beta~pNum*lNum*phase*group+(1|participant)')

LLoIp = compare(lmocI,lmpp,'CheckNesting',true,'NSim',100) %new vali p .27
LL0c = compare(lmoc0,lmocI,'CheckNesting',true,'NSim',100) %new val.029-.06
LL0g = compare(lmoc0,lmiG,'CheckNesting',true,'NSim',100) %p = .16
%LLo0i =compare(lmo0,lmoI,'CheckNesting',true,'NSim',100) %p = .29
LLca =compare(lmocI,lmiA,'CheckNesting',true,'NSim',100) %new valp = .079
LL0cp = compare(lmocI,lmop,'CheckNesting',true,'NSim',100) % new vali p = .86
LLcig = compare(lmcG,lmiG,'CheckNesting',true,'NSim',100) %BF10: 0.64(R)
LL0cg = compare(lmocI,lmiG,'CheckNesting',true,'NSim',100) % new vali p = .05-.08
LL0gai = compare(lmiG,lmiA,'CheckNesting',true,'NSim',100) % new val p=.188
% LLo0ia = compare(lmocI,lmiA,'CheckNesting',true,'NSim',100) % new val
%posthoc
motbeta1 = median(oldtable.beta(oldtable.TLabel == 1))%PnLn
motbeta2 = median(oldtable.beta(oldtable.TLabel == 2))
motbeta3 = median(oldtable.beta(oldtable.TLabel == 3))
motbeta4 = median(oldtable.beta(oldtable.TLabel == 4))
qotbeta1 = iqr(oldtable.beta(oldtable.TLabel == 1))
qotbeta2 = iqr(oldtable.beta(oldtable.TLabel == 2))
qotbeta3 = iqr(oldtable.beta(oldtable.TLabel == 3))
qotbeta4 = iqr(oldtable.beta(oldtable.TLabel == 4))


betaOL = oldtable.beta(find(contains(oldtable.phase,'Learn')));
betaOT = oldtable.beta(find(contains(oldtable.phase,'Trans')));
[pLT,hLT] = signrank(betaOL,betaOT)

lmotI = fitlme(ottable,'beta~PrType*LType+(1|participant)')
lmot = fitlme(ottable,'beta~PrType+LType+(1|participant)')
LLoT = compare(lmot,lmotI,'CheckNesting',true,'NSim',100) %new vali p .27


motbeta1 = median(ottable.beta(ottable.TLabel == 1))
motbeta2 = median(ottable.beta(ottable.TLabel == 2))
motbeta3 = median(ottable.beta(ottable.TLabel == 3))
motbeta4 = median(ottable.beta(ottable.TLabel == 4))
qotbeta1 = iqr(ottable.beta(ottable.TLabel == 1))
qotbeta2 = iqr(ottable.beta(ottable.TLabel == 2))
qotbeta3 = iqr(ottable.beta(ottable.TLabel == 3))
qotbeta4 = iqr(ottable.beta(ottable.TLabel == 4))
% LL_to0I = compare(lmot0,lmotI,'CheckNesting',true,'NSim',100)

%% create table for statistics of new trials 
szpn = [nSj*1 16];
pre_new = table('Size',szpn,'VariableTypes',varTypeso,'VariableNames',varNameso);%size 95X16
pre_learn = table('Size',szo,'VariableTypes',varTypeso,'VariableNames',varNameso);

pre_new(:,1:11) = newtable(:,1:11); %all transfer "new" trials
pre_new(:,14:16) = newtable(:,12:14);%only these are testing
pre_new.PhLabel = (repmat(2,nSj,1));%real transfer 
pre_new.phase(find(pre_new.PhLabel == 2)) = "Trans";

pre_learn(:,1:9) = LearnTable(:,1:9);%all learn trials
pre_learn(:,11) = LearnTable(:,10); %group
pre_learn(:,16) = LearnTable(:,11);%participant
pre_learn(:,14:15) = LearnTable(:,12:13);
pre_learn.GLabel(find(contains(pre_learn.group,'serial'))) = 1;
pre_learn.GLabel(find(contains(pre_learn.group,'parallel'))) = 2;
pre_learn.PhLabel = (repmat(1,nSj*2,1));
pre_learn.phase(find(pre_learn.PhLabel == 1)) = "Learn";

ontable4s = [pre_learn; pre_new]; %all data 
%%
on_l0 = fitlme(ontable4s,'beta~pNum+lNum+(1|participant)');
on_lp = fitlme(ontable4s,'beta~pNum+lNum+phase+(1|participant)'); %p = .85
%on_lg = fitlme(ontable4s,'beta~pNum+lNum+group+(1|participant)'); %p = .98
%on_lpg = fitlme(ontable4s,'beta~pNum+lNum+phase+group+(1|participant)'); %p = .99
%on_lI = fitlme(ontable4s,'beta~pNum*lNum*phase+group+(1|participant)') %p = .07
on_lip = fitlme(ontable4s,'beta~phase*pNum*lNum+(1|participant)'); %p = .03
on_l0C = fitlme(ontable4s,'beta~PrType*LType+(1|participant)') 
on_lipC = fitlme(ontable4s,'beta~PrType*LType*phase+(1|participant)') 
on_laiC = fitlme(ontable4s,'beta~PrType*LType*phase*group+(1|participant)') 

LLa0ipC = compare(on_l0C,on_lipC,'CheckNesting',true,'NSim',100) %p =.019
LLa0iaC = compare(on_lipC,on_laiC,'CheckNesting',true) %p =.019

%LLa0p = compare(on_l0,on_lp,'CheckNesting',true,'NSim',100) %p = .85
%LLa0g = compare(on_l0,on_lg,'CheckNesting',true,'NSim',100) %p = .98
%LLa0pg = compare(on_l0,on_lpg,'CheckNesting',true,'NSim',100)  %p = .99
%LLa0i =compare(on_l0,on_lI,'CheckNesting',true,'NSim',100) %p = .07
LLa0ip = compare(on_l0,on_lip,'CheckNesting',true,'NSim',100) %p =.03

% fitlme(pre_learn,'beta~pNum+lNum+(1|participant)')
% ln0 = fitlme(pre_new,'beta~pNum+lNum+(1|participant)')
% fitlme(pre_new,'beta~PrLabel+LLabel+(1|participant)')
% lna = fitlme(pre_new,'beta~group*PrType*LType+(1|participant)')
ln0Pc = fitlme(pre_new,'beta~PrType+(1|participant)')
ln0PLc = fitlme(pre_new,'beta~PrType+LType+(1|participant)')

lngC = fitlme(pre_new,'beta~group+PrType+LType+(1|participant)')
lnaC = fitlme(pre_new,'beta~group*PrType*LType+(1|participant)')
LL21C = compare(ln0Pc,ln0PLc,'CheckNesting',true,'NSim',100) %%validation .85
LLn0gC = compare(ln0PLc,lngC,'CheckNesting',true,'NSim',100) %%validation .85
LLn0aC = compare(ln0PLc,lnaC,'CheckNesting',true,'NSim',100) %dicoveryLR=6.48,p=.27 %validation .21
%%
betaTrLn = newtable.beta(find(newtable.LLabel == 1)); %narrow likelihood
betaTrLw = newtable.beta(find(newtable.LLabel == 2));
betaTrPn = newtable.beta(find(newtable.PrLabel == 1)); %narrow pr
betaTrPw = newtable.beta(find(newtable.PrLabel == 2));
[pTrL,~] = ranksum(betaTrLn,betaTrLw, "tail","right")
[pTrP,~] = ranksum(betaTrPn,betaTrPw, "tail","left")
bNS = newtable.beta(find(newtable.LLabel == 1 & newtable.GLabel == 1)); %Ln serial
bWS = newtable.beta(find(newtable.LLabel == 2 & newtable.GLabel == 1)); %Ln serial
bNP = newtable.beta(find(newtable.LLabel == 1 & newtable.GLabel == 2));%Ln serial
bWP = newtable.beta(find(newtable.LLabel == 2 & newtable.GLabel == 2)); %Ln serial
medbNS = median(newtable.beta(find(newtable.LLabel == 1 & newtable.GLabel == 1))) %Ln serial
medbWS =median(newtable.beta(find(newtable.LLabel == 2 & newtable.GLabel == 1))) %Ln serial
medbNP = median(newtable.beta(find(newtable.LLabel == 1 & newtable.GLabel == 2)))%Ln serial
medbWP = median(newtable.beta(find(newtable.LLabel == 2 & newtable.GLabel == 2))) %Ln serial
median(betaTrLn) %validation: .81
median(betaTrLw) %validation: .77
%%
find(contains(oldtable.phase,'Trans'))
pre_old = oldtable(find(contains(oldtable.phase,'Trans')),:);
n1 = 'on' %add old new trials before merge
n2 = 'prgr' %4 groups for plot based on pr/serial vs parallel 
n3 = 'plon' %inside 4 subplot, 6 based on LL and ontrials

pre_learn.(n1) = ones(size(pre_learn,1),1);
pre_old.(n1) = ones(size(pre_old,1),1);
pre_new.(n1) = ones(size(pre_new,1),1)*2;
ontable4p = [pre_learn;pre_old;pre_new];

ontable4p.(n2) = repmat(0,size(ontable4p,1),1);
ontable4p.(n3) = repmat(0,size(ontable4p,1),1);

r1 = find(ontable4p.GLabel == 1 & ontable4p.PrLabel == 1);%serial narrow%pr&grp based 
r2 = find(ontable4p.GLabel == 1 & ontable4p.PrLabel == 2);%serial wide%label
r3 = find(ontable4p.GLabel == 2 & ontable4p.PrLabel == 1);%parallel narrow 
r4 = find(ontable4p.GLabel == 2 & ontable4p.PrLabel == 2);%parallel wide
prgr1 = repmat(1,size(r1,1),1);prgr2 = repmat(2,size(r2,1),1);
prgr3 = repmat(3,size(r3,1),1);prgr4 = repmat(4,size(r4,1),1);

ontable4p.prgr(r1) = prgr1;
ontable4p.prgr(r2) = prgr2;
ontable4p.prgr(r3) = prgr3;
ontable4p.prgr(r4) = prgr4;

rL1 = find(ontable4p.PhLabel == 1 & ontable4p.LLabel == 1);%Learn-Ln
rL2 = find(ontable4p.PhLabel == 1 & ontable4p.LLabel == 2);%Learn-Lw
rL3 = find(ontable4p.PhLabel == 2 & ontable4p.LLabel == 1 & ontable4p.on == 1);%Tr-LnO
rL4 = find(ontable4p.PhLabel == 2 & ontable4p.LLabel == 2 & ontable4p.on == 1);%Tr-LwO
rL5 = find(ontable4p.PhLabel == 2 & ontable4p.LLabel == 1 & ontable4p.on == 2);%Tr-LnN
rL6 = find(ontable4p.PhLabel == 2 & ontable4p.LLabel == 2 & ontable4p.on == 2);%Tr-LwN
pL1 = repmat(1,size(rL1,1),1);pL2 = repmat(2,size(rL2,1),1);
pL3 = repmat(3,size(rL3,1),1);pL4 = repmat(4,size(rL4,1),1);
pL5 = repmat(5,size(rL5,1),1);pL6 = repmat(6,size(rL6,1),1);

ontable4p.plon(rL1) = pL1;
ontable4p.plon(rL2) = pL2;
ontable4p.plon(rL3) = pL3;
ontable4p.plon(rL4) = pL4;
ontable4p.plon(rL5) = pL5;
ontable4p.plon(rL6) = pL6;

nPrGr = max(unique(ontable4p.prgr))
nplon = max(unique(ontable4p.plon))

%% %% 3.box plot of learn/trans-old/trans-new - prior oriented
figure; %Weights
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m, n, p, [0.01 0.01], [0.11 0.08], [0.10 0.10]);
if ~make_it_tight,  clear subplot;  end

myhexvalues = ['#FF6700';'#ADE25D';'#F42272';'#7030A0']; %PnLn PnLw PwLn PwLw
myrgb = hex2rgb(myhexvalues)


% Titles = ["narrow prior","wide prior"];
XT1 = ["learning","transfer-","transfer-"];
XT2 = ["","Old","New"];

YT = ["Serial","Parallel"];
xpos = [2 3 5 6 8 9];
sfacec = [0.9 0.59 0.33;0.4 0.2 0.11];

Nnote = zeros(4,6);

for ij  = 1:4%:4 %subplot number
%1 1 2 3 4 5 6 Ln-L Lw-L Ln-To Lw-To Ln-TN Lw-TN   
u2idx = find(ontable4p.prgr == ij); %subplot by pgrp label
use = ontable4p(u2idx,:);
subplot(2,2,ij)
hold on;

yy1 = [Wop(1) Wop(1)];
yy2 = [Wop(2) Wop(2)];
yy3 = [Wop(3) Wop(3)];
yy4 = [Wop(4) Wop(4)];
xmin = 0.5; xmax = 10.5;
xlim([xmin xmax]);
xx = [xmin xmax];
if mod(ij,2) == 0 %wide prior 
    plot(xx,yy1,'k--','LineWidth',1);
    plot(xx,yy2,'k--','LineWidth',1);
    plot(xx,yy3,'--','Color',[myrgb(3,:) .5],'LineWidth',3);
    plot(xx,yy4,'--','Color',[myrgb(4,:) .5],'LineWidth',3);
    set(gca,'Yticklabel',[])
else
    plot(xx,yy1,'--','Color',[myrgb(1,:) .5],'LineWidth',3);
    plot(xx,yy2,'--','Color',[myrgb(2,:) .5],'LineWidth',3);
    plot(xx,yy3,'k--','LineWidth',1);
    plot(xx,yy4,'k--','LineWidth',1);
    ylabel({YT((ij-1)/2+1);'sensory weight'});
end
hold on;
for jj =1:nplon 
    y1 = use.beta(use.plon == jj);
    x1 = repmat(xpos(jj),size(y1,1),1);
    Nnote(ij,jj) = size(y1,1);
        sc(jj) = scatter(x1(:),use.beta(use.plon == jj),12,sfacec(2,:),...
        'filled','MarkerFaceAlpha',0.85,...
        'jitter','on','jitterAmount',0.25);
    if mod(ij,2) == 1 && mod(jj,2) == 1 % PnLn
        fcvalue = myrgb(1,:);
    elseif mod(ij,2) == 1 && mod(jj,2) == 0 %PnLw
        fcvalue = myrgb(2,:);
    elseif mod(ij,2) == 0 && mod(jj,2) == 1 % PwLn
        fcvalue = myrgb(3,:);
    elseif mod(ij,2) == 0 && mod(jj,2) == 0 %PwLw
        fcvalue = myrgb(4,:);
    end
%     erp = errorbar(xpos(jj),median(y1),...
%         quantile(y1,0.25)-median(y1),quantile(y1,0.75)-median(y1));
%     erp.Color = [0.2 0.2 0.2];
%     erp.LineStyle = 'none';
%     erp.LineWidth = 2;
    vp(jj) = violin(y1,'x',[xpos(jj)],'facecolor',fcvalue,...
        'edgecolor','k','plotlegend',0);%,'BoxWidth',.75,'LineWidth',1.5)
end
%figure setup apply to all 
ylim([-0.2 1.5]);
xticks([2 3 5 6 8 9]);
box on;
% axis square;
grid on;
ax = gca;
ax.FontSize = 16;


% if ij > 2 && mod(ij,2) == 1
% xticklabels({'PnLn','PnLw','PnLn','PnLw','PnLn','PnLw'});
% xtickangle(45)
% elseif ij > 2 && mod(ij,2) == 0
% xticklabels({'PwLn','PwLw','PwLn','PwLw','PwLn','PwLw'});
% xtickangle(45)
% else
set(gca,'Xticklabel',[]) 
% tl = title(Titles(ij));
% tl.FontSize = 20;
% end

if mod(ij,2) == 0 %wide prior 
    set(gca,'Yticklabel',[])
else
    yticks([0 0.5 1]);
    ylabel('slope')
end

txt11 = text(xpos(1)-.2,1.45,XT1(1),'FontSize', 12);
txt21 = text(xpos(3)-.2,1.45,XT1(2),'FontSize', 12);
txt31 = text(xpos(5)-.2,1.45,XT1(3),'FontSize', 12);
txt12 = text(xpos(1),1.35,XT2(1),'FontSize', 12);
txt22 = text(xpos(3)+.2,1.35,XT2(2),'FontSize', 12);
txt32 = text(xpos(5)+.2,1.35,XT2(3),'FontSize', 12);
end
%% new violin plot
%3.box plot of learn/trans-old/trans-new - prior oriented
figure(10); %Weights
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m, n, p, [0.02 0.02], [0.08 0.08], [0.10 0.10]);
if ~make_it_tight,  clear subplot;  end

myhexvalues = ['#FF6700';'#ADE25D';'#F42272';'#7030A0']; %PnLn PnLw PwLn PwLw
myrgb = hex2rgb(myhexvalues)

% Titles = ["narrow prior","wide prior"];
XT1 = ["learning","transfer-","transfer-"];
XT2 = ["","Old","New"];

YT = ["Serial","Parallel"];
xpos = [2 3 5 6];
sfacec = [0.9 0.59 0.33;0.4 0.2 0.11];

Nnote = zeros(4,6);

for ij  = 1:6%N subplot 1:Serial-learn 2 Serial transfer-Pn 3 serial transfer-Pw 
%subplot 1&4 1PnLn 2PnLw 3PwLn 4PwLw 
%subplot 2&5 1PnLn-old 2PnLw-old 3PnLn-new 4 PnLw-new 
%subplot 3&6 1PwLn-old 2PwLw-old 3PwLn-new 4 PwLw-new 
if ij == 1 %1 serial-learn
    u2idx = intersect(find(ontable4p.GLabel == 1),find(ontable4p.PhLabel == 1)); %subplot by pgrp label
elseif ij == 2 %serial-transfer-Pn
    tmpi1 = intersect(find(ontable4p.GLabel == 1),find(ontable4p.PhLabel == 2)); %subplot by pgrp label
    tmpi2 = intersect(find(ontable4p.PhLabel == 2),find(ontable4p.PrLabel == 1)); %subplot by pgrp label
    u2idx = intersect(tmpi1,tmpi2); %subplot by pgrp label
elseif ij == 3 %serial-transfer-Pw
    tmpi1 = intersect(find(ontable4p.GLabel == 1),find(ontable4p.PhLabel == 2)); %subplot by pgrp label
    tmpi2 = intersect(find(ontable4p.PhLabel == 2),find(ontable4p.PrLabel == 2)); %subplot by pgrp label
    u2idx = intersect(tmpi1,tmpi2); %subplot by pgrp label
elseif ij == 4 %1 parallel-learn
    u2idx = intersect(find(ontable4p.GLabel == 2),find(ontable4p.PhLabel == 1)); %subplot by pgrp label
elseif ij == 5 %parallel-transfer-Pn
    tmpi1 = intersect(find(ontable4p.GLabel == 2),find(ontable4p.PhLabel == 2)); %subplot by pgrp label
    tmpi2 = intersect(find(ontable4p.PhLabel == 2),find(ontable4p.PrLabel == 1)); %subplot by pgrp label
    u2idx = intersect(tmpi1,tmpi2); %subplot by pgrp label
elseif ij == 6 %parallel-transfer-Pw
    tmpi1 = intersect(find(ontable4p.GLabel == 2),find(ontable4p.PhLabel == 2)); %subplot by pgrp label
    tmpi2 = intersect(find(ontable4p.PhLabel == 2),find(ontable4p.PrLabel == 2)); %subplot by pgrp label
    u2idx = intersect(tmpi1,tmpi2); %subplot by pgrp label
end
use = ontable4p(u2idx,:);

subplot(2,3,ij)
hold on;

if mod(ij,3) == 0 %trans-wide prior
    plot(xx,yy1,'k--','LineWidth',1);
    plot(xx,yy2,'k--','LineWidth',1);
    plot(xx,yy3,'--','Color',[myrgb(3,:) .5],'LineWidth',3);
    plot(xx,yy4,'--','Color',[myrgb(4,:) .5],'LineWidth',3);
    for jj =1:4
        if jj == 1 %old-PwLn
            uj = intersect(find(use.on == 1),find(use.LLabel == 1));
            fcvalue = myrgb(3,:);
        elseif jj == 2 %old =PwLw
            uj = intersect(find(use.on == 1),find(use.LLabel == 2));
            fcvalue = myrgb(4,:);
        elseif jj == 3
            uj = intersect(find(use.on == 2),find(use.LLabel == 1));
            fcvalue = myrgb(3,:);
        elseif jj == 4
            uj = intersect(find(use.on == 2),find(use.LLabel == 2));
            fcvalue = myrgb(4,:);
        end
        y1 = use.beta(uj);
        x1 = repmat(xpos(jj),size(y1,1),1);
        Nnote(ij,jj) = size(y1,1);
        sc(jj) = scatter(x1(:),y1(:),27,sfacec(2,:),...
            'filled','MarkerFaceAlpha',0.95,...
            'jitter','on','jitterAmount',0.25);
        vp(jj) = violin(y1,'x',[xpos(jj)],'facecolor',fcvalue,...
            'edgecolor','k','plotlegend',0);%,'BoxWidth',.75,'LineWidth',1.5)
        erp = errorbar(xpos(jj),median(y1),...
            quantile(y1,0.25)-median(y1),quantile(y1,0.75)-median(y1));
        erp.Color = [0.1 0.1 0.1];
        erp.LineStyle = 'none';
        erp.LineWidth = 1;
    end

elseif mod(ij,3) == 2 %trans-narrow prior
    plot(xx,yy1,'--','Color',[myrgb(1,:) .5],'LineWidth',3);
    plot(xx,yy2,'--','Color',[myrgb(2,:) .5],'LineWidth',3);
    plot(xx,yy3,'k--','LineWidth',1);
    plot(xx,yy4,'k--','LineWidth',1);
    for jj =1:4
        if jj == 1 %old-PwLn
            uj = intersect(find(use.on == 1),find(use.LLabel == 1));
            fcvalue = myrgb(1,:);
        elseif jj == 2 %old-PwLw
            uj = intersect(find(use.on == 1),find(use.LLabel == 2));
            fcvalue = myrgb(2,:);
        elseif jj == 3 %new-PwLn
            uj = intersect(find(use.on == 2),find(use.LLabel == 1));
            fcvalue = myrgb(1,:);
        elseif jj == 4 %new-PwLw
            uj = intersect(find(use.on == 2),find(use.LLabel == 2));
            fcvalue = myrgb(2,:);
        end
        y1 = use.beta(uj);
        x1 = repmat(xpos(jj),size(y1,1),1);
        Nnote(ij,jj) = size(y1,1);
        sc(jj) = scatter(x1(:),y1(:),27,sfacec(2,:),...
            'filled','MarkerFaceAlpha',0.95,...
            'jitter','on','jitterAmount',0.25);
        vp(jj) = violin(y1,'x',[xpos(jj)],'facecolor',fcvalue,...
            'edgecolor','k','plotlegend',0);
            erp = errorbar(xpos(jj),median(y1),...
                quantile(y1,0.25)-median(y1),quantile(y1,0.75)-median(y1));
            erp.Color = [0.1 0.1 0.1];
            erp.LineStyle = 'none';
            erp.LineWidth = 1;
    end
elseif mod(ij,3) == 1 %learn
    plot(xx,yy1,'--','Color',[myrgb(1,:) .5],'LineWidth',3);
    plot(xx,yy2,'--','Color',[myrgb(2,:) .5],'LineWidth',3);
    plot(xx,yy3,'--','Color',[myrgb(3,:) .5],'LineWidth',3);
    plot(xx,yy4,'--','Color',[myrgb(4,:) .5],'LineWidth',3);
    for jj =1:4
        fcvalue = myrgb(jj,:);
        y1 = use.beta(use.TLabel == jj);
        x1 = repmat(xpos(jj),size(y1,1),1);
        Nnote(ij,jj) = size(y1,1);
        sc(jj) = scatter(x1(:),y1(:),27,sfacec(2,:),...
            'filled','MarkerFaceAlpha',0.95,...
            'jitter','on','jitterAmount',0.25);
        vp(jj) = violin(y1,'x',[xpos(jj)],'facecolor',fcvalue,...
            'edgecolor','k','plotlegend',0);
        erp = errorbar(xpos(jj),median(y1),...
            quantile(y1,0.25)-median(y1),quantile(y1,0.75)-median(y1));
        erp.Color = [0.1 0.1 0.1];
        erp.LineStyle = 'none';
        erp.LineWidth = 1;
    end
end

%figure setup apply to all 
ax = gca;
ax.FontSize = 32;
ax.LineWidth = 1;
xlim([0.5 7.5])
ylim([-0.2 1.5]);
xticks([2 3 5 6]);
box off; %axis square;
grid on;
set(gca,'Xticklabel',[]) 


if mod(ij,3) == 1 %wide prior 
    yticks([0 0.5 1]);
    ylabel('slope')
else
    set(gca,'Yticklabel',[]);
end
end
%% 5. statistics of old trials 
%normality test SWtest null (H0) normality
nType = size(unique(oldtable.TType),1)
nCond = size(unique(oldtable.phase),1)
HvalueBO = nan(nCond,nType);
pvBO = nan(nCond,nType);
HvalueBN = nan(1,nType);
pvBN = nan(1,nType);
for i = 1:nCond
    for j  = 1: nType
        if i == 1
        swind =  intersect(find(oldtable.TLabel == j),find(contains(oldtable.phase,'Learn')))
        else
        swind =  intersect(find(oldtable.TLabel == j),find(contains(oldtable.phase,'Trans')))
        end
            [HvalueBO(i,j), pvBO(i,j)] = swtest(oldtable.beta(swind)) 
        end          
end
for jj  = 1: nType
    swind =  find(newtable.TLabel == j);
    [HvalueBN(1,jj), pvBN(1,jj)] = swtest(newtable.beta(swind))
end
%% pooling all data 
swtest(oldtable.beta(:)) 
swtest(newtable.beta(:)) 

%%fit mixed linear model 
%prepare old_new trans phase data 
tmpOT = oldtable(oldtable.PhLabel == 2,1:16);
tmpOT= removevars(tmpOT,{'PhLabel','phase'}); %transfer phase all old cond

n1 = 'on'
tmpOT.(n1) = ones(95,1);
newtable.(n1) = 2*ones(95,1);
ontable = [tmpOT; newtable];

lmeMO = fitlme(oldtable,'beta~PrLabel+LLabel+GLabel+phase');
lme2IO = fitlme(oldtable,'beta~(PrLabel*GLabel)+(LLabel*GLabel)+(phase*GLabel)');
lmeMN = fitlme(newtable,'beta~PrLabel+LLabel+GLabel');
lme2IN = fitlme(newtable,'beta~PrLabel+LLabel+(PrLabel*GLabel)+(LLabel*GLabel)');

lmeMon = fitlme(ontable,'beta~PrLabel+LLabel+GLabel+on');
lmeMonr = fitlme(ontable,'beta~PrLabel+LLabel+GLabel+on+(1|participant)');
lme2on = fitlme(ontable,'beta~PrLabel+LLabel+(PrLabel*on)+(LLabel*on)+(PrLabel*GLabel)+(LLabel*GLabel)');
lme2onr = fitlme(ontable,'beta~PrLabel+LLabel+(PrLabel*on)+(LLabel*on)+(PrLabel*GLabel)+(LLabel*GLabel)+(1|participant)');
%%
function [s1,p1] = betaplot(x,y,p)
figure;
hold on; 
s1=scatter(x,y,'k',filled');
y1 = polyval(p,x);
p1=plot(x,y1,'LineWidth',2,'Color','#203864')
hold off
axis equal
xlabel('centroid of splashes')
ylabel('net positions')
end 