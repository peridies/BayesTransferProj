%Read Transfer task data
%parallel/serial data (which has 2 prior 2 LL in learning/transfer)
%%house keeping
%%output ..4stats_Good.mat
clear all;
close all;

addpath('C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Code')
dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp1_discovery\'
% dir = 'C:\Users\CSLIN\OneDrive - The University of Melbourne\Documents\Project_BayesTransfer\Data\exp1_validation\'

cd(dir);
%% setting up some basics
phase = 2; %0: LL only 1 learning 2 transfer
pickAllGood = 1; % 0: all trs 1: all good trs 2:find Mat isoutlier
%read csv
T = readtable('SP Transfer Task.csv');
load("SubjData.mat");
%%
% notes: changed Ex to P In to S and LLType to TType
% if pickAllGood ~= 2
nGrp = 2; %serial/parallel
subj = T.Participant;
Group = T.Group;
type = T.TType;

uSubj = unique(subj);
uGroup = unique(Group);
uType = flip(unique(type)); %pl pL Pl PL
nSubj = size(uSubj,1);
nType = size(uType,1);
nTrial = size(subj,1)/size(uSubj,1);

nGdSubj = nSubj-size(badSubj,1);
% get Pr, LL,, and optimal weighting
PrSD = unique(T.PriorSD);
nPr = size(PrSD,1);
LLSD = unique(T.LikelihoodSD);
nLL = size(LLSD,1);
WoT = zeros(nPr,nLL);
%%
if pickAllGood == 1
    nSj = nGdSubj;
else
    nSj = nSubj;
end

%1Wpl 2WpL 3WPl 4WPL 5OIpl 6OIpL 7OIPl 8PL
%9Group 10subj
sz = [nSj*2 13];
varTypes = ["double","double","double","double",...
    "string","double","string","double","string",...
    "string","double","double","double"];
varNames = ["beta","OI","OIc","TLabel",...
    "TType","PLabel","PType","LLabel","LType",...
    "group","GLabel","participant","Tr"];
plttable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

% if phase == 1
for i = 1:nSj
    for k = 1:nType %1pl 2pL 3Pl 4PL%not related to Var ABCD1234 
        if i <= nSS
            if k == 1 %pl
                useTmp = intersect(find(subj == SSubj(i)),find(contains(T.TType,'pl')));
                tmptp = "pl"; tmppr = "p"; tmpl = "l";
                tmppnum = 1; tmplnum = 1;
            elseif k == 2 %pL
                useTmp = intersect(find(subj == SSubj(i)),find(contains(T.TType,'pL')));
                tmptp = "pL"; tmppr = "p";tmpl = "L";
                tmppnum = 1; tmplnum = 2;
            elseif k == 3 %Pl
                useTmp = intersect(find(subj == SSubj(i)),find(contains(T.TType,'Pl')));
                tmptp = "Pl"; tmppr = "P"; tmpl = "l";                
                tmppnum = 2; tmplnum = 1;
            elseif k == 4 %PL
                useTmp = intersect(find(subj == SSubj(i)),find(contains(T.TType,'PL')));
                tmptp = "PL";tmppr = "P";tmpl = "L";
                tmppnum = 2; tmplnum = 2;
            end

        elseif i > nSS
            plotData(i,9) = 2;
            plotData(i,10) = PSubj(i-nSS);
            if k == 1 %narrow/narrow
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(T.TType,'pl')));
                tmptp = "pl"; tmppr = "p"; tmpl = "l";
                tmppnum = 1; tmplnum = 1;
            elseif k == 2 %n/W
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(T.TType,'pL')));
                tmptp = "pL"; tmppr = "p";tmpl = "L";
                tmppnum = 1; tmplnum = 2;
            elseif k == 3 %W/n
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(T.TType,'Pl')));
                tmptp = "Pl"; tmppr = "P"; tmpl = "l";                
                tmppnum = 2; tmplnum = 1;
            elseif k == 4 %W/W
                useTmp = intersect(find(subj == PSubj(i-nSS)),find(contains(T.TType,'PL')));
                tmptp = "PL";tmppr = "P";tmpl = "L";
                tmppnum = 2; tmplnum = 2;
            end
        end

        if isempty(useTmp)
            continue
        end

        if i <= nSS
            if any(plttable.participant == SSubj(i),'all') %1row already got data from the currrent subj
                plttable.group(i*2) = "serial"; %group/cond serial parallel
                plttable.GLabel(i*2) = 1; %group/cond serial parallel
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
                plttable.GLabel((i-1)*2+1) = 1; %group/cond serial parallel
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
                plttable.GLabel(i*2) = 2; %group/cond serial parallel
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
                plttable.GLabel((i-1)*2+1) = 2; %group/cond serial parallel
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

        if pickAllGood == 0
            useTr = useTmp;
        else
            useTr = intersect(useTmp,find(contains(T.Outlier,'no')));%take out outliers
        end
        splashX = table2array(T(useTr,10)); % slash pos = LLMean(10)
        coinX = table2array(T(useTr,12));% coin pos ()
        netX = table2array(T(useTr,18));% coin pos ()

        fitdata = polyfit(splashX,netX,1);%linear fit, single individual
%       [b1,s1] = betaplot(splashX,netX,fitdata)

        oitmp = table2array(T(useTr,21));%OI (21)
        oictmp = table2array(T(useTr,22));%OIc (22)
        %         Intercept(i,k) = fitdata(2);

        if tmpmk == 1
            plttable.beta((i-1)*2+1) = fitdata(1); %group/cond serial parallel
            plttable.OI((i-1)*2+1) = median(oitmp);
            plttable.OIc((i-1)*2+1) = median(oictmp);
            plttable.Tr((i-1)*2+1) = size(useTr,1);
        elseif tmpmk == 2
            plttable.beta(i*2) = fitdata(1); %group/cond serial parallel
            plttable.OI(i*2) = median(oitmp);
            plttable.OIc(i*2) = median(oictmp);
            plttable.Tr(i*2) = size(useTr,1);
        end
    end
end
%% save data
clearvars ans coinX splashX dir fitdata...
    i k type ii jj xx yy useTmp...
    tmpmk tmptp temps sz
%% old code saving data, with 96 (48 serial, 47 parallel)
if pickAllGood == 2
    if phase == 1 %
        outlT = nan(nSj,6); %col1:outlier(0/1)PL1 col1:outlier(0/1)PL2 col3:Grp col4:subj
        sbj = reshape(plttable.participant,2,nSj);
        outlT(:,6)  = sbj(1,:)';
        gp = reshape(plttable.GLabel,2,nSj);
        outlT(:,5)  = gp(1,:)';
    elseif phase == 2
        outlT = nan(nSj,6); %col1:outlier(0/1)PL1 col1:outlier(0/1)PL2 col3:Grp col4:subj
        sbj = reshape(plttable.participant,2,nSj);
        outlT(:,6)  = sbj(1,:)';
        gp = reshape(plttable.GLabel,2,nSj);
        outlT(:,5)  = gp(1,:)';
    end
    for tt =1:1 %1:comb1 2:comb2
        for pp=1:nSj%subj            
            outlT(pp,2+tt) = plttable.TLabel((pp-1)*2+tt);
            tmpi = intersect(find(plttable.GLabel == plttable.GLabel((pp-1)*2+tt)),find(plttable.TLabel == outlT(pp,2+tt)));
            tmp = isoutlier(plttable.beta(tmpi));
            podr = plttable.participant(tmpi);%isub of tmpi
            ii = find(podr == outlT(pp,6));
            outlT(pp,tt) = tmp(ii)       
        end
    end
% [ro,~]= find(outlT(:,1:2));
% sort(ro)
% U = unique(ro);
% orow = U(hist(ro,U)>1)
% outliers = outlT(orow,6)
end
%% save data
%        save('Transfer4Stats_GoodG.mat','plttable')
%%
function [s1,p1] = betaplot(x,y,p)
figure;
hold on; 
s1=scatter(x,y,'k','filled');
x1 = -.75:.01:.75;
y1 = polyval(p,x1);
p1=plot(x1,y1,'LineWidth',4,'Color','#094183')%'#203864'
p1=plot(x1,x1,'k--','LineWidth',0.5)%'#203864'
hold off
axis square
xlabel('centroid of splashes')
ylabel('net position')
xlim([-0.75 0.75])
ylim([-0.75 0.75])
yline(0,'--')
end 