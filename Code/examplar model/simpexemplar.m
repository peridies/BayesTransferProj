function net = simpexemplar(LLE,CoinE,LSD,N)
%estimation of net is a weighted sum of exemplars
%weights based on two factors: 
% 1. euclidan distance between exemplar & current splash abs(Splash-LLE)
% 2. splash spreading: LSD (can also use LSD^2)
%change to all based on pdf p(prior sample|LL)
net =zeros(size(LLE,1),1);
for i = 1:size(LLE,1)
    Splash = LLE(i);
    if i == 1
        net(i) = Splash/2;
    elseif i>1 && i<=N %(trial 2-10)%take all trials
        w = pdf('Normal',Splash-LLE(1:i-1),0,LSD(i)/5);%w1 LL based w
        w = w/sum(w);
%         tuning = 1./(LLE(1:i-1)-CoinE(1:i-1)).^2;%dist between LLMean and Coin 
%         w = dist.*tuning/sum(dist.*tuning);
        net(i) = sum(CoinE(1:i-1).*w);
    else
        w = pdf('Normal',Splash-LLE(i-N:i-1),0,LSD(i)/5);%w1 LL based w
        w = w/sum(w);
%         dist = 1./abs(Splash-LLE(i-N:i-1)); %10X1
%         tuning = 1./(LLE(i-N:i-1)-CoinE(i-N:i-1)).^2;%10x1
%         w = dist.*tuning/sum(dist.*tuning);
        net(i) = sum(CoinE(i-N:i-1).*w);
    end
end