function net = simpexemplarT(LLE,CoinE,LLMean,sLSD,N)
%LLE: exemplar LLMean
%CoinE: Exemplar coinX
%LLMean: transfer trial LL centre
%sLSD: subject specific LSD  
%estimation of net is a weighted sum of exemplars
%weights based on two factors:
% 1. euclidan distance between exemplar & current splash abs(Splash-LLE)
% 2. splash spreading: LSD (can also use LSD^2)
%change all based on pdf p(prior sample|LL)
s = RandStream('mlfg6331_64');
% Y = datasample(s,1:200,5,'Replace',true)
% datasample(1:200)
net =zeros(size(LLMean,1),1);
for i = 1:size(LLMean,1)
    Splash = LLMean(i);
    [Coin,ind] = datasample(s,CoinE,N); %Y = datasample(s,X,5,2,'Replace',false)
    LL=LLE(ind);
    w = pdf('Normal',Splash-LL,0,sLSD(i)/5);%w1 LL based w
    w = w/sum(w);
    %         dist = 1./abs(Splash-LLE(i-N:i-1)); %10X1
    %         tuning = 1./(LLE(i-N:i-1)-CoinE(i-N:i-1)).^2;%10x1
    %         w = dist.*tuning/sum(dist.*tuning);
    net(i) = sum(Coin.*w);
end
end