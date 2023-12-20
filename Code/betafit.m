function slope = betafit(OGamma,NGamma,beta)
P = polyfit(OGamma,beta,1);
slope = polyval(P,NGamma);
end