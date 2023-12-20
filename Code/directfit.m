function slope = directfit(OLGamma,OPGamma,NLGamma,NPGamma,Obeta)
%only care about a map between previous N trials centroid-coin
%find n using learning data
%force b(1) (intercept) = 0
x1 = OLGamma;
x2 = OPGamma;    
x1fit = NLGamma;
x2fit = NPGamma;    
y = Obeta;
%Compute the regression coefficients for a linear model with no interaction term.
X = [ones(size(x1)) x1 x2];
b = regress(y,X);    % Removes NaN data
slope = b(1) + b(2)*x1fit + b(3)*x2fit;
end