function [chi2] = chi2Val(obsVectors)
%CHI2 Calculates chi2 values  
%Each row points to one observation
Exp=expectedVal(obsVectors);
chi2=sum(((obsVectors-Exp).^2)./Exp,'all');
end

