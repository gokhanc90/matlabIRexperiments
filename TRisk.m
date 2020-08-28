function [tRisk] = TRisk(base,run,alpha)
    n=numel(base);
    diff = run - base;
    diff(diff<0)=(alpha+1) * diff(diff<0);
    meanDiff = mean(diff);
    
    sum1 = sum((diff - meanDiff) .* (diff - meanDiff) );
    sum2 = sum(diff - meanDiff);
    
    varianceDifference = (sum1 - (sum2 * sum2 / n)) / (n - 1);
    tRisk = meanDiff / sqrt(varianceDifference / n);
end

