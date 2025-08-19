function [coeffs, R2]=regression_fit_values(gt,est)
    
    coeffs = polyfit(gt, est, 1);
    yfit = polyval(coeffs,gt);
    SStot= sum( ( est-mean(est) ).^2 );
    SSres=sum( ( est-yfit ).^2 );
    R2 = 1-SSres/SStot;

end
