function [wgpa] = WGPA(A,W,lambda)
%WGPA Weighted Generalized Power Average operator
%   Detailed explanation goes here
    if lambda == 0, lambda = 0.000000001; end
    T=zeros(length(A),1);
    
    for a=1:length(A)
        
        S=0;
        for j=1:length(A)
            if j==a, continue, end
            support = 0.5*exp(-2*(A(a)-A(j))^2);
            %support = abs(log(A(i)) - log(A(j)));
            S=S+W(j)*support;
        end
        T(a)=S;
    end
    
    over=0;
    below=0;
    for a=1:length(A)
        over=over+W(a)*(1+T(a))*A(a)^lambda;
        below=below+W(a)*(1+T(a));
    end
    
    wgpa=(over/below)^(1/lambda);
end

