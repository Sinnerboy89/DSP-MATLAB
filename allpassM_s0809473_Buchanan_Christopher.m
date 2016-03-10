function [ y ] = allpassM_s0809473_Buchanan_Christopher(x,a,M)
%ALLPASS This is a function to implement an Mth order allpass filter with
% coefficient a on an input x

if size(x,2)==1    

    N = length(x);
    y = zeros(N,1);
    
    g = a;
    
    %initial conditions
    y(1:M) = g*x(1:M);
    y(M+1:end) = g*x(M+1:end)+x(1:end-M);
    
    % loop for feedback component
    for n=M+1:N
        y(n) = y(n) - g*y(n-M);
    end

else
    
    N = length(x);
    y = zeros(N,2);
    
    g = a;
    
    %initial conditions
    y(1:M,1) = g*x(1:M,1);
    y(1:M,2) = g*x(1:M,2);
    y(M+1:end,1) = g*x(M+1:end,1)+x(1:end-M,1);
    y(M+1:end,2) = g*x(M+1:end,2)+x(1:end-M,2);
    
    % loop for left feedback component
    for n=M+1:N
        y(n,1) = y(n,1) - g*y(n-M,1);
    end
    % loop for right feedback component
    for n=M+1:N
        y(n,2) = y(n,2) - g*y(n-M,2);
    end

end

