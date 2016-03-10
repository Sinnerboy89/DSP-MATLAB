function [ y ] = feedforcombex_s0809473_Buchanan_Christopher(x,b,M)
%FEEDFORCOMB This is a function to implement P feed forward
% comb filters with P coefficients in vector b and P corresponding
% order vector M entries on an input x in parallel (adding the outputs)

N = length(x);
P = length(M);
y = P*x; % until lowest order filter can be switched on, output will just be gained input

Ms = sort(M); % sort M entries into ascending order
sortorder = zeros(1,P);
for n=1:length(M);
    sortorder(n) = find(M==Ms(n)); % to ensure correct b-M matchup, we need to track what the sorting did
end

for k=1:P
    for i=Ms(k)+1:Ms(P)
        y(i) = y(i)+b(sortorder(k))*x(i-Ms(k)); % this loops deals with the initialization period of all filters (as we move through x in time, more filters will begin contributing)
    end
end

for p=1:P
    y(Ms(P)+1:N) = y(Ms(P)+1:N)+b(sortorder(p))*x(Ms(P)+1-Ms(p):N-Ms(p)); % this loop deals with the 'main' block, where all filters are feeding forward
end

