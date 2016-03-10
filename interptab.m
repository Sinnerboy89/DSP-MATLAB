function [M] = interptab(N,Q)
%INTERPTAB Summary of this function goes here
%   Detailed explanation goes here
M = zeros(Q,N); % M is the lookup table/matrix required for interpolation

% below are 'manual' constructions of the matrix M for specific cases N=3
% and N=4. These were used to compare against an attempt of constructing M
% generally for any positive integer N.

% M2 = M;
% if N == 4;
%     for i=1:Q
%         a = -((Q/2)+i-1)/Q;
%         M2(i,:) = [(a+0.5)*(a-0.5)*(a-1.5)/-6 (a+1.5)*(a-0.5)*(a-1.5)/2 (a+0.5)*(a-1.5)*(a+1.5)/-2 (a+0.5)*(a-0.5)*(a+1.5)/6];
%     end
% elseif N == 3;
%     for i=1:Q
%         a = -((Q/2)+i-1)/Q;
%         M2(i,:) = [a*(a-1)/2 (a+1)*(a-1)/-1 a*(a+1)/2];
%     end
% end

r = (0.5-(N/2):1:-0.5+(N/2)); % contains all possible roots of each component polynomial in ascending order
r2 = r; % required for resetting r within loop

c = zeros(1,N); % this will become our vector of 'corrections', to ensure polynomial=1 at non-root
for i=1:floor(N/2)
    c(i) = factorial(N-i)*factorial(i-1); % we only need to calculate first half of c and use inherent symmetry
end
if mod(N,2) ~= 0
    c(ceil(N/2)) = (factorial((N-1)/2))^2; % calculate 'middle' correction value for N odd
end
c(ceil(N/2)+1:end) = c(floor(N/2):-1:1); % symmetric property of c utilized here
if mod(N,2) == 0
    c(1:2:end) = -c(1:2:end); % signs alternate for each ascending element of c...
else
    c(2:2:end) = -c(2:2:end); % ...sign of first element is dependent on whether N is even or odd
end

for i=1:Q
    a = -((Q/2)+i-1)/Q; % set of allowed alpha values
    for j=1:N
        r2(N-j+1)=[]; % component polynomials are of (N-1) order, hence r needs the appropriate factor removed for each j loop
        M(i,j) = prod(a*ones(1,N-1)+r2)/c(j); % component polynomials calculated for each alpha value
        r2=r; % reset r vector before calculation of next polynomial
    end
end

end
