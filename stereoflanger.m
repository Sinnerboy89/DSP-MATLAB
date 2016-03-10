%|----------PREAMBLE----------|

M0 = 50; % flange depth (larger values turn flanger into a more chorus-like effect)
f0 = 2; % flange frequency
g = 0.7; % flange strength
Q = 100; % Lagrangian interpolation density
N = 5; % Lagrangian interpolation order

%|----------AUDIO READ-IN----------|

[x, Fs]=audioread('CF_clip_2.wav');
if size(x,2) == 1
    x = repmat(x,2); % this script works with stereo input; if mono is given, it'll be copied into L and R channels here
end
Ts = 1/Fs;
L = length(x);
f = Fs*(0:(L/2))/L;

%|----------LAGRANGIAN INTERPOLATION TABLE GENERATION----------|

M = interptab_Buchanan_s0809473(N,Q);
LGMl = zeros(Q,N);
LGMr = LGMl;

%|-----------STEREO LAGRANGIAN FEEDFORWARD COMB FILTER----------|

yffsl = x(:,1);
for i=2*M0+ceil(N/2)+1:L-floor(N/2) % optimization needed, reduce order of interpolation at edges
    DxIl = i-M0*(1+sin(2*pi*f0*i*Ts)); % exact 'delayed' x location
    fracl = DxIl-round(DxIl); % 'displacement' to nearest sample
    DxIQl = round(((floor(N/2)/(N-1))+(fracl/(N-1)))*Q); % nearest Q to DxI
    if mod(N,2) == 0
        for j=1:N
            LGMl(:,j) = x(floor(DxIl)-(N/2)+j,1)*M(:,j);
        end
    else
        for j=1:N
            LGMl(:,j) = x(round(DxIl)-floor(N/2)+j-1,1)*M(:,j);
        end
    end
    LGl = sum(LGMl,2);
    yffsl(i) = x(i,1) + LGl(DxIQl);
end
yffsr = x(:,2);
for i=2*M0+ceil(N/2)+1:L-floor(N/2) % optimization needed, reduce order of interpolation at edges
    DxIr = i-M0*(1+cos(2*pi*f0*i*Ts)); % exact 'delayed' x location
    fracr = DxIr-round(DxIr); % 'displacement' to nearest sample
    DxIQr = round(((floor(N/2)/(N-1))+(fracr/(N-1)))*Q); % nearest Q to DxI
    if mod(N,2) == 0
        for j=1:N
            LGMr(:,j) = x(floor(DxIr)-(N/2)+j,2)*M(:,j);
        end
    else
        for j=1:N
            LGMr(:,j) = x(round(DxIr)-floor(N/2)+j-1,2)*M(:,j);
        end
    end
    LGr = sum(LGMr,2);
    yffsr(i) = x(i,2) + LGr(DxIQr);
end

%|-----------STEREO LAGRANGIAN FEEDBACK COMB FILTER----------|

yfbsl = x(:,1);
for i=2*M0+ceil(N/2)+1:L-floor(N/2) % optimization needed, reduce order of interpolation at edges
    DxIl = i-M0*(1+sin(2*pi*f0*i*Ts)); % exact 'delayed' x location
    fracl = DxIl-round(DxIl); % 'displacement' to nearest sample
    DxIQl = round(((floor(N/2)/(N-1))+(fracl/(N-1)))*Q); % nearest Q to DxI
    if mod(N,2) == 0
        for j=1:N
            LGMl(:,j) = g*yfbsl(floor(DxIl)-(N/2)+j)*M(:,j);
        end
    else
        for j=1:N
            LGMl(:,j) = g*yfbsl(round(DxIl)-floor(N/2)+j-1)*M(:,j);
        end
    end
    LGl = sum(LGMl,2);
    yfbsl(i) = x(i,1) - LGl(DxIQl);
end
yfbsr = x(:,2);
for i=2*M0+ceil(N/2)+1:L-floor(N/2) % optimization needed, reduce order of interpolation at edges
    DxIr = i-M0*(1+cos(2*pi*f0*i*Ts)); % exact 'delayed' x location
    fracr = DxIr-round(DxIr); % 'displacement' to nearest sample
    DxIQr = round(((floor(N/2)/(N-1))+(fracr/(N-1)))*Q); % nearest Q to DxI
    if mod(N,2) == 0
        for j=1:N
            LGMr(:,j) = g*yfbsr(floor(DxIr)-(N/2)+j)*M(:,j);
        end
    else
        for j=1:N
            LGMr(:,j) = g*yfbsr(round(DxIr)-floor(N/2)+j-1)*M(:,j);
        end
    end
    LGr = sum(LGMr,2);
    yfbsr(i) = x(i,2) - LGr(DxIQr);
end

%|----------STEREO COMBINATION----------|

yffs = zeros(L,2);
yffs(:,1) = yffsl; % L channel of feedforward flanger
yffs(:,2) = yffsr; % R channel of feedforward flanger
yfbs = zeros(L,2);
yfbs(:,1) = yfbsl; % L channel of feedback flanger
yfbs(:,2) = yfbsr; % R channel of feedback flanger

%|----------DEBUG - DxIQ CHECK----------|

% if DxIQ > (N/2*(N-1))*Q
%     disp('DxI goes higher than centre-range at:')
%     disp(i)
% elseif DxIQ < ((N-2)/2*(N-1))*Q
%     disp('DxI goes lower than centre-range at:')
%     disp(i)
% end