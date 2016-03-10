clear all

%This script is a memory-optimized version of
%pv_ChristopherBuchanan_s0809473.m such that the I/O signals are now the
%largest files (i.e. all other datasets previously larger have been reduced, and workspace cleaned up).

tic;

%|----------PREAMBLE (2-3)----------|
N = 1024; % Block size
HA = N/4; % Analysis hop size 
Q = 1.25; % Time stretch factor (>1 slows it down)
HS = ceil(Q*HA); % Synthesis hop size (both hops need to be integers)
Q_i = HS/HA; % Redefines exact stretch factor after integer correction

%|---------HANN WINDOWING FUNCTION (4)----------|

w = 0.5*(1-cos(2*pi*((0:N-1)')/N)); % function taken from lecture 7

%|----------AUDIO READ-IN (5)----------|
[x, SR]=audioread('TestGuitarPhrasephase_stereo.wav'); % reads in stereo .wav and sample rate
x = (x(:,1)+x(:,2))/2; % brutally converts to mono
[~,~,x] = find(x); % chuck out the zeros at start and end
xs = length(x); % total number of samples before padding

%|----------ZERO-PADDING (5)----------|
r = mod((xs-N),HA); % compute remainder when chopping x into OL frames
x = [x;zeros((N-r),1)]; % ensures 'next' frame would be 0s

%|----------NUMBER OF ANALYSIS FRAMES (5)----------|

xsp = length(x); % number of samples after padding
NF = ((xsp-N)/HA)+1; % number of analysis frames

%|----------PHASE VOCODER LOOP INITIALIZATION----------|
X = zeros(N,2);
X(:,1) = x(1:N).*w;
D = zeros(N,2);
v = zeros(N,1);
v(1:((N/2))) = 2*pi*HA*((0:(N/2)-1)')/N; 
v(end:-1:(N/2)+2) = -v(2:((N/2)));
YFP = zeros(N,2);
ys = ((NF-1)*HS)+N; %number of samples in synthesis output y
y = zeros(ys,1);
for i=0:NF-1
    if i~=0
         X(:,1) = X(:,2);
    end
    if i~=NF-1
        X(:,2) = x((1+((i+1)*HA)):(N+((i+1)*HA))).*w;
    end
    if i==NF-1
        X(:,2) = 0;
    end
    
    %----------DFT (7-9)----------|
    XF = fft(X); % DFT of X
    XFM = abs(XF); % magnitudes of X
    XFP = angle(XF); % complex phases of X
    XFP(1,:) = 0;
    XFP(513,:) = 0;
    if i==0
            YFP(:,1) = XFP(:,1); % set first frame of output DFT equal to input DFT
    end
    
    %|----------PHASE DIFFERENCES (10 b)----------|
    D(:,1) = D(:,2);
    D(:,2) = XFP(:,2) - XFP(:,1);
    
    %|----------PHASE INCREMENT (10 c)----------|
    D = D-v(:, ones(1,2)); % subtract v from every column of D
    
    %|----------WRAP TO -/+PI (10 d)----------|
    D(2:end,:) = D(2:end,:) - (2*pi*floor((D(2:end,:)+pi)/(2*pi))); % standard functions only.
    D(1,:) = 0; % first row needn't be touched with clumsy numerical hands - it should contain only 0s at this point
    
    %|----------PHASE RESCALING (10 e)----------|
    D = D + v(:, ones(1,2)); % add back v to every column of D
    D = D*Q_i; % 'stretching' the phase component by our time factor
    
    %|----------PHASE DIFFERENCES -> ABS PHASE (10 f)----------|
    if i~=0
        YFP(:,1) = YFP(:,2);
    end
    YFP(:,2) = YFP(:,1) + D(:,2);
    
    %|----------IDFT (11-12)----------|
    YF = XFM.*exp(YFP*1i); % combining magnitude and phase components via exponential form
    Y = ifft(YF); % back to time domain
    
    %|----------SYNTHESIS PHASE (13)----------|
    y((1+((i)*HS)):(N+((i)*HS))) = y((1+((i)*HS)):(N+((i)*HS))) + Y(:,2).*w; % overlap-add procedure
end
yr = real(y); % numerical inaccuracies cause small imaginary components to remain
toc;

%|----------DEBUG - REAL/IMAG RATIO----------|

% yr_l1 = sum((real(y)).^2);
% yi_l1 = sum((imag(y)).^2);
% disp('L1 norm of reals / L1 norm of imags is...')
% disp(yi_l1/yr_l1)
% c1 = abs(imag(y))./abs(real(y));
% cmax = max(c1);
% c1 = c1/cmax;
% disp('the maximum Im/Re is...')
% disp(cmax)

%|----------DEBUG - PHASE VECTOR ANTI-SYMMETRY CHECK----------|

% vcheck=zeros(1,(N/2)+1);
% for i=2:(N/2)+1
%     if v(i)==-v(N+2-i);
%         vcheck(i)=1;
%     end
% end
% if sum(vcheck)==511
%     disp('phase vector v: good')
% else
%     disp('something is wrong with v.')
%     disp('The following entries do not have an anti-symmetric partner:')
%     disp(find(vcheck<1))
% end

%|----------PLOTS----------|

% hold on
% plot(x, 'k') %after clipping
% plot(xpadd, 'r') %after padding
% plot(ytut7) %after simple time-stretch
% plot(yr) %after PV
% plot(X(:,5),'k') %specific frame at analysis phase
% plot(Y(:,5)) %specific frame after PV
% plot (abs(real(y)))
% plot (abs(imag(y)))
% plot(v,'k')
% plot(c)

%|----------MEMORY USAGE----------|

memory
whos % All matrices have been significantly reduced in size, x or y should now be the largest dataset.




