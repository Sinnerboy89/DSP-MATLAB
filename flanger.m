tic

%|----------PREAMBLE----------|

M0 = 50; % flange depth (larger values turn flanger into a more chorus-like effect)
f0 = 2; % flange frequency
g = 0.7; % flange strength
Q = 100; % Lagrangian interpolation density
N = 5; % Lagrangian interpolation order
N2 = N; % required for any stencil size 'resetting' in Lagrangian flangers

%|----------AUDIO READ-IN----------|

[x, Fs]=audioread('CF_clip_2.wav'); % an exciting glimpse into the murky world of wehuntjackals, solo project of Christopher Buchanan
if size(x,2)>1
    x = x(:,1); % if loading another clip (and I struggle to think why you would), this deals with stereo and mono signals
end
Ts = 1/Fs; % time between adjacent samples
L = length(x);
f = Fs*(0:(L/2))/L; % frequency range of input in Hz - needed for plots

%|-----------VARIABLE (BRUTE-FORCE TRUNCATION) FEEDFORWARD COMB FILTER----------|

yff = x;
for i=1:L
    if (i-round(M0*(1+sin(2*pi*f0*i*Ts)))) >= 1 % this ensures delayed x exists
        yff(i) = x(i) + x(i-round(M0*(1+sin(2*pi*f0*i*Ts)))); % 'squidginess' audible upon playback
    end
end

%|-----------VARIABLE (BRUTE-FORCE TRUNCATION) FEEDBACK COMB FILTER----------|

yfb = x;
for i=1:L
    if (i-round(M0*(1+sin(2*pi*f0*i*Ts)))) >= 1 % this ensures delayed yfb exists
        yfb(i) = x(i) - g*yfb(i-round(M0*(1+sin(2*pi*f0*i*Ts)))); % deepened flange, but still with 'squidginess'
    end
end

%|----------LAGRANGIAN INTERPOLATION TABLE GENERATION----------|

M = interptab_Buchanan_s0809473(N,Q);
M2 = M; % required for resetting
LGM = zeros(Q,N);
LGM2 = LGM; % required for resetting

%|-----------LAGRANGIAN FEEDFORWARD COMB FILTER----------|

yffl = x;
for i=2*M0+ceil(N/2)+1:L-floor(N/2) % i defines position of 'stencil', whose size is governed by N
    DxI = i-M0*(1+sin(2*pi*f0*i*Ts)); % exact 'delayed' x location
    frac = DxI-round(DxI); % 'displacement' to nearest sample
    DxIQ = round(((floor(N/2)/(N-1))+(frac/(N-1)))*Q); % nearest Q value to DxI's relative position within N-stencil
    if mod(N,2) == 0 % dealing with odd and even N cases seperately
        for j=1:N
            LGM(:,j) = x(floor(DxI)-(N/2)+j)*M(:,j); % specific Langrangian polynomials for each point in stencil calculated here
        end
    else
        for j=1:N
            LGM(:,j) = x(round(DxI)-floor(N/2)+j-1)*M(:,j); % specific Langrangian polynomials for each point in stencil calculated here
        end
    end
    LG = sum(LGM,2); % polynomial which passes through N samples within stencil obtained here
    yffl(i) = x(i) + LG(DxIQ); % read off polynomial at specific 'alpha' = DxI (rounded to nearest Q value)
end

% below is an attempt at extending the temporal range of the flanger, by
% reducing the order of interpolation at the 'edges' such that it is dependent on
% proximity of 'delayed' x to start or end (that is, the closer to an edge,
% the lower the order of interpolation ---> smaller the stencil)

yfflex = x;
for i=2:L-1 % i defines position of undelayed 'stencil', whose size is governed by N
    DxI = i-M0*(1+sin(2*pi*f0*i*Ts)); % exact 'delayed' x location
    if DxI>1
        if DxI < ceil(N/2)
            N = 2*floor(DxI); % stencil size restricted at start edge
            M = interptab(N,Q);
            LGM = zeros(Q,N);
        elseif DxI > L-ceil(N/2)
            N = 2*(L-ceil(DxI)); % stencil size restricted at end edge
            M = interptab(N,Q);
            LGM = zeros(Q,N);
        end
        frac = DxI-round(DxI); % 'displacement' to nearest sample
        if N == 2
            DxIQ = round(frac*Q); % linear interpolation
        else
            DxIQ = round(((floor(N/2)/(N-1))+(frac/(N-1)))*Q); % nearest Q value to DxI's relative position within N-stencil
        end
        if mod(N,2) == 0 % dealing with odd and even N cases seperately
            for j=1:N
                LGM(:,j) = x(floor(DxI)-(N/2)+j)*M(:,j); % specific Langrangian polynomials for each point in stencil calculated here
            end
        else
            for j=1:N
                LGM(:,j) = x(round(DxI)-floor(N/2)+j-1)*M(:,j); % specific Langrangian polynomials for each point in stencil calculated here
            end
        end
        LG = sum(LGM,2); % polynomial which passes through N samples within stencil obtained here
        yfflex(i) = x(i) + LG(DxIQ); % read off polynomial at specific 'alpha' = DxI (rounded to nearest Q value)
    end
    N=N2;
    M=M2;
    LGM=LGM2; % resetting in case of 'edge' calculations
end

%|-----------LAGRANGIAN FEEDBACK COMB FILTER----------|

yfbl = x;
for i=2*M0+ceil(N/2)+1:L-floor(N/2) % i defines position of 'stencil', whose size is governed by N
    DxI = i-M0*(1+sin(2*pi*f0*i*Ts)); % exact 'delayed' x location
    frac = DxI-round(DxI); % 'displacement' to nearest sample
    DxIQ = round(((floor(N/2)/(N-1))+(frac/(N-1)))*Q); % nearest Q value to DxI's relative position within N-stencil
    if mod(N,2) == 0
        for j=1:N
            LGM(:,j) = g*yfbl(floor(DxI)-(N/2)+j)*M(:,j); % specific Langrangian polynomials for each point in stencil calculated here
        end
    else
        for j=1:N
            LGM(:,j) = g*yfbl(round(DxI)-floor(N/2)+j-1)*M(:,j); % specific Langrangian polynomials for each point in stencil calculated here
        end
    end
    LG = sum(LGM,2); % polynomial which passes through N samples within stencil obtained here
    yfbl(i) = x(i) - LG(DxIQ); % read off polynomial at specific 'alpha' = DxI (rounded to nearest Q value)
end


% below is an attempt at extending the temporal range of the flanger, by
% reducing the order of interpolation at the 'edges' such that it is dependent on
% proximity of 'delayed' x to start or end (that is, the closer to an edge,
% the lower the order of interpolation ---> smaller the stencil)

yfblex = x;
for i=2:L-1 % i defines position of undelayed 'stencil', whose size is governed by N
    DxI = i-M0*(1+sin(2*pi*f0*i*Ts)); % exact 'delayed' x location
    if DxI>1
        if DxI < ceil(N/2)
            N = 2*floor(DxI); % stencil size restricted at start edge
            M = interptab(N,Q);
            LGM = zeros(Q,N);
        elseif DxI > L-ceil(N/2)
            N = 2*(L-ceil(DxI)); % stencil size restricted at end edge
            M = interptab(N,Q);
            LGM = zeros(Q,N);
        end
        frac = DxI-round(DxI); % 'displacement' to nearest sample
        if N == 2
            DxIQ = round(frac*Q); % linear interpolation
        else
            DxIQ = round(((floor(N/2)/(N-1))+(frac/(N-1)))*Q); % nearest Q value to DxI's relative position within N-stencil
        end
        if mod(N,2) == 0 % dealing with odd and even N cases seperately
            for j=1:N
                LGM(:,j) = g*yfblex(floor(DxI)-(N/2)+j)*M(:,j); % specific Langrangian polynomials for each point in stencil calculated here
            end
        else
            for j=1:N
                LGM(:,j) = g*yfblex(round(DxI)-floor(N/2)+j-1)*M(:,j); % specific Langrangian polynomials for each point in stencil calculated here
            end
        end
        LG = sum(LGM,2); % polynomial which passes through N samples within stencil obtained here
        yfblex(i) = x(i) - LG(DxIQ); % read off polynomial at specific 'alpha' = DxI (rounded to nearest Q value)
    end
    N=N2;
    M=M2;
    LGM=LGM2; % resetting in case of 'edge' calculations
end

toc

% %|---------DFTs (FOR PLOTTING PURPOSES)----------|
% 
% XF = fft(x);
% XFM2sided = abs(XF/L);
% XFM1sided = XFM2sided(1:L/2+1);
% XFM1sided(2:end-1) = 2*XFM1sided(2:end-1);
% YF_FF = fft(yff);
% YF_FFM2sided = abs(YF_FF/L);
% YF_FFM1sided = YF_FFM2sided(1:L/2+1);
% YF_FFM1sided(2:end-1) = 2*YF_FFM1sided(2:end-1);
% YF_FB = fft(yfb);
% YF_FBM2sided = abs(YF_FB/L);
% YF_FBM1sided = YF_FBM2sided(1:L/2+1);
% YF_FBM1sided(2:end-1) = 2*YF_FBM1sided(2:end-1);
% YF_FFL = fft(yffl);
% YF_FFLM2sided = abs(YF_FFL/L);
% YF_FFLM1sided = YF_FFLM2sided(1:L/2+1);
% YF_FFLM1sided(2:end-1) = 2*YF_FFLM1sided(2:end-1);
% YF_FFLEX = fft(yfflex);
% YF_FFLEXM2sided = abs(YF_FFLEX/L);
% YF_FFLEXM1sided = YF_FFLEXM2sided(1:L/2+1);
% YF_FFLEXM1sided(2:end-1) = 2*YF_FFLEXM1sided(2:end-1);
% YF_FBL = fft(yfbl);
% YF_FBLM2sided = abs(YF_FBL/L);
% YF_FBLM1sided = YF_FBLM2sided(1:L/2+1);
% YF_FBLM1sided(2:end-1) = 2*YF_FBLM1sided(2:end-1);
% YF_FBLEX = fft(yfblex);
% YF_FBLEXM2sided = abs(YF_FBLEX/L);
% YF_FBLEXM1sided = YF_FBLEXM2sided(1:L/2+1);
% YF_FBLEXM1sided(2:end-1) = 2*YF_FBLEXM1sided(2:end-1);
% 
% %|----------PLOTS----------|
% 
% hold on
% subplot(3,3,2)
% plot(f,XFM1sided)
% title('Amplitude spectrum of x (upto 2kHz)')
% xlabel('f (Hz)')
% ylabel('|x|')
% xlim([0 2000])
% subplot(3,3,4)
% plot(f,YF_FFM1sided)
% title('Amplitude spectrum of yff (upto 2kHz)')
% xlabel('f (Hz)')
% ylabel('|yff|')
% xlim([0 2000])
% subplot(3,3,5)
% plot(f,YF_FBM1sided)
% title('Amplitude spectrum of yfb (upto 2kHz)')
% xlabel('f (Hz)')
% ylabel('|yfb|')
% xlim([0 2000])
% subplot(3,3,6)
% plot(f,YF_FFLM1sided)
% title('Amplitude spectrum of yffl (upto 2kHz)')
% xlabel('f (Hz)')
% ylabel('|yffl|')
% xlim([0 2000])
% subplot(3,3,7)
% plot(f,YF_FFLEXM1sided)
% title('Amplitude spectrum of yfflex (upto 2kHz)')
% xlabel('f (Hz)')
% ylabel('|yfflex|')
% xlim([0 2000])
% subplot(3,3,8)
% plot(f,YF_FBLM1sided)
% title('Amplitude spectrum of yfbl (upto 2kHz)')
% xlabel('f (Hz)')
% ylabel('|yfbl|')
% xlim([0 2000])
% subplot(3,3,9)
% plot(f,YF_FBLEXM1sided)
% title('Amplitude spectrum of yfblex (upto 2kHz)')
% xlabel('f (Hz)')
% ylabel('|yfblex|')
% xlim([0 2000])
% 
% %|----------DEBUG - DxIQ CHECK (INSERTED WITHIN FOR LOOP AFTER DxIQ)----------|
% 
% % if DxIQ > (N/2*(N-1))*Q
% %     disp('DxI goes higher than centre-range at:')
% %     disp(i)
% % elseif DxIQ < ((N-2)/2*(N-1))*Q
% %     disp('DxI goes lower than centre-range at:')
% %     disp(i)
% % end