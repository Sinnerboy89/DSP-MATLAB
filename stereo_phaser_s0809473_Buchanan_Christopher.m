function [] = stereo_phaser_s0809473_Buchanan_Christopher(~)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% 
% This is a function to read in a file (supplied by the user after a file browser prompt) and perform
% a quadrature stereo phaser effect using 4 concatenated 2nd order allpass filters.
%
% This function will work with all audio files audioread/audiowrite
% supports.
%
% NOTE: if input is mono, crude stereo duplication will occur before
% phasing
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Have user browse for a file, starting at current location.
currentFolder = pwd;

% Get the full path of the audio file that the user wants to phase.
defaultfilename = fullfile(currentFolder, '*.*');
[filename, folder] = uigetfile(defaultfilename, 'Pick something to Eddie Van Halen-ize');
if filename == 0  % user clicked cancel (maybe because of stage fright).
    return;
end
filepathin = fullfile(folder, filename);

[x,Fs] = audioread(filepathin);
if size(x,2)== 1
    x = repmat(x,2); % this script works with stereo input; if mono is given, it'll be copied into L and R channels here
end

% Nyq = Fs/2;

% N = length(x);

rate = 1;

% Pole location

R = 0.9;
notchfreq = 300;
theta = 2*pi*notchfreq/Fs;
delta = 3*theta;

% Apply time-varying allpass filters

yl1 = allpassTV_s0809473_Buchanan_Christopher(x(:,1),R,theta,delta,rate,Fs);
yr1 = allpassTV_quad_s0809473_Buchanan_Christopher(x(:,2),R,theta,delta,rate,Fs);

R = 0.98;
notchfreq = 800;
theta = 2*pi*notchfreq/Fs;

yl2 = allpassTV_s0809473_Buchanan_Christopher(yl1,R,theta,delta,rate,Fs);
yr2 = allpassTV_quad_s0809473_Buchanan_Christopher(yr1,R,theta,delta,rate,Fs);

R = 0.8;
notchfreq = 1000;
theta = 2*pi*notchfreq/Fs;

yl3 = allpassTV_s0809473_Buchanan_Christopher(yl2,R,theta,delta,rate,Fs);
yr3 = allpassTV_quad_s0809473_Buchanan_Christopher(yr2,R,theta,delta,rate,Fs);

R = 0.9;
notchfreq = 4000;
theta = 2*pi*notchfreq/Fs;

yl4 = allpassTV_s0809473_Buchanan_Christopher(yl3,R,theta,delta,rate,Fs);
yr4 = allpassTV_quad_s0809473_Buchanan_Christopher(yr3,R,theta,delta,rate,Fs);

g=1;

ylfinal = yl4+g.*x(:,1);
yrfinal = yr4+g.*x(:,2);

yfinal = [ylfinal,yrfinal];

yfinal = yfinal/norm(yfinal,Inf);

soundsc(yfinal,Fs)

[~,name,ext] = fileparts(filename);

audiowrite([name,'_phased',ext],yfinal,Fs);