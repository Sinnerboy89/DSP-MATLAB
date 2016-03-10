function [] = phaser_s0809473_Buchanan_Christopher(~)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% 
% This is a function to read in a file (supplied by the user after a file browser prompt) and perform
% a phaser effect using 4 concatenated 2nd order allpass filters.
%
% This function will work with all audio files audioread/audiowrite
% supports.
%
% NOTE: if input is stereo, crude conversion to mono will occur before
% phasing and output
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
if size(x,2)>1
    x = (x(:,1)+x(:,2))/2;
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

y1 = allpassTV_s0809473_Buchanan_Christopher(x,R,theta,delta,rate,Fs);

R = 0.98;
notchfreq = 800;
theta = 2*pi*notchfreq/Fs;

y2 = allpassTV_s0809473_Buchanan_Christopher(y1,R,theta,delta,rate,Fs);

R = 0.8;
notchfreq = 1000;
theta = 2*pi*notchfreq/Fs;

y3 = allpassTV_s0809473_Buchanan_Christopher(y2,R,theta,delta,rate,Fs);

R = 0.9;
notchfreq = 4000;
theta = 2*pi*notchfreq/Fs;

y4 = allpassTV_s0809473_Buchanan_Christopher(y3,R,theta,delta,rate,Fs);

g=1;
yfinal = y4+g.*x;

yfinal = yfinal/norm(yfinal,Inf);

soundsc(yfinal,Fs)

[~,name,ext] = fileparts(filename);

audiowrite([name,'_mono_phased',ext],yfinal,Fs);