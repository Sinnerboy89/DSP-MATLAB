clear all;

% This script adds mono reverb to an audio file of the user's choosing,
% outputting it with a '_reverbed' suffix into same folder as input. 4 FFC
% filters are implemented here as specified in block diagram.

% Have user browse for a file, starting at current location.
currentFolder = pwd;

% Get the full path of the audio file that the user wants to Schroeder.
defaultfilename = fullfile(currentFolder, '*.*');
[filename, folder] = uigetfile(defaultfilename, 'Pick something');
if filename == 0  % user clicked cancel (maybe because of stage fright).
    return;
end
filepathin = fullfile(folder, filename);

[x,Fs] = audioread(filepathin);

if size(x,2)>1
    x = (x(:,1)+x(:,2))/2;
    disp('input is stereo - converted to mono')
end

Nyq = Fs/2;

%AP1

y1 = allpassM_s0809473_Buchanan_Christopher(x,0.7,1051);

%AP2

y2 = allpassM_s0809473_Buchanan_Christopher(y1,0.7,337);

%AP3

y3 = allpassM_s0809473_Buchanan_Christopher(y2,0.7,113);

%FFCF1

y4 = feedforcomb_s0809473_Buchanan_Christopher(y3,0.742,4799);

%FFCF2

y5 = feedforcomb_s0809473_Buchanan_Christopher(y3,0.733,4999);

%FFCF3

y6 = feedforcomb_s0809473_Buchanan_Christopher(y3,0.715,5399);

%FFCF4

y7 = feedforcomb_s0809473_Buchanan_Christopher(y3,0.697,5801);

%Combination

y8 = y4 + y5 + y6 + y7;

%normalize and output

yfinal = y8/norm(y8,Inf);

[~,name,ext] = fileparts(filename);

audiowrite([name,'_reverbed',ext],yfinal,Fs);