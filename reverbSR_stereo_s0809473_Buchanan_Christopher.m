clear all;

% stereo implementation, starting from reverbSR_ex method

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

if size(x,2)== 1
    x = repmat(x,2); % this script works with stereo input; if mono is given, it'll be copied into L and R channels here
end

Nyq = Fs/2;

%AP1

y1 = allpassM_s0809473_Buchanan_Christopher(x,0.7,1051);

%AP2

y2 = allpassM_s0809473_Buchanan_Christopher(y1,0.7,337);

%AP3

y3 = allpassM_s0809473_Buchanan_Christopher(y2,0.7,113);

%FFCF1 - FFCF4 %stereo aspect implemented here with differing M values)

y4(:,1) = feedforcombex_s0809473_Buchanan_Christopher(y3(:,1),[0.742,0.733,0.715,0.697],[4500,5000,5200,5800]);
y4(:,2) = feedforcombex_s0809473_Buchanan_Christopher(y3(:,2),[0.742,0.733,0.715,0.697],[5000,4500,5500,5400]);

%normalize and output

yfinal4 = y4/norm(y4,Inf);

[~,name,ext] = fileparts(filename);

audiowrite([name,'_st_reverbed',ext],yfinal4,Fs);

