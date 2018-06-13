%Anders Sejr Hansen
%Oct, 2015
clear; clc;
%This script reads in Single-Particle Tracking data from the SLIMfast
%program and re-formats the tab-delimited file.

filename = 'inmask';

%Save in a local MATLAB folder as well
MATLAB_path = '/Volumes/Data2/JF_NikonScope/170707/SJ_Halo-FUS_PAJF646_2nM_YFP-FUS_SlowTracking_cell5_30C/Data/';
MATLAB_name = [filename, '.mat'];
save_val = 1;

%Specify the microscope integration time in seconds:
exposure = 0.50; %Frame rate (time in seconds)
resolution = 0.160; %um/pixel

%Paths:
Path = '/Volumes/Data2/JF_NikonScope/170707/SJ_Halo-FUS_PAJF646_2nM_YFP-FUS_SlowTracking_cell5_30C/';
%Path = '/Users/AndersSejrHansen/Documents/DATA/Microscopy/20160104_JM8N4_MFM_JF549_500pM/H2b-HALO_mESC/';
%Path = '/Users/AndersSejrHansen/Dropbox/Data/Microscopy/Test_Data/20151009_mESC_Clone87_500ms/Using_10-6_Error/';

%Column 1: x-coordinate (pixels)
%Column 2: y-coordinate (pixels)
%Column 3: Frame number
%Column 4: Trajectory number

%data = dlmread([Path, '20151117_HALOmCTCF_C87_JF549_L20_#11_800ms_.rpt_tracked.txt']);
data = dlmread([Path, filename, '.txt']);

%Save two version:
%Save one for InferenceMAP using this format:  
%   Column 1: trajectory number
%   Column 2: x-coordinate (micro meter)
%   Column 3: y-coordinate (micro meter)
%   Column 4: Time stamp (s)

%Save another one for your own analysis code
InfMAP_data = zeros(size(data,1),4);
for i=1:size(data,1)
    InfMAP_data(i,1) = data(i,4);
    InfMAP_data(i,2) = data(i,1)*resolution;
    InfMAP_data(i,3) = data(i,2)*resolution;
    InfMAP_data(i,4) = data(i,3)*exposure;
end

%write file for InferenceMAP
dlmwrite([Path, filename, 'InfMap_data.txt'], InfMAP_data, 'delimiter','\t');



%Now save the data in a MATLAB convenient format
trackedPar = struct;
iter = 1;
for i=1:size(data,1)
    new = 0;
    %Determine whether the new trajectory is equal to the previous one
    NewTraj = data(i,4); %Trajectory number
    if i == 1
        CurrTraj = NewTraj;
        new = 1;
    else
        if NewTraj ~= CurrTraj;
            CurrTraj = NewTraj;
            iter = iter + 1;
            new = 1;
        end
    end
    
    %Save the information
    
    if new == 1
        trackedPar(iter).xy = [data(i,1)*resolution, data(i,2)*resolution];
        trackedPar(iter).TimeStamp = data(i,3)*exposure;
        trackedPar(iter).Frame = data(i,3);
    else
        trackedPar(iter).xy = [trackedPar(iter).xy; data(i,1)*resolution, data(i,2)*resolution];
        trackedPar(iter).TimeStamp = [trackedPar(iter).TimeStamp, data(i,3)*exposure];
        trackedPar(iter).Frame = [trackedPar(iter).Frame data(i,3)];
    end
end

%Calculate the average length of trajectories within or outside clusters
numtraj = sum(~cellfun(@isempty,{trackedPar.xy}));
trajleng = zeros(numtraj,1);
for j = 1:numtraj
    trajleng(j) = length(trackedPar(j).TimeStamp);
end
avetrajleng = mean(trajleng);

save([Path, filename, 'TrackedParticles.mat'], 'trackedPar');


%Save to your local MATLAB folder
if save_val == 1
if exist([MATLAB_path, MATLAB_name])
    temp = load([MATLAB_path, MATLAB_name]);
    temp = temp.trackedPar;
    trackedPar = horzcat(temp, trackedPar);
    save([MATLAB_path, MATLAB_name], 'trackedPar');
else
    save([MATLAB_path, MATLAB_name], 'trackedPar');
end
end