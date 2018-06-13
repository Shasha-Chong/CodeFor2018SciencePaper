%Shasha Chong
%Sept, 2016
clear; clc;
%This script sorts all the single-particle trajectories according to their
%xy positions relative to EWS/FLI1 cluster mask generated from the 561 channel image stack

%Locate the name and path of the 561 channel image stack
clfname = '561nm_corrected.tif_NucleusClusterMasked';
clpath = '/Volumes/Data2/JF_NikonScope/170721/KO116_HE-EFS_PAJF646_20nM_CPJF549_200nM_cell1/';

%Extract the number of images in the 561 channel image stack
clfinfo=imfinfo([clpath, clfname, '.tif']);
stacklength=numel(clfinfo);

%Specify the microscope integration time in seconds:
exposure = 0.50; %Frame rate (time in seconds)
resolution = 0.160; %um/pixel

climage = struct;
for i = 1:stacklength
    climage(i).xyf = imread([clpath, clfname, '.tif'], i); 
    climage(i).TimeStamp = i*exposure;
    climage(i).Frame = i;
end

%Load the 633 channel trajectories
trajname = 'All';
trajpath = '/Volumes/Data2/JF_NikonScope/170721/KO116_HE-EFS_PAJF646_20nM_CPJF549_200nM_cell1/Data/';
load([trajpath, trajname, '.mat']);

%Categorize trajectories according their locations relative to clusters
numtraj = sum(~cellfun(@isempty,{trackedPar.xy}));
cltraj = struct;
cltrajnum = 0;
ncltraj = struct;
ncltrajnum = 0;

%Extract the information of each trajectory 
for j = 1:numtraj
    trajlength = length(trackedPar(j).TimeStamp);
    trajxy = trackedPar(j).xy./resolution;
    trajframe = trackedPar(j).Frame;
    clfreq = 0;
    climintensity = zeros(trajlength,1);
    for k = 1:trajlength
        trajx = round(trajxy(k,1));
        trajy = round(trajxy(k,2));
        climintensity(k) = climage(trajframe(k)).xyf(trajy,trajx);
    end
    clfreq = sum(climintensity);
    %A cluster trajectory is defined to interact with a cluster for trajlength/2
    %times
    if clfreq > trajlength*0.5;
       cltrajnum = cltrajnum+1;
       cltraj(cltrajnum).xy = trackedPar(j).xy;
       cltraj(cltrajnum).TimeStamp = trackedPar(j).TimeStamp;
       cltraj(cltrajnum).Frame = trackedPar(j).Frame;
       cltraj(cltrajnum).Freq = clfreq;
       cltraj(cltrajnum).TrajIndex = j;
    else if clfreq < trajlength*0.05;
            ncltrajnum = ncltrajnum+1;
            ncltraj(ncltrajnum).xy = trackedPar(j).xy;
            ncltraj(ncltrajnum).TimeStamp = trackedPar(j).TimeStamp;
            ncltraj(ncltrajnum).Frame = trackedPar(j).Frame;
            ncltraj(ncltrajnum).Freq = clfreq;
            ncltraj(ncltrajnum).TrajIndex = j;
        end
    end
end

%Calculate the average length of trajectories within or outside clusters
cltrajleng = zeros(cltrajnum,1);
ncltrajleng = zeros(ncltrajnum,1);
for g = 1:cltrajnum
    cltrajleng(g) = length(cltraj(g).TimeStamp);
end
avecltrajleng = mean(cltrajleng);

for g = 1:ncltrajnum
    ncltrajleng(g) = length(ncltraj(g).TimeStamp);
end
avencltrajleng = mean(ncltrajleng);

%Save the information of in-cluster and out-of-cluster trajectories
save([trajpath, trajname, '_ClusterTraj_v4.mat'], 'cltraj');
save([trajpath, trajname, '_NonClusterTraj_v4.mat'], 'ncltraj');

