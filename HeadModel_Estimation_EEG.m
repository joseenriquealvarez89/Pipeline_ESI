%% cleaning all...
clc;
clear all;
close all;
warning off;

cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));

%% loading Anatomical (T1) and MEG data...
disp('******************************************');
disp('*   loading anatomical and EEG data...   *');
disp('******************************************');
[filename_channel, pathname] = uigetfile({'dataEEG\'},'Channel file for the EEG data:');
pathElectPos=[pathname,filename_channel];
% load([pathname,filename_channel]);      % Channels data
% elect = [Channel.Loc]';

MriFile = uigetdir('dataEEG\','Pick the Freesurfer Output folder:');

% MriFile='E:\ComparingInverseProblem\LeadfieldPipeline\PipelinePedroArielOriginal(copiaTania)\Pipeline ESI\anatomy\MC0000010_t13d_anatVOL_20060115002658_2.nii_out';
% 
% pathElectPos='channel_ASA_10-05_343.mat';


nVertices        = 10000;                                % must be defined by the user...
resamplingMethod = 'reducepatch';
erodeFactor      = 1; %%????
fillFactor       = 2; %%????
headVertices     = 1922;
isEEG            = 1;
SnrFixed         = 3;
NoiseReg         = 0.1;
conductivity     = 0.0125;
showFigure       = 1;



[Gain,InverseS] = pipeline_brainstorm_fs(MriFile,pathElectPos,nVertices,...
    resamplingMethod,erodeFactor,fillFactor,headVertices,...
    isEEG,SnrFixed,NoiseReg,conductivity,showFigure);

