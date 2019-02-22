function [Gain,InverseSolution]=pipeline_brainstorm_fs(MriFile,electrodPos, nVertices,resamplingMethod,erodeFactor,fillFactor,headVertices, isEEG, SnrFixed,NoiseReg,conductivity, showFigure )
%Pipeline ESI Readme
% INPUT:
% -MriFile               : Mri Freesufer File Path
% -electrodPos           : Electrode Position Matrix
% -nVertices             : Number of vertices for the cortex
% -resamplingMethod      : Resampling Method options: 'reducepatch' ||  'iso2mesh'
% -erodeFactor           : Parameter for convexifing the scalp surface, ranges from 0 (min) to 3 (max)(default: 1)
% -fillFactor            : Parameter for filling holes in the scalp surface ranges from 0 (min) to 3 (max)(default: 2)
% -headVertices          : Number of vertices on the estimated scalp surface(default: 1922)
% -isEEG                 : 1 for EEG or 0 for MEG (default: 1)
% -conductivity          : parameter modelling the conductivity of the skull surface, relative to the conductivities of the scalp and the cortex surfaces (default: 0.0125)
% -NoiseReg              : Noise Covariance Regularization (default:0.1)
% -SnrFixed              :Signal to noise ratio 1/lambda (default:3)
% -showFigure            :Produce Graphics
%End Pipeline ESI Readme

Brainstorm_route = cd;
my_addjava(Brainstorm_route);

%Validation Parameters
if (nargin < 12) || isempty(showFigure)  
    showFigure = 0;
end

if (nargin < 11) || isempty(conductivity)  
    conductivity = 0.0125;
end

if (nargin < 10) || isempty(NoiseReg)
    NoiseReg = 0.1;
end

if (nargin < 9) || isempty(SnrFixed)
      SnrFixed = 3;
end

if (nargin < 8) || isempty(isEEG) || isEEG>0 
    isEEG = 1;
end

if (nargin < 7) || isempty(headVertices)
    headVertices = 1922;
end

if (nargin < 6) || isempty(fillFactor)
    fillFactor = 2;
end

if (nargin < 5) || isempty(erodeFactor)
    erodeFactor =1;
end

if (nargin < 4) || isempty(resamplingMethod) || (not(strcmp(resamplingMethod,'reducepatch' ))&& not(strcmp(resamplingMethod,'iso2mesh')))
    resamplingMethod = 'reducepatch';
end

if (nargin < 3) || isempty(nVertices)
   nVertices = 10000;
end

if (nargin < 2) || isempty(electrodPos) || isempty(MriFile)
   error('Error, electrodPos or/and MriFile isempty, you must define the path to load it')
end

load(electrodPos)

%Fiducial point estimation
sMri = my_bst_normalize_mni(MriFile);

%Importing the anatomy
[sMri, Cortex ] = my_import_anatomy_fs( MriFile,sMri, nVertices, 0, [], 0, resamplingMethod);

%Cleaning the surfaces
sHead = my_tess_isohead(sMri, 10000, erodeFactor, fillFactor);

%Skull Interpolation
nvert =[headVertices 1922 1922];
thickness= [7 4 3];
[ sInner,sOuter,sHeadBEM ] = bem_surfaces_brainstorm(sMri,sHead,Cortex,nvert,thickness,Brainstorm_route);

elect = [Channel.Loc]';
%Electrode fitting
if(isEEG)
    [locsx] = channel_project_scalp(sHeadBEM.Vertices,elect);
else
    locsx = elect;
end
 

%Checking electrode positions
if(showFigure)
  head.vc = sHeadBEM.Vertices;
  head.tri = sHeadBEM.Faces;
  figure;showsurface(head,[],locsx);
end


%BEM
OPTIONS.GridLoc = Cortex.Vertices;
OPTIONS.isEeg = isEEG;
OPTIONS.isMeg = not(OPTIONS.isEeg);
OPTIONS.BemSurf = cell(3,1);
OPTIONS.BemSurf{3} = sInner;
OPTIONS.BemSurf{2}= sOuter;
OPTIONS.BemSurf{1} = sHeadBEM;


OPTIONS.BemNames = {'Scalp'    'Skull'    'Brain'};
OPTIONS.BemCond = [1.0000    conductivity    1.0000];

OPTIONS.isAdjoint = 1;
OPTIONS.isAdaptative = not(OPTIONS.isAdjoint);
OPTIONS.GridOrient = Cortex.VertNormals;

OPTIONS.Channel.Loc = locsx;
iVertInside = find(inpolyhd(locsx', sInner.Vertices, sInner.Faces));
Gain = my_bst_openmeeg(OPTIONS);
iBad = find(any(isnan(OPTIONS.GridOrient),2) | any(isinf(OPTIONS.GridOrient),2) | (sqrt(sum(OPTIONS.GridOrient.^2,2)) < eps));
if ~isempty(iBad)
    OPTIONS.GridOrient(iBad,:) = repmat([1 0 0], length(iBad), 1);
end

%Inverse Solution
OPTIONS = [];
if (isEEG==1)
  NoiseCov = diag(1e-10*ones(size(locsx,1),1));
  OPTIONS.ChannelTypes = repmat({'EEG'},[1,size(locsx,1),]);
  OPTIONS.DataTypes = {'EEG'};
else
    NoiseCov = diag(1e-15*ones(size(locsx,1),1));
    OPTIONS.ChannelTypes = repmat({'MEG'},[1,size(locsx,1),]);
    OPTIONS.DataTypes = {'MEG'};
end    

OPTIONS.NoiseCovMat.NoiseCov = NoiseCov;
OPTIONS.NoiseCovMat.nSamples  = [];
OPTIONS.NoiseCovMat.FourthMoment = [];
OPTIONS.InverseMethod = 'minnorm';
OPTIONS.InverseMeasure = 'sloreta';
OPTIONS.SourceOrient = {'fixed'};
OPTIONS.ComputeKernel = 1;
OPTIONS.Loose = 0.2;
OPTIONS.UseDepth = 0;
OPTIONS.WeightExp = 0.5;
OPTIONS.WeightLimit = 10;
OPTIONS.NoiseMethod = 'reg';
OPTIONS.NoiseReg = NoiseReg;
OPTIONS.SnrMethod = 'fixed';
OPTIONS.SnrRms = 1e-06;
OPTIONS.SnrFixed = SnrFixed;


HeadModel.Gain = Gain;
HeadModel.GridOrient = Cortex.VertNormals;
HeadModel.GridLoc = Cortex.Vertices;
HeadModel.HeadModelType = 'surface';


[Results, OPTIONS] = bst_inverse_linear_2016(HeadModel, OPTIONS);
InverseSolution=Results.ImagingKernel;

end

 
