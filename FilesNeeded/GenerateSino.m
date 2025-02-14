function [A,ReconData] = GenerateSino(A,Myphantom)
%% Initialization for PET System

FOV_Width = 300;              %Width of the Field of View
nR = 256;
DetWidth = 4;                 %Width of each detector
nPhi = 288;                   %Number of Projection angles
TotalCounts = 6.8E6;          %Total number of event counts(True+Random+Scatter)
%High counts: 6.8E6;   Low counts: 6.8E5; 1.7E6

PixelWidth = FOV_Width/nR;
nDet = 2*nPhi;
RF = 0.25;                    %Random Fraction
SF = 0.25;                    %Scatter Fraction
AC = 0.096;                   %Attenuation Coefficient

[ LORdist,LORwidth ] = LOR_DistWidth(nDet,DetWidth,FOV_Width);
nLOR = numel(LORdist);

%% Simulation of Spatial Resolution:
%  1.Detector Size; 2.Positron Range;
%  3.Noncolinearity: Residual momentum of the positron
%  4.Decoding: The ability of detector to identify the location of photon
%  5.Penetration: Gamma rays penetrate some distance into the detector ring

SRR = (DetWidth/2)/tan(pi/nDet);    %Short Ring Radius
LRR = (DetWidth/2)/sin(pi/nDet);    %Long Ring Radius
FWHM_Det = DetWidth/2;              %FWHM of Detector Size
FWHM_PosiRange = 0.5;               %FWHM of Positron Range
FWHM_Noncolinear = 0.0044*SRR;      %FWHM of Noncolinearity
FWHM_Decoding = DetWidth/3;         %FWHM of Decoding
FWHM_Penetration = (12.5*FOV_Width/4)/sqrt((FOV_Width/4)^2+LRR^2);   %FWHM of Penetration
%FWHM = sqrt(FWHM_Det^2+FWHM_PosiRange^2+FWHM_Noncolinear^2)/PixelWidth;
FWHM = sqrt(FWHM_Det^2+FWHM_PosiRange^2+FWHM_Noncolinear^2+FWHM_Decoding^2+FWHM_Penetration^2)/PixelWidth;
FWHM = 2*FWHM;
%FWHM = 0;
sigma = FWHM/(2*sqrt(2*log(2)));    %Gaussian FWHM=2*sqrt(2*ln(2))*sigma
fSizeOdd = max(3*ceil(sigma),5);
fSizeOdd = mod(fSizeOdd+1,2) + fSizeOdd;
if FWHM>0
    PSF = fspecial('gaussian',[fSizeOdd fSizeOdd],sigma);
    PSFphantom = imfilter(Myphantom,PSF,'same');
else
    PSFphantom = Myphantom;
end

sino = A*PSFphantom(:);

%% Attenuation Simulation

AttThickness = LOR_Attenuation(Myphantom,nPhi,LORdist,PixelWidth);
AttenMat = reshape(exp(-(AC*AttThickness/10)),[nLOR nPhi]); %Attenuation Matrix
%The projection matrix for attenuation is different from that for PET.
AttenMat(AttenMat>1) = 1;
AttenMat = sparse(AttenMat);

%% Random and Scatter Coincidence Simulation

RandCounts = RF*TotalCounts;
ScatCounts = SF*(TotalCounts-RandCounts);
TrueCounts = TotalCounts-RandCounts-ScatCounts;
ScatFWHM = 2/3*nR;
ScatMat = max(ceil(3*ScatFWHM),3);
if (mod(ScatMat,2)==0)
    ScatMat = ScatMat+1;
end
ScatKernel = fspecial('gaussian',[ScatMat ScatMat],ScatFWHM/(2*sqrt(2*log(2))));
PTscatter = imfilter(PSFphantom,ScatKernel,'replicate','same','conv');
S = reshape(A*PTscatter(:),[nLOR nPhi]);
R = ones(size(S));
G = reshape(sino,[nLOR nPhi]).*AttenMat;
G = full(G);

%% DATA for Reconstuction
G = TrueCounts*G/sum(G(:));
S = ScatCounts*S/sum(S(:));
R = RandCounts*R/sum(R(:));
RS = R+S;

seed = rng(0);
ReconData.sinogram = poissrnd(G+RS);
ReconData.sinogram(ReconData.sinogram<0) = 0;
ReconData.gamma = reshape(RS,[nLOR*nPhi 1]);
FOV = reshape(FieldofView(nR),[nR*nR 1]);
G_NoAtten = full(G./AttenMat);
p = sum(G_NoAtten(:));
ReconData.TMC = p/nPhi/sum(FOV(:));
ReconData.PerfPhantom = p/nPhi*Myphantom/sum(Myphantom(:));

% MRI = p/nPhi*MRI/sum(MRI(:));
%% A_attenuation*A_geometry

tic;
load('Ablur_2.8sigma.mat');
A = sparse(diag(AttenMat(:))*A*Ablur);
toc;
fprintf('It costs %f seconds to calculate A_atten * A_geom * A_blur\n',toc);

ReconData.InitIM = ReconData.TMC*FOV;

end