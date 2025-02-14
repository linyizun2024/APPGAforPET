clear;
close all;
clc;

%% Choose the Phantom
PhanNum = 1;

if PhanNum == 1
    load('Brain256.mat');
    Myphantom = Brain256;
elseif PhanNum == 2
    load('UniformPhantom.mat');
    Myphantom = UniformPhantom;
else
    error('PhanNum can only be 1 or 2!!')
end

%% Parameters for Reconstruction
ReconParam.ITER = 100;
ReconParam.beta = 1;            % beta=0.1 for Uniform Phantom
ReconParam.epsilon = 0.001;
ReconParam.nSubsets = 8;
ReconParam.s = 2;
ReconParam.t = 8;


if PhanNum == 1
    ReconParam.lambda1 = 0.04; 
    ReconParam.lambda2 = 0.04;
else
    ReconParam.lambda1 = 0.4;
    ReconParam.lambda2 = 0;
end


%% Load the geometry projection matrix A
matName = 'Amat_256x256_77LORx288P_300FOV_4DW_32Ray.mat';
isFile = exist(matName, 'file');
if (isFile == 2)
    tic;
    load(matName);
    LoadA=toc;
    fprintf('It costs %f seconds to load the system matrix A\n',LoadA);
else
    error('!!! No matrix available !!!')
end

%% Generate system matrix and sinogram
[A,ReconData] = GenerateSino(A,Myphantom);
%%
[HOTVKMIM,HOTVKMFOM] = ReconHOTVPKMA(A,ReconData,ReconParam);
[SHOTVPPGAIM,SHOTVPPGAFOM] = ReconSHOTVPPGA(A,ReconData,ReconParam);
[SHOTVAPPGAIM1,SHOTVAPPGAFOM1] = ReconSHOTVAPPGA(A,1/4,ReconData,ReconParam);
[SHOTVAPPGAIM2,SHOTVAPPGAFOM2] = ReconSHOTVAPPGA(A,1/2,ReconData,ReconParam);
[SHOTVAPPGAIM3,SHOTVAPPGAFOM3] = ReconSHOTVAPPGA(A,3/4,ReconData,ReconParam);
[SHOTVAPPGAIM4,SHOTVAPPGAFOM4] = ReconSHOTVAPPGA(A,1,ReconData,ReconParam);

%% Show the reconstructed images and plots
maxIM = max(ReconData.PerfPhantom(:));
figure;imshow(ReconData.PerfPhantom,[0 maxIM]);title('Original Phantom');
figure;imshow(HOTVKMIM,[0 maxIM]);title(['HOTV-PKMA (' num2str(ReconParam.ITER) ' Iterations)']);
figure;imshow(SHOTVPPGAIM,[0 maxIM]);title(['SHOTV-PPGA (' num2str(ReconParam.ITER) ' Iterations)']);
figure;imshow(SHOTVAPPGAIM1,[0 maxIM]);title(['SHOTV-APPGA(\omega=1/4) (' num2str(ReconParam.ITER) ' Iterations)']);
figure;imshow(SHOTVAPPGAIM2,[0 maxIM]);title(['SHOTV-APPGA(\omega=1/2) (' num2str(ReconParam.ITER) ' Iterations)']);
figure;imshow(SHOTVAPPGAIM3,[0 maxIM]);title(['SHOTV-APPGA(\omega=3/4) (' num2str(ReconParam.ITER) ' Iterations)']);
figure;imshow(SHOTVAPPGAIM4,[0 maxIM]);title(['SHOTV-APPGA(\omega=1) (' num2str(ReconParam.ITER) ' Iterations)']);

plotfunc(HOTVKMFOM,SHOTVPPGAFOM,SHOTVAPPGAFOM1,SHOTVAPPGAFOM2,SHOTVAPPGAFOM3,SHOTVAPPGAFOM4,ReconParam.ITER,PhanNum);