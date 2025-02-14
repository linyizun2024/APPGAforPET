function [IM,FOM] = ReconHOTVPKMA(A,ReconData,ReconParam)
%% Input Variables
%A:The System Matrix;
%Data: 1.g(sinogram); 2.gamma(the sum of Random and Scatter);
%      3.InitIM;      4.PerfPhantom(Perfect Phantom).
%ReconParam: 1.ITER; 2.lambda1; 3.lambda2.

%% Output Variables
%IM:The reconstucted image
%FOM(Figure of Merit):
%1.OFV(Objective function Value)
%2.RE(Relative Error)
%3.MSE(Mean Squar Error)
%4.PSNR(Peak Signal to Noise Ratio)
%5.RMSE(Root Mean Squar Error)

%% Initialization
AT1 = A'*ones([size(A,1) 1]);
f = ReconData.InitIM;
dimf = numel(f);
nR = sqrt(dimf);
nLOR = size(ReconData.sinogram,1);
nPhi = size(ReconData.sinogram,2);
g = reshape(ReconData.sinogram,[nLOR*nPhi 1]);
gamma = ReconData.gamma;
object = reshape(ReconData.PerfPhantom,[dimf 1]);
ITER = ReconParam.ITER;
beta = ReconParam.beta;
MAX = max(object);
lambda1 = ReconParam.lambda1;
lambda2 = ReconParam.lambda2;
delta = 0.1;
eta = 0.9;

% calculate Phi0
% Af= A*f;
% penalty1 = FirstOrderITV(f);
% penalty2 = SecondOrderITV(f);
% penalty  = lambda1*penalty1+lambda2*penalty2;
% OFV0 = sum(Af)-sum(g.*log(Af+gamma))+penalty;

preAT1 = AT1;
preAT1(preAT1<=0) = 1;

Phi_ref = -3.647042927589290e+07;
Phi0 = -3.285870424788608e+07;

b = zeros(2*dimf,1);
c = zeros(4*dimf,1);
FOM.OFV = zeros(ITER,1);    %Objective Funtion Value
FOM.RE = zeros(ITER,1);     %Relative Error
FOM.MSE = zeros(ITER,1);    %Mean Squar Error
FOM.RMSE = zeros(ITER,1);   %Root Mean Squar Error
FOM.PSNR = zeros(ITER,1);   %Peak Signal to Noise Ratio
FOM.IterTime = zeros(ITER,1);

%% Reconstruction
Af = A*f;
for k = 1:ITER
    tic;
%--------------------------------------------------------------------------
%Update f
%-------------------------------------------------------------------------- 
    ORI_f = f;
    g_Af_gmma = g./(Af+gamma);
    q = AT1-A'*g_Af_gmma;
    subgrad = q+FirstOrderDiffTrans(b)+SecondOrderDiffTrans(c);
    precond = f./preAT1;
    f = ORI_f - beta*precond.*subgrad;
    f(f<0) = 0;

%-------------------------------------------------------------------------- 
%Update b and c
%--------------------------------------------------------------------------
    MaxS = max(precond);
    rho1 = 1/(2*8*MaxS);
    rho2 = 1/(2*64*MaxS);
    
    ORI_b = b;
    ORI_c = c;
    
    B1f = FirstOrderDiff(2*f-ORI_f);
    b_B1f = 1/rho1*b+B1f;
    z1 = b_B1f(1:dimf);
    z2 = b_B1f(dimf+1:2*dimf);
    z3 = sqrt(z1.^2+z2.^2);
    b_prox = zeros(2*dimf,1);
    ii = find(z3>lambda1/rho1);
    b_prox(ii) = (1-lambda1/rho1*(1./z3(ii))).*z1(ii);
    b_prox(dimf+ii) = (1-lambda1/rho1*(1./z3(ii))).*z2(ii);
    b = rho1*(b_B1f - b_prox);
    
    
    B2f = SecondOrderDiff(2*f-ORI_f);
    c_B2f = 1/rho2*c+B2f;     
    z1 = c_B2f(1:dimf);
    z2 = c_B2f(dimf+1:2*dimf);
    z3 = c_B2f(2*dimf+1:3*dimf);
    z4 = c_B2f(3*dimf+1:4*dimf);
    z5 = sqrt(z1.^2+z2.^2+z3.^2+z4.^2);   
    c_prox = zeros(4*dimf,1);
    ii = find(z5>lambda2/rho2);
    c_prox(ii) = (1 - lambda2/rho2*(1./z5(ii))).*z1(ii);
    c_prox(dimf+ii) = (1 - lambda2/rho2*(1./z5(ii))).*z2(ii);
    c_prox(2*dimf+ii) = (1 - lambda2/rho2*(1./z5(ii))).*z3(ii);
    c_prox(3*dimf+ii) = (1 - lambda2/rho2*(1./z5(ii))).*z4(ii);
    c = rho2*(c_B2f - c_prox);
    
    theta = eta*(k-1)/(k+delta-1);
    f = f+theta*(f-ORI_f);
    b = b+theta*(b-ORI_b);
    c = c+theta*(c-ORI_c);
%--------------------------------------------------------------------------        
    Af = A*f;
    
    FOM.IterTime(k)=toc;
    if mod(k,50)==0
        fprintf('HOITV-PKMA: The %dth iteration costs %f seconds\n',k,FOM.IterTime(k));
    end
    
    penalty1 = lambda1*FirstOrderITV(f);
    penalty2 = lambda2*SecondOrderITV(f);
    FOM.OFV(k) = sum(Af)-sum(g.*log(Af+gamma))+penalty1+penalty2;
    FOM.NOFV(k) = (FOM.OFV(k)-Phi_ref)/(Phi0-Phi_ref);
    FOM.RE(k) = norm(f-ORI_f,2)/norm(f,2);
    FOM.MSE(k) = (norm(object-f,2)^2)/dimf;
    FOM.PSNR(k) = 10*log10(MAX*MAX/FOM.MSE(k));
    FOM.RMSE(k) = sqrt(FOM.MSE(k));
    FOM.NRMSE(k) = (norm(object-f,2))/(norm(object,2));
end
fprintf('The total time of %d HOITV-PKMA iterations is: %f seconds\n',ITER,sum(FOM.IterTime(:)));
fprintf('The average time of HOITV-PKMA for each iteration is: %f seconds\n',sum(FOM.IterTime(:))/ITER);
IM = reshape(f,[nR nR]);
% OFVf = min(FOM.OFV);   %calculate Phi_ref (we use the value when k=1000)
end

