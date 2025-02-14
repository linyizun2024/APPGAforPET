function [IM,FOM] = ReconSHOTVAPPGA(A,w,ReconData,ReconParam)
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
object = ReconData.PerfPhantom(:); 
MAX = max(object);
ITER = ReconParam.ITER;
lambda1 = ReconParam.lambda1;
lambda2 = ReconParam.lambda2;
beta = ReconParam.beta;
epsilon = ReconParam.epsilon;

% calculate Phi0
% Af= A*f;
% penalty1 = FirstOrderHuberITV(f,epsilon);
% penalty2 = SecondOrderHuberITV(f,epsilon);
% penalty  = lambda*(penalty1+penalty2);
% OFV0 = sum(Af)-sum(g.*log(Af+gamma))+penalty;

Phi_ref = -3.647043531565344e+07;
Phi0 = -3.285870429522608e+07;

preAT1 = AT1;
preAT1(preAT1<=0) = 1;


FOM.OFV = zeros(ITER,1);    
FOM.NOFV = zeros(ITER,1);  
FOM.RE = zeros(ITER,1);     %Relative Error
FOM.MSE = zeros(ITER,1);    %Mean Squar Error
FOM.RMSE = zeros(ITER,1);   %Root Mean Squar Error
FOM.PSNR = zeros(ITER,1);   %Peak Signal to Noise Ratio 
FOM.IterTime = zeros(ITER,1);

%% Reconstruction
 
a = 1/8;
b = 1;
t = 1; 
ORI_f = f;
for k = 1:ITER
    tic;
%--------------------------------------------------------------------------
%Update f
%-------------------------------------------------------------------------- 
    ORI_t = t;
    t = a*(k^w)+b; 
    theta_k = (ORI_t-1)/t;
    ftd_k = f+theta_k*(f-ORI_f);
    ORI_f = f;
    Af = A*ftd_k;    
    g_Af_gmma = g./(Af+gamma);
    q = AT1-A'*g_Af_gmma;      
    B1f = FirstOrderDiff(ftd_k); 
    B2f = SecondOrderDiff(ftd_k);
    gradg_atB1f = gradSITV_B1f(B1f,epsilon);
    gradg_atB2f = gradSITV_B2f(B2f,epsilon);
    subgrad = q+lambda1*FirstOrderDiffTrans(gradg_atB1f)+lambda2*SecondOrderDiffTrans(gradg_atB2f);
    precond =ftd_k./preAT1;   
    f = ftd_k - beta*precond.*subgrad;
    f(f<0) = 0;
 
%------------------------------ --------------------------------------------      
    FOM.IterTime(k)=toc;
    if mod(k,50)==0
        fprintf('SHOITV-APPGA: The %dth iteration costs %f seconds\n',k,FOM.IterTime(k));
    end
        
    Af = A*f;
    penalty1 = FirstOrderSITV(f,epsilon);
    penalty2 = SecondOrderSITV(f,epsilon);
    penalty  = lambda1*penalty1+lambda2*penalty2;
    FOM.OFV(k) = sum(Af)-sum(g.*log(Af+gamma))+penalty;
    FOM.NOFV(k) = (FOM.OFV(k)-Phi_ref)/(Phi0-Phi_ref);
    FOM.RE(k) = norm(f-ORI_f,2)/norm(f,2);
    FOM.MSE(k) = (norm(object-f,2)^2)/dimf;
    FOM.PSNR(k) = 10*log10(MAX*MAX/FOM.MSE(k));
    FOM.RMSE(k) = sqrt(FOM.MSE(k));
    FOM.NRMSE(k) = (norm(object-f,2))/(norm(object,2));
end
fprintf('The total time of %d SHOITV-APPGA iterations is: %f seconds\n',ITER,sum(FOM.IterTime(:)));
fprintf('The average time of SHOITV-APPGA for each iteration is: %f seconds\n',sum(FOM.IterTime(:))/ITER);
IM = reshape(f,[nR,nR]);
%OFVf = min(FOM.OFV);           %calculate Phi_ref (we use the value when k=1000)
end
