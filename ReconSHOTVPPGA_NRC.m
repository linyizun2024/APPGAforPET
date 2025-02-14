function [IM,FOM] = ReconSHOTVPPGA_NRC(A,ReconData,ReconParam)
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
ITER = ReconParam.ITER;
lambda1 = ReconParam.lambda1;
lambda2 = ReconParam.lambda2;
beta = ReconParam.beta;
epsilon = ReconParam.epsilon;

preAT1 = AT1;
preAT1(preAT1<=0) = 1;

FOM.CNR1 = zeros(ITER,1);     %Contrast to Noise Ratio
FOM.CNR2 = zeros(ITER,1);
FOM.RC1 = zeros(ITER,1);      %Relative Contrast
FOM.RC2 = zeros(ITER,1);
FOM.NRC1 = zeros(ITER,1);     %Normalized Relative Contrast
FOM.NRC2 = zeros(ITER,1);
FOM.IR = zeros(ITER,1);       %Image Roughness

FOM.IterTime = zeros(ITER,1);
center = ceil(nR/2);
DiskRadius = ceil(3*nR/8);
Hot_dist = ceil(DiskRadius/2);

ROIradius1 = 4;
HotRow = center + ceil(sqrt(3)*Hot_dist/2);
HotCol = center + ceil(Hot_dist/2);
ROI1 = FieldofView(2*ROIradius1);
ROI1_PixelNum = numel(find(ROI1(:)));
HotROI1_RowFirst = HotRow-ROIradius1+1;
HotROI1_RowLast = HotRow+ROIradius1;
HotROI1_ColFirst = HotCol-ROIradius1+1;
HotROI1_ColLast = HotCol+ROIradius1;

ROIradius2 = 14;
ROI2 = FieldofView(2*ROIradius2);
ROI2_PixelNum = numel(find(ROI2(:)));
HotROI2_RowFirst = center-ROIradius2+1;
HotROI2_RowLast = center+ROIradius2;
HotROI2_ColFirst = center+Hot_dist-ROIradius2+1;
HotROI2_ColLast = center+Hot_dist+ROIradius2;

BackROI1_RowFirst = center-ROIradius1+1;
BackROI1_RowLast = center+ROIradius1;
BackROI1_ColFirst = center-ROIradius1+1;
BackROI1_ColLast = center+ROIradius1;

BackROI2_RowFirst = center-ROIradius2+1;
BackROI2_RowLast = center+ROIradius2;
BackROI2_ColFirst = center-ROIradius2+1;
BackROI2_ColLast = center+ROIradius2;
%% Reconstruction
Af= A*f;
for k = 1:ITER
    tic;
%--------------------------------------------------------------------------
%Update f 
%--------------------------------------------------------------------------    
    g_Af_gmma = g./(Af+gamma);
    q = AT1-A'*g_Af_gmma;     
    B1f = FirstOrderDiff(f); 
    B2f = SecondOrderDiff(f);
    gradg_atB1f = gradSITV_B1f(B1f,epsilon);
    gradg_atB2f = gradSITV_B2f(B2f,epsilon);
    subgrad = q+lambda1*FirstOrderDiffTrans(gradg_atB1f)+lambda2*SecondOrderDiffTrans(gradg_atB2f);  
    precond =f./preAT1; 
    f = f - beta*precond.*subgrad;
    f(f<0) = 0;
%--------------------------------------------------------------------------      
    Af = A*f;
    
    FOM.IterTime(k)=toc;
    if mod(k,50)==0
        fprintf('SHOITV-PPGA: The %dth iteration costs %f seconds\n',k,FOM.IterTime(k));
    end
     
    IM = reshape(f,[nR nR]);
    SampleHot1 = zeros(ROI1_PixelNum,1);
    NumCount = 0;
    for i = HotROI1_RowFirst:HotROI1_RowLast
        for j = HotROI1_ColFirst:HotROI1_ColLast
            if (ROI1(i-HotROI1_RowFirst+1,j-HotROI1_ColFirst+1)>0)
                NumCount = NumCount+1;
                SampleHot1(NumCount) = IM(i,j);
            end
        end
    end
    
    SampleHot2 = zeros(ROI2_PixelNum,1);
    NumCount = 0;
    for i = HotROI2_RowFirst:HotROI2_RowLast
        for j = HotROI2_ColFirst:HotROI2_ColLast
            if (ROI2(i-HotROI2_RowFirst+1,j-HotROI2_ColFirst+1)>0)
                NumCount = NumCount+1;
                SampleHot2(NumCount) = IM(i,j);
            end
        end
    end
    
    SampleBack1 = zeros(ROI1_PixelNum,1);
    NumCount = 0;
    for i = BackROI1_RowFirst:BackROI1_RowLast
        for j = BackROI1_ColFirst:BackROI1_ColLast
            if (ROI1(i-BackROI1_RowFirst+1,j-BackROI1_ColFirst+1)>0)
                NumCount = NumCount+1;
                SampleBack1(NumCount) = IM(i,j);
            end
        end
    end
    
    SampleBack2 = zeros(ROI2_PixelNum,1);
    NumCount = 0;
    for i = BackROI2_RowFirst:BackROI2_RowLast
        for j = BackROI2_ColFirst:BackROI2_ColLast
            if (ROI2(i-BackROI2_RowFirst+1,j-BackROI2_ColFirst+1)>0)
                NumCount = NumCount+1;
                SampleBack2(NumCount) = IM(i,j);
            end
        end
    end
    MeanHotROI1 = sum(SampleHot1(:))/ROI1_PixelNum;
    MeanHotROI2 = sum(SampleHot2(:))/ROI2_PixelNum;
    MeanBackROI1 = sum(SampleBack1(:))/ROI1_PixelNum;
    MeanBackROI2 = sum(SampleBack2(:))/ROI2_PixelNum;
    SDBackROI1 = std(SampleBack1);
    SDBackROI2 = std(SampleBack2);
    
    FOM.CNR1(k) = (MeanHotROI1-MeanBackROI1)/SDBackROI1;
    FOM.CNR2(k) = (MeanHotROI2-MeanBackROI2)/SDBackROI2;
    FOM.RC1(k) = (MeanHotROI1-MeanBackROI1)/MeanBackROI1;
    FOM.RC2(k) = (MeanHotROI2-MeanBackROI2)/MeanBackROI2;
    FOM.NRC1(k) = FOM.RC1(k)/3;
    FOM.NRC2(k) = FOM.RC2(k)/3;
    FOM.IR(k) = sum((SampleBack2-MeanBackROI2).^2)/(ROI2_PixelNum-1);
end
fprintf('The total time of %d SHOITV-PPGA iterations is: %f seconds\n',ITER,sum(FOM.IterTime(:)));
fprintf('The average time of SHOITV-PPGA for each iteration is: %f seconds\n',sum(FOM.IterTime(:))/ITER);  
IM = reshape(f,[nR,nR]);
end
