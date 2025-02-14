function [LOR_Att] = LOR_Attenuation(Myphantom,NumAngles,LORdist,PixelWidth)

ProjAngles = -(0:180/NumAngles:180*(1-1/NumAngles));
NumLOR = numel(LORdist);
LOR_Att = zeros(NumLOR*NumAngles,1);
ImSize = size(Myphantom,1);
center = ImSize*PixelWidth/2;
LORcolumn = ceil((LORdist+center)/PixelWidth);
for angle = 1:NumAngles
    RotatedIm = imrotate(Myphantom,ProjAngles(angle),'bilinear','crop');
    for i=1:NumLOR
        ColVector = RotatedIm(:,LORcolumn(i));
        FirstNonzero = 1;
        while ColVector(FirstNonzero)==0&&FirstNonzero<ImSize
            FirstNonzero = FirstNonzero+1;
        end
        if FirstNonzero<ImSize
            LastNonzero = ImSize;
            while ColVector(LastNonzero)==0
                LastNonzero = LastNonzero-1;
            end
            LOR_Att((angle-1)*NumLOR+i) = (LastNonzero-FirstNonzero+1)*PixelWidth;
        else
            LOR_Att((angle-1)*NumLOR+i) = 0;
        end
    end
end

end