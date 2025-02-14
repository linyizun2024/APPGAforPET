function [ LORdist,LORwidth ] = LOR_DistWidth(NumDet,DetWidth,ImWidth)

angle = 2*pi/NumDet;
RingRadius = (DetWidth/2)/(tan(angle/2));  %The distance from the center of the Ring to the center of the detector
LongRadius = (DetWidth/2)/(sin(angle/2));  %The distance from the center of the Ring to the left or right boundary of the detector
MaxImWidth = 2*RingRadius/sqrt(2);
if ImWidth<=MaxImWidth
    tmpLORdist = zeros(1,((NumDet-2)/2-1)/2);
    tmpLORwidth = zeros(1,((NumDet-2)/2-1)/2);
    count = 0;
    dist = 0;
    while dist<=(ImWidth/2)
        tmpLORdist(count+1) = RingRadius*sin((count+1)*angle);
        dist = LongRadius*sin((count+1)*angle+angle/2);
        tmpLORwidth(count+1) = 2*(dist-tmpLORdist(count+1));
        count = count+1;
    end
    count = count-1;
    LORdist = zeros(1,2*count+1);
    LORwidth = zeros(1,2*count+1);
    for i=1:count
        j = count+1-i;
        LORdist(i) = -tmpLORdist(j);
        LORwidth(i) = tmpLORwidth(j);
    end
    LORdist(count+1) = 0;
    LORwidth(count+1) = DetWidth;
    LORdist(count+2:2*count+1) = tmpLORdist(1:count);
    LORwidth(count+2:2*count+1) = tmpLORwidth(1:count);
else
    LORdist=0;
    LORwidth =0;
    fprintf('Error: The width of the image to be reconstructed is too large\n');
end

end
    
