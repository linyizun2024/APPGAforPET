IM = SHOTVAPPGAIM4;
minIM = min(IM(:));
IM  = IM-minIM;
IM = IM/max(IM(:));
IM = 1-IM;
figure,imshow(IM,[0.1,1]);