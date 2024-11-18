function [DistanceImg] = distanceFunction1 (c1,c2,R,sizeI,ro);
 

[Img] = DrawCircle(R,c1,c2,sizeI(1,1),sizeI(1,2));

FilledImg = fillingProc(c1,c2,Img);

DistanceImg = 2*Img - FilledImg; 

Io = find ( DistanceImg ==0);

I1 = find( DistanceImg > 0);

I_1 = find( DistanceImg <0);
    
DistanceImg =ro* ones(sizeI);


DistanceImg(I1) = 0;

DistanceImg(I_1) = -ro;
