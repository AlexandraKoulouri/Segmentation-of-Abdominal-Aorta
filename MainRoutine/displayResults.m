function displayResults(Img,x,c,R,ROI);
x = round(x);
[Img1] = DrawCircle(x(3,1),x(1,1),x(2,1),512,512,4);
[x1,map1] = gray2ind (Img1,12) ;
Img1 = ind2rgb(x1,map1) ;
Img1(:,:,1) = 0 ;
Img1(:,:,3) = 0 ;
%[Img2] = DrawCircle(R, c(1,1),c(1,2),512,512,4);
%[x2,map2] = gray2ind (Img2,12) ;
%Img2 = ind2rgb(x2,map2) ;
%Img2(:,:,2) = 0 ;
%Img2(:,:,3) = 0 ;
[x3, map3] = gray2ind(Img,120);
Img3 = ind2rgb(x3,map3) ;
Frame  = zeros (512,512);
ROI = round(ROI);
Frame(ROI(1,1):ROI(1,1)+3,ROI(2,1):ROI(2,2))=255;
%Frame(ROI(1,2),ROI(2,1):ROI(2,2))=255;
Frame(ROI(1,2):ROI(1,2)+3,ROI(2,1):ROI(2,2))=255;
Frame(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,1)+3)=255;
%Frame(ROI(1,1):ROI(1,2),ROI(2,1)+1)=255;
Frame(ROI(1,1):ROI(1,2),ROI(2,2):ROI(2,2)+3)=255;
%Frame(ROI(1,1):ROI(1,2),ROI(2,2)+1)=255;
[x1,map1] = gray2ind (Frame ,12);
Frame= ind2rgb(x1,map1) ;
Frame(:,:,3) = 0 ;
Frame(:,:,2) = 0 ;

imshow(Img3+Img1+Frame)
hold on;
plot(x(2,1),x(1,1),'o')%,c(1,2),c(1,1),'xr')
%text(c(1,2),c(1,1),num2str(R))
 hold off