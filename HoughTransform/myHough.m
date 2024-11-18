function[Radius, c1,c2] = myHough(image,R,sigma,thres,ROI,control)
Img = image;%rgb2gray((image));
InImg = Img;

ROI = round(ROI);
if nargin ==5
%Pre-processing
%determine Region of Interest(aortic region) for a [512x512]CT
%image(170:389,170:350)
Img = Img(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));
%filtering using median filter
%Img = medfilt2 ( Img, [9,9] ) ;
end
    
%keep only the higher intensity values
%Img(Img<max(max(Img))-2)=0;
%Egde detetion using canny

[EImg,filterSize]= edge1(Img,'canny',thres,sigma);%,[80, 200],3);

%%%%%%%%%%%%Hough Transform%%%%%%%%%%%%%%%%%%%%
%imshow(EImg)
%pause(0.1)
%Determine Edges

[x,y] = find(EImg);
[sx,sy] = size(EImg);

totalpix = length(x);


%Number of possible Radius values
Rnum = length(R);
R_2(1,1,:) = R.^2;
R_2 = repmat(R_2,[sy,totalpix,1]);
y = repmat(y',[sy,1,Rnum ]);
x = repmat(x',[sy,1,Rnum ]);

% Preparing all the matrices for the computation without "for-loop"
b = 1:sy;
a = zeros(sy,totalpix);

b = repmat(b',[1,totalpix,Rnum]);
b2 = b;
%
% The equation for the circle (estimate cirlce centrers)
a = (round(x - sqrt(R_2 - (y - b).^2)));
a2 = (round(x + sqrt(R_2 - (y - b2).^2)));
%a = a1;
a = cat(1,a,a2);
b=cat(1,b,b2);
  
clear b2 a2
[xx,yy,zz]=find((imag(a)==0 & a>0 & a<sx));
% Removing all the invalid value in matrices a and b

b = b(imag(a)==0 & a>0 & a<sx);
a = a(imag(a)==0 & a>0 & a<sx);

%ind1 = sub2ind([sy,sx],b1,a1);%Position of a and b in the image(circle coordinates)
%ind2 = sub2ind([sy,sx],b2,a2);

%ind = [ind1; ind2];
ind = sub2ind([sy,sx],b,a);
ind = ind + (ceil(yy/totalpix)-1)*sx*sy;
clear b a xx yy zz x y totalpix
%Layers of Hough transformed Image
HImg = zeros(sy*sx*Rnum ,1);
val = ones(length(ind),1);%
data = accumarray(ind,val);
HImg(1:length(data)) = data;
HImg2 = reshape(HImg,[sy,sx,Rnum ]);
clear HImg Img EImg
%find the higher value of each layer for every radius
cnt = 1;
H(cnt) = max(max(HImg2(:,:,cnt)));
cnt = cnt + 1;
H(cnt) = max(max(HImg2(:,:,cnt)));
%while  H(cnt-1)<H(cnt)+2 && cnt <Rnum;
 % cnt = cnt +1;
  % H(cnt) = max(max(HImg2(:,:,cnt)));
%end
for i =1:length(R)
   H(i) = max(max(HImg2(:,:,i)));
   %imshow(HImg2(:,:,i),[])
   %pause
end

%if length(I)>1
 %     H(I(2:length(I)))=0;
%end
%figure
%plot(R,H,'*-')
[maxval, maxind] = max(H);
[B,A] = find(HImg2(:,:,maxind)==maxval);
clear HImg2 

%bwh = im2bw(HImg2(:,:,maxind), 0.67);
%bwh = imdilate(bwh, ones(3, 3));
%bwh = imerode(bwh, ones(3, 3));
%bwh = bwlabel(bwh);
%figure(4)
%imshow(uint8(HImg2(:,:,maxind)))
SizeI = size(InImg);
%center of cicle
c1 =round(ROI(1,1)+ mean(A));%x
c2 =round(ROI(2,1)+ mean(B));%y
Radius = R(maxind);
[CircleImg] = DrawCircle(Radius,c1,c2,SizeI(1,1),SizeI(1,2),2);

%Represent the results

[x1,map1] = gray2ind (CircleImg,12) ;
rgb1 = ind2rgb(x1,map1) ;
rgb1(:,:,1) = 0 ;
rgb1(:,:,3) = 0 ;
[x2, map2] = gray2ind(InImg ,120);
rgb2 = ind2rgb(x2,map2) ;
Frame  = zeros (SizeI(1,1),SizeI(1,2));
ROI = round(ROI);
Frame(ROI(1,1),ROI(2,1):ROI(2,2))=255;
Frame(ROI(1,2),ROI(2,1):ROI(2,2))=255;
Frame(ROI(1,1):ROI(1,2),ROI(2,1))=255;
Frame(ROI(1,1):ROI(1,2),ROI(2,2))=255;
[x1,map1] = gray2ind (Frame ,12);
rgb3 = ind2rgb(x1,map1) ;
rgb3(:,:,3) = 0 ;
rgb3(:,:,2) = 0 ;

%subplot(1,2,1)
%imshow( rgb2+rgb1+rgb3 );
%pause(0.5)
%hold;plot(c2,c1,'xr')
%text(c2,c1,num2str(Radius),'color','green')
%hold;
%subplot(1,2,2)
%imshow( rgb1 );
%hold;
%plot(c2,c1,'xr')
%text(c2,c1,num2str(maxind),'color','green')
%hold;