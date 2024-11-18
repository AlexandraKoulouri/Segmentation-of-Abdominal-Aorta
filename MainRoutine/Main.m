clear
close all

x= input('enter the data directory-name  :','s'); %the directory which has the tomographies (e.g. data\test3)

struct=dir(x);
cd (x);
tic
% Preallocate the array

% Load image sequence

for i = 5:length(struct)-1
    k = struct(i,1).name;
    image = rgb2gray(imread(k));
    ImgSequence(:,:,i-4) = image;
end
ImgSequenceSize = size(ImgSequence);
 
%Define directories
%cd('..')
cd('..')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRE-PROCESSING
% Process image sequence using 2 threshold and doing rescaling
i_lower = double(0.2*max(max(ImgSequence(:,:,1))));
i_upper = double(0.85*max(max(ImgSequence(:,:,1))));

for i = 1:ImgSequenceSize(1,3)
    
    image=double(ImgSequence(:,:,i));
    image(image>i_upper) = 255;
    image(image<i_lower) = 0;
    indeces = find(image>i_lower & image<i_upper);
    image(indeces) = (double(image(indeces))-i_lower)/(i_upper-i_lower)*255;
    image = medfilt2 ( image, [9,9] ) ;
    ImgSequence(:,:,i) =image;
    indeces = 0;
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1st stage: KALMAN PROCESS%
% Kalman filter initialization

R=[3 0 0;0 3 0;0 0 8];
%R=[1800 0 0; 0 1800 0; 0 0 6000];
H=[[1,0,0]',[0,1,0]',[0,0,1]',[0,0,0]',[0,0,0]',[0,0,0]'];
P = eye(6);
dt=0;%dt=1;
Q=5*dt/6 * kron ([2*dt^2 3*dt;3*dt 6],eye(3));
Q=[3 0 0 0 0 0;0 3 0 0 0 0;0 0 15 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0];
A=[[1,0,0,0,0,0]',[0,1,0,0,0,0]',[0,0,1,0,0,0]',[dt,0,0,1,0,0]',[0,dt,0,0,1,0]',[0,0,dt,0,0,1]'];
x=zeros(size(ImgSequence,3),6);
c = zeros(size(ImgSequence,3),2);
Radius = zeros(size(ImgSequence,3),1);

ROI = [150,350; 150,330];
RR = 18:45;
sigma = 1.3;
thres = [0.2,0.5];
[Radius1, c1,c2] = myHough(ImgSequence(:,:,1),RR,sigma,thres,ROI);
clear RR
c(1,1)=c1;
c(1,2)=c2;
x(1,1) = c1;
x(1,2)= c2;
x(1,3)= Radius1;
Radius(1)=Radius1;

xp=A*x(1,:)';% initial estimation

PP = P;
Lumen = zeros(size( ImgSequence));

%Kalman filtering for each image
 for i = 2 : size(ImgSequence,3)
     ROI = [xp(1,1)-round(2*xp(3,1)),xp(1,1)+round(1.2*xp(3,1));xp(2,1)-round(2*xp(3,1)),xp(2,1)+round(2*xp(3,1))];
     RR =round(xp(3,1))-4:round(xp(3,1))+3;
     
     [Radius(i), c(i,1),c(i,2)] = myHough(ImgSequence(:,:,i),RR,sigma,thres,ROI);%Hough transform->Kalman measurement
     K = PP*H'*inv(H*PP*H'+R);                              %kalman gain 
     x(i,:) = (xp + K*([c(i,1),c(i,2),Radius(i)]' - H*xp))';%update estimation
     P = (eye(6)-K*H)*PP;                                   %update error covariance matrix
     xp=A*x(i,:)';                                          %Project Ahead step
     PP = A*P*A' + Q;                                       %Estimate of the error covariance matrix
     %figure(1)
     %displayResults(ImgSequence(:,:,i),xp,c(i,:),Radius(i),ROI);
     
 
 end

 
%2nd Stage: LEVEL SET METHOD
 
for i = 1 : size(ImgSequence,3)

    [u] = MainFast2DlevelSet(ImgSequence(:,:,i),round(x(i,1)),round(x(i,2)),round(x(i,3)/2));
    Indeces = find (smooth(u)<=0);
    tempLumen = zeros(size(u));
    tempLumen(Indeces) = 255;
    Lumen(:,:,i) =  edge(tempLumen,'canny');
end
clear tempLumen;
%3rd Stage: SHAPE CORRECTION

 for i= 2 : size(ImgSequence,3)
     
     Img=ControlAortaCurvatur(Lumen(:,:,i), round(x(i,2)),round(x(i,1)));
   
     if length(Img)>1
         Lumen(:,:,i)=Img;
         
     end
     clear Img
     finalImg(:,:,i-1) = fillingProc(round(x(i,1)),round(x(i,2)),Lumen(:,:,i));

 end
toc
%Create_isosurface( finalImg);
% for i = 1 : size(ImgSequence,3)
% 
% [sensitivity(i),specificity(i),DICE(i)]=ValidationFunction(finalImg(:,:,i),ImgSequence(:,:,i+1),1);
% end
 

