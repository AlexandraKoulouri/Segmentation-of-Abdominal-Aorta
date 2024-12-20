function   ImgNew = ControlAortaCurvature(Img,c1,c2)
%function[ImgNew] = ControlAortaCurvatur(Img,c1,c2)
%Input parameters
%Img = filtered Image(aorta contour) 
%[c1,c2] approximation of aorta's centre
%Ouput parameter
%ImgNew : this is the aorta's contour after correction of its morphology
%
[x , y] = ImgEdgesInOrder(Img,c1,c2); 
%[c1 , c2] = centroid (Img);
x = x - c2;
y = y - c1;
PointsNum1=20;
PointsNum2=3;
addx=x(1:PointsNum1-1); 
%xNew = cat(2,x(length(x)-PointsNum2:length(x))',cat(2,x',addx'));
xNew = cat(2,x',addx');
%plot(x,y)
%pause(1)
addy=y(1:PointsNum1-1);
%yNew = cat(2,y(length(y)-PointsNum2:length(y))',cat(2,y',addy'));
yNew = cat(2,y',addy');
dot_prod = zeros(length(x),1);
Vectorr = zeros(length(x),3);PolyGradiet = zeros(length(x),3);

for i =1:length(x)
    
      % limit=180/pi*atan2(yNew(2*PointsNum2+i-1) , xNew(2*PointsNum2+i-1));
       %if limit <0
        %  limit=limit+360;
       %end          
       

       %if (limit >= 45 && limit <= 135) || (limit>=225 && limit<=315) 
           
        %     [k_y,Angle1(i)]=InterpolationPolyn2degree(xNew(i:2*PointsNum2+i-1),yNew(i:2*PointsNum2+i-1));
             
         %    if (limit>=180+45 && limit<=360-45) 
         
          %       PolyGradient(i,:)=3*[cos(Angle1(i)),sin(Angle1(i)),0] ;
             
           %  else

            %     PolyGradient(i,:)=3*[cos(Angle1(i)),sin(Angle1(i)),0];
            
             %end
             %k_x = xNew(i:2*PointsNum2+i-1);
              
       %else

        %   [k_x,Angle1(i)]=InterpolationPolyn2degree(yNew(i:2*PointsNum2+i-1),xNew(i:2*PointsNum2+i-1));
             
         %  if (limit<45 || limit>(360-45) )
            
               %   PolyGradient(i,:) = 3*[sin(Angle1(i)),cos(Angle1(i)),0]; 
              %else
                              
               % PolyGradient(i,:) = 3*[sin(Angle1(i)),cos(Angle1(i)),0]; 
              
              %end
              %k_y = yNew(i:2*PointsNum2+i-1);
       %end

  
      % u = i + (PointsNum2+1)/2;
       
       Vector1(i,:) = [xNew(i+4)-xNew(i) yNew(i+4)-yNew(i) 0];%./MagnitudeVector1
       AugmentedVector(i,:) = [xNew(i+PointsNum1-1)-xNew(i) yNew(i+PointsNum1-1)-yNew(i) 0];
       
       Vectorr(i,:) = [xNew(i) yNew(i) 0];
       %VectorrxPolyGradient(i,:) = cross(Vectorr(i,:),PolyGradient(i,:));
       %VectorrxPolyGradient(i,:)=VectorrxPolyGradient(i,:)/abs(VectorrxPolyGradient(i,3));
      
      % AugmentedVectorxVector1(i,:) = cross (AugmentedVector(i,:),Vector1(i,:));
      % VectorrxVector1(i,:) = cross (Vectorr(i,:),Vector1(i,:));
       theta_Vector1_Vectorr(i)=acos((dot(Vectorr(i,:),Vector1(i,:)))/sqrt(( Vector1(i,1)^2+Vector1(i,2)^2)*( Vectorr(i,1)^2+Vectorr(i,2)^2)));
     %plot(xNew(:),yNew(:),'r',[0,Vectorr(i,1)],[0,Vectorr(i,2)],[xNew(i+4),xNew(i)], [yNew(i+4),yNew(i)],'m')
     %pause(0.05)
     %  plot(xNew,yNew,'b', [xNew(1,u) xNew(1,u+4)] , [yNew(1,i) yNew(1,i+4)],'r',[xNew(1,i) xNew(i+PointsNum1-1)],[yNew(1,i) yNew(i+PointsNum1-1)],'g',[0 xNew(1,u)],[0 yNew(1,i)],'m' , [0,PolyGradient(1,1)],[0, PolyGradient(1,2)],'k',k_x, k_y,'g')
      
       %mult_Cross_prod(i) = VectorrxVector2(3) .* Vector1xVector2(3);% cross product
       %Angle_Dot_prod(i) = dot( Vector, Vector2)/sqrt(( Vector1(1)^2+Vector1(2)^2)*( Vector2(1)^2+Vector2(2)^2));
end
clear i
IndexOfCheckPoint = find ((theta_Vector1_Vectorr*180/pi)<45  |  (theta_Vector1_Vectorr*180/pi)>150) ;
KeepIndeces = 0 ; i=1; xFinal=0; yFinal=0;
if length(IndexOfCheckPoint)>0
    
   while  length(KeepIndeces) <=5   
      
    [KeepIndeces] = CheckPoint(x,y,IndexOfCheckPoint(i));
    
       if length(KeepIndeces)<2*length(x)/3 && length(IndexOfCheckPoint)<2*length(x)/3 && length(IndexOfCheckPoint)>5
                [xFinal, yFinal] = CurveCorrection(x,y,KeepIndeces);  
                %ImgNew(round(x+c1),round(y+c2))=255
                %imshow(ImgNew)
                % ControlAortaCurvatur(ImgNew,c1,c2)
                %end
                i = i + 1;
                xFinal = xFinal + c2;
                yFinal = yFinal +c1;
              
        end
    
    end
end


clear i

if xFinal(1)~=0 &&  yFinal(1)~=0;
    
    ImgNew=zeros(size(Img));
    for i = 1:length(xFinal)
        ImgNew(round(xFinal(i)),round(yFinal(i))) = 1;%255;
       % imshow(ImgNew)
        %pause(0.05)
    end
    figure(1)
    imshow(Img)
    figure(2)
    imshow(ImgNew)
else
    ImgNew=0;
    %disp('no artifacts correction is needed')
end


   
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Function                                             %
%           Put the points in an order                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x , y ] =ImgEdgesInOrder(Img,c1,c2);
%
%this function puts the points of the contour in an order
%input parameters
%Img==>in edged image
%[c1,c2] centre of the contour(centre of the aortic lumen)
%Output parameters
%[x(i),y(i)] coordinates of each point(put them in anticlock order)
%
[x,y] = find (Img>0);
c1=round(c1);c2=round(c2);
IndexMinDistance = FunctionMinDistance(x,y,c1,c2);
%put the edges of the Img in order accoding to the distance between each other
%sizeImg = size(Img); ImgNew = zeros(sizeImg);

%[x_minX,y_minX]=find(x<=min(x));
curx = x(IndexMinDistance(1));%x_minX(1));
cury = y(IndexMinDistance(1));%y_minX(1)); 
x(IndexMinDistance(1)) = 5000;
y(IndexMinDistance(1)) = 5000;
%ImgNew(curx,cury)=255;
NextIndex = [];
if  curx <= c1 && cury <= c2
   
    NextIndex = find (y>=cury & y<=cury+1 & x<=curx & x>=curx-1);
  
end

if  (curx <= c1  && cury >= c2 && length(NextIndex) == 0)
    
     NextIndex = find (y>=cury & y<=cury+1 & x>=curx & x<=curx+1);
end

if curx >= c1 && cury >= c2 &&  length(NextIndex) ==0
  
    NextIndex = find (y<=cury & y>=cury-1 & x>=curx & x<=curx+1);
end

if  curx >=c1 &&  cury <= c2 && length(NextIndex) ==0
    
    NextIndex = find (y<=cury & y>=cury-1 & x<=curx & x>=curx-1);
end


if length(NextIndex)>1
    
    for i = 1 :length(NextIndex)
        
        distance(i) = sqrt(((curx-x(NextIndex(i)))^2)+(cury-y(NextIndex(i)))^2);
    
    end
    
    index = find (distance<=min(distance));
    
    NextPoint = NextIndex (index(1));
else
    
     NextPoint = NextIndex ;
    
end
    
Indeces = zeros(length(x),1);
curx = x(NextPoint);
cury = y(NextPoint);
x(NextPoint) = 5000;
y(NextPoint) = 5000;

Indeces(1,1) = IndexMinDistance(1);
Indeces(2,1) = NextPoint;


distance = zeros(length(1:length(x)),1);
%ImgNew(curx,cury)=255; 
%imshow(ImgNew)



      




for k = 2 :length(x)
    
    
    for i = 1:length(x)
    
        distance(i) = sqrt(((curx-x(i))^2)+(cury-y(i))^2);
     
    end
        index = find (distance<=min(distance));
    
        if length(index)>1
        
            Indeces(k+1) = index(1);
    
        else

            Indeces(k+1) = index;
    
        end
    
    %ImgNew(curx,cury) = 255;
    %imshow(ImgNew)
    %pause(0.5)
    curx = x(Indeces(k+1));
    cury = y(Indeces(k+1));
    x(Indeces(k+1)) = 50000;
    y(Indeces(k+1)) = 50000; 
     
end

[x,y] = find (Img>0);
x = x(Indeces);
y = y(Indeces);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       FunctionMinDistance                     %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IndexMinDistance = FunctionMinDistance(I,J,c1,c2)
%calculates the distance between each contour point and
%the centre of the contour
%and finds the edge with the  minimum distance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input Parameters
%[I,J] points coordinates
%[c1,c2] coordinates of the centre of the contour

MinDistance = zeros(length(I),1);

MinDistance(:) = sqrt((I(:)-c1).^2+(J(:)-c2).^2);

IndexMinDistance =find (MinDistance <= min (MinDistance));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Interpolation                                     %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,Angle1] =InterpolationPolyn2degree(a,b)
%length of 'a' is an odd number
v = length(a);
y=zeros(length(a),1);
Sa = sum(a);
Sa_2 = sum(a.^2);
Sa_3 = sum(a.^3);
Sa_4 = sum(a.^4);
Sb =sum(b);
Sba = sum(a.*b);
Sba_2 = sum(b.*(a.^2));

D=det([v Sa Sa_2;Sa Sa_2 Sa_3;Sa_2 Sa_3 Sa_4]);

D1 = det([Sb Sa Sa_2;
          Sba Sa_2 Sa_3;
          Sba_2 Sa_3 Sa_4]);

D2 = det([v Sb Sa_2;
          Sa Sba Sa_3;
          Sa_2 Sba_2 Sa_4]);

D3 = det( [v Sa Sb;
           Sa Sa_2 Sba;
           Sa_2 Sa_3 Sba_2]);

if D>0 
    C1=D1/D;
    C2=D2/D;
    C3=D3/D;
    grad = C2+C3*2*a(ceil(v/2));
    Angle1 =atan(grad);
   
    %y = C1+C2*a+C3*a.^2;
   
else
  a=a-a(1);
  b=b-b(1);
  Angle1=(atan2(b(ceil(v/2)),a(ceil(v/2))))
  y(:,1)= mean(b);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Function Check Point                               % 
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[KeepIndeces] = CheckPoint(X,Y,IndexOfCheckPoint);
%function[KeepIndeces] = CheckPoint(X,Y,IndexOfCheckPoint);
%
%this function check the curvature of a contour in specific area 
%Input parameters
% [X,Y] contour coordinates
%IndexOfCheckPoint====> [X(IndexOfCheckPoint),Y(IndexOfCheckPoint)]
%of a specific point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CurveLength = length(X);

ChangeSequenceX = zeros(CurveLength,1);
ChangeSequenceY = zeros(CurveLength,1);

if IndexOfCheckPoint>1
    
    ChangeSequenceX (1:CurveLength-IndexOfCheckPoint+1) = X (IndexOfCheckPoint:CurveLength);
    ChangeSequenceX (CurveLength-IndexOfCheckPoint+2:CurveLength) = X (1:IndexOfCheckPoint-1);


    ChangeSequenceY (1:CurveLength-IndexOfCheckPoint+1) = Y (IndexOfCheckPoint:CurveLength);
    ChangeSequenceY (CurveLength-IndexOfCheckPoint+2:CurveLength) = Y (1:IndexOfCheckPoint-1);
else
    
    ChangeSequenceX = X;
    ChangeSequenceY = Y;
end
   
VectorR = [ChangeSequenceX(1),ChangeSequenceY(1),0];
Vectorn = [ChangeSequenceX(4)-ChangeSequenceX(1),ChangeSequenceY(4)-ChangeSequenceY(1),0];
VectorK = -Vectorn;
%theta_Rn =  180/pi* acos(dot (VectorR,Vectorn)/sqrt((VectorR(1)^2+VectorR(2)^2)*(Vectorn(1)^2+Vectorn(2)^2)));
Rxn = cross(VectorR,Vectorn);
i = 4; n = 0;
Kxn =- Rxn;
KeepIndeces = 0;
Theta_Kn =180/pi* acos(dot (Vectorn,Vectorn)/sqrt((Vectorn(1)^2+Vectorn(2)^2)*(Vectorn(1)^2+Vectorn(2)^2)));
k=1;
while (i<CurveLength) && k==1 %Theta_Kn<100

    i=i+1;
    
    VectorK = [ChangeSequenceX(i)-ChangeSequenceX(1), ChangeSequenceY(i)-ChangeSequenceY(1), 0];
   Vectorn_= [ChangeSequenceX(i+2)-ChangeSequenceX(i),ChangeSequenceY(i+2)-ChangeSequenceY(i),0];
   VectorR_=[ChangeSequenceX(i),ChangeSequenceY(i),0];
   Kxn = cross(VectorK,Vectorn);
   
    Theta_Kn =180/pi* acos(dot (VectorK,Vectorn)/sqrt((VectorK(1)^2+VectorK(2)^2)*(Vectorn(1)^2+Vectorn(2)^2)));
    k=0;
% plot(ChangeSequenceX,ChangeSequenceY, [ChangeSequenceX(1),ChangeSequenceX(i)], [ChangeSequenceY(1),ChangeSequenceY(i)],'r',[ChangeSequenceX(1),ChangeSequenceX(4)],[ChangeSequenceY(1),ChangeSequenceY(4)],'m')
 %  pause(1)
   Theta_RK =180/pi* acos(dot (VectorK,VectorR)/sqrt((VectorK(1)^2+VectorK(2)^2)*(VectorR(1)^2+VectorR(2)^2)));   
   Theta_Rn_=180/pi* acos(dot (Vectorn_,VectorR_)/sqrt((Vectorn_(1)^2+Vectorn_(2)^2)*(VectorR_(1)^2+VectorR_(2)^2)));  
   
   if (Rxn(1,3)*Kxn(1,3)>0) || (Theta_Kn <20 && (Rxn(1,3)*Kxn(1,3))<=0 )|| Theta_RK<35 || (Theta_Rn_<50 || Theta_Rn_>145)
        n=n+1;
        KeepIndeces(n)=i;
        k=1;
    end
            
     
end

 if KeepIndeces(1) == 5;
     KeepIndeces= cat(2,[1:4],KeepIndeces);
 end

KeepIndeces = KeepIndeces+IndexOfCheckPoint-1;

I = find (KeepIndeces > CurveLength );
KeepIndeces(I) =KeepIndeces(I)- CurveLength ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%                        Function CurveCorrection                          % 
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            
 function[xxx, yyy] = CurveCorrection(x,y,KeepIndeces);
%
%
%
xxx=y;yyy=x;
LengthInd = length(KeepIndeces);

if KeepIndeces(1) >=2 && KeepIndeces(length(KeepIndeces))<length(x)-1
  
    Radius1 = sqrt (x(KeepIndeces(1)-1)^2+ y(KeepIndeces(1)-1)^2);
    
    Radius2 = sqrt (x(KeepIndeces(LengthInd)+2)^2+ y(KeepIndeces(LengthInd)+2)^2);
    
    a = acosd( dot ([x(KeepIndeces(1)-1) y(KeepIndeces(1)-1)],[x(KeepIndeces(LengthInd)+1) y(KeepIndeces(LengthInd)+1)])/(Radius1* Radius2)); 
    theta1 = 180/pi*atan2(y(KeepIndeces(1)-1),x(KeepIndeces(1)-1));
    theta2 = 180/pi*atan2 (y(KeepIndeces(LengthInd)+2),x(KeepIndeces(LengthInd)+2));

elseif KeepIndeces(1) <=1 && KeepIndeces(length(KeepIndeces))<length(x)-1

   Radius1 = sqrt (x(length(x)-1)^2+ y(length(x)-1)^2);
  
   Radius2 = sqrt (x(KeepIndeces(LengthInd)+2)^2+ y(KeepIndeces(LengthInd)+2)^2);
   
    a = acosd( dot ( [x(length(x)-1) y(length(x)-1)],[x(KeepIndeces(LengthInd)+2) y(KeepIndeces(LengthInd)+2)])/(Radius1* Radius2)); 
    theta1 = 180/pi*atan2(y(length(x)-1),x(length(x)-1));
    theta2 = 180/pi*atan2 (y(KeepIndeces(LengthInd)+2),x(KeepIndeces(LengthInd))+2);
    
%elseif KeepIndeces(1)> 2 && KeepIndeces(length(KeepIndeces))>length(x)-2
    
    
end


if a < 85
    
    if theta1 <0  
        Temptheta1 = 360 + theta1;
    else
        Temptheta1 = theta1;
    end
    if theta2 <0  
        Temptheta2 = 360 + theta2;
    else
        Temptheta2 = theta2;
    
    end
    
    
    if Temptheta1<Temptheta2
      
        Angle1 = Temptheta2 - Temptheta1;
        
        if Angle1>360-Angle1
            %step1= (360-Angle1)/50;
            %step1 = (360-Angle1)/LengthInd;
             % theta = linspace(theta2,theta1,150);
            % theta = (Temptheta2:abs(step1):Temptheta2+(LengthInd-1)*abs(step1)) ; 
            % if Temptheta2 ==  360 + theta2;
                 
             %    theta = theta -360;
             %end
                
                theta = linspace(theta2,theta1,150);
             
        else
            % step1= (360-Angle1)/50;
             %step1 = Angle1/LengthInd; 
            % theta = (Temptheta1:abs(step1):Temptheta1+(LengthInd-1)*abs(step1)) ;  
             % theta = linspace(Temptheta1,Temptheta2,150);
             %if Temptheta1 ==  360 + theta1;
              
              %   theta = theta -360;
                 if (theta1 >=0  && theta2>=0) || (theta1 <=0 && theta2<=0 ) 
                            theta = linspace(theta1,theta2,150);
               end
               if (theta2<0 && theta1>0)
                          theta2 = theta2+360;
                         theta = linspace(theta1,theta2,150);
               end         
              
              
              
              %end
        end
            
    else
        Angle1 = Temptheta1 - Temptheta2;
       
        if Angle1>360-Angle1
             % step1= (360-Angle1)/50;
            %step1 = (360-Angle1)/LengthInd;  
             %theta = (Temptheta1:abs(step1):Temptheta1+(LengthInd-1)*abs(step1)) ;
            %theta = linspace(theta1,theta2,150);
             %if Temptheta1 ==  360 + theta1;
              
              %  theta = theta -360;
               %else
                %   if Temptheta2 == 360+theta2;
                 %      Ind1 = find(theta>180);
                  %     theta(Ind1) = theta(Ind1)- 360;
                  % end
             %end
              
               theta=linspace(theta1,theta2,150);
 
             
             
        else

               if (theta1 >=0  && theta2>=0) || (theta1 <=0 && theta2<=0 ) 
                            theta = linspace(theta2,theta1,150);
               end
               if (theta1<0 && theta2>0)
                           theta1 = theta1+360;
                           theta = linspace(theta2,theta1,150);
               end

                        
               %  step1= (360-Angle1)/50;
            % step1 = Angle1/LengthInd;
            % theta = (Temptheta2:abs(step1):Temptheta2+(LengthInd-1)*abs(step1)) ;
            
               %if Temptheta2 ==  360 + theta2;
              
                % theta = theta -360;
               %else
                %   if Temptheta1 == 360+theta1;
                 %      Ind1 = find(theta>180);
                  %     theta(Ind1) = theta(Ind1)- 360;
                  % end
         
        end

    end


   % Exceed360 = find(theta>360);
    %theta(Exceed360)=theta(Exceed360)-360;
    %[p] = CreateCurve (Radius1,Radius2,theta1,theta2);
%if theta1>0
   [a,b] = CreateCurve (Radius1,Radius2,theta1,theta2);
 %   if  round(a) == round(360 +( theta2 - theta1 ))
    
  %      theta2 = 360 + theta2;
      
   % end
   
 %end

%if theta1<0
    
 %   if  round(a) == round((360-(theta2-theta1)))
       
  %      theta1 = 360 + theta1;
   % end
 
%end
    %Radius = polyval(p,theta);
    Radius = a * theta + b;
    n = Radius'.* cosd (theta)';
    z = Radius'.*sind (theta)';
  %  x(KeepIndeces)=round(n);
   % y(KeepIndeces)=round(z);
   x(KeepIndeces)=50000;
   y(KeepIndeces)=50000;
   Indx_0 =find(x~=50000);
   Indy_0 = find(x~=50000);
   xx=x(Indx_0);
   yy=y(Indy_0);
   xxx=cat(1,xx,round(n));

   yyy=cat(1,yy,round(z));
  
 end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%function[p] = CreateCurve (R1,R2,theta1,theta2);
function[a,b] = CreateCurve (R1,R2,theta1,theta2);

%y=ax+b
%linear case
a = (R1 - R2)/(theta1 - theta2);

b = R1 - a *theta1;
%y=ax^2+bx+c
%R3 = (R1+R2)/2+std([R1,R2])/2;%0.01*R2;
%R=[R1,R3,R2]; 
%theta = [theta1,(theta1+theta2)/2,theta2];
%p = polyfit(theta,R,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
