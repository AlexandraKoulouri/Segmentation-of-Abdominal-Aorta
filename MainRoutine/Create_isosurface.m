 % Create an isosurface
 function []=Create_isosurface(Lumen);
 
 
 sizeLumen = size(Lumen);
 I = find(Lumen>0);
 Lumen(I) = 255;
 clear I J I_min J_min
 for i = 2:sizeLumen(1,3)
           
        [I,J] = find(Lumen(:,:,i)>0);
        I_min(i-1) = min(I); I_max(i-1) = max(I);
        J_min(i-1) = min(J); J_max(i-1) = max(J);
        J = 0; I = 0; 
        
 end

 I_min = min(I_min)-3;
 I_max = max(I_max)+3;
 J_min = min(J_min)-3;
 J_max = max(J_max)+3;
       
TightLumens = zeros(length(I_min:I_max),length(J_min:J_max),length(2:sizeLumen(1,3)));
TightLumens (:,:,:) = Lumen(I_min:I_max,J_min:J_max,2:sizeLumen(1,3));
%contourslice( Lumen, [], [], [2:4:ImgSequenceSize(1,3)-1] , 3 )
%save test2Results
%save 'Tomes' TightLumens

%clear all

%load Tomes
TightLumens = smooth3(TightLumens,'gaussian',[9,9,15],30);
TightLumens =smooth3 (TightLumens);
sizeTightLumens = size(TightLumens);      
      OpThLumen = zeros(size(TightLumens)); 
     
for i = 1:sizeTightLumens(1,3)

    OpThLumen(:,:,sizeTightLumens(1,3)-i+1) =uint8(TightLumens(:,:,i));
    
end

%clear TightLumens

     % X=[0.5:0.5:512*0.5];
      %Y=X;
      %Z=[1:61];%sizeLumen(1,3)];
%      [X,Y,Z] = meshgrid(X,Y,Z);
 
 LumenSurface = isosurface(  OpThLumen,'LineWidth' ,'0.2');
%       
%       % Display the surface
figure;
%subplot(1,2,1);
hiso = patch('Vertices',LumenSurface.vertices,...
                    'Faces',LumenSurface.faces,...
                   'FaceColor','r',...
                   'EdgeColor','none',...
                   'FaceLighting','phong',...
                   'EdgeLighting','phong');
                   %'LineWidth', '3');
 view(-45,5) 
 axis tight 
 daspect([1,1,.4])
lightangle(45,30);
       
set(gcf,'Renderer','zbuffer'); lighting phong; 
%set(hiso,'SpecularColorReflectance',5,'SpecularExponent',50)
view(-70,10)
lightangle(-70,30);
        
 %Reconstruct the volume and display it as montage
%OV = surface2volume(LumenSurface,[],1);
%nDims = size(OV);
%figure
%montage(reshape(OV,nDims(1),nDims(2),1,nDims(3)),[0 1]);      
      