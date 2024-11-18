%function [finalImg]=fillAorta(i,j,aortaValue,filteredImg)
%function finalImg=fillAorta(i,j,filteredImg,inputImg)
function finalImg = fillingProc(centre_i,centre_j,filteredImg);


%maxPixelValue=250;
setValue=max(max(filteredImg));
finalImg=zeros(size(filteredImg));
 
%threshold choice
threshold=max(max(filteredImg));%mean(mean(filteredImg));

centre_i=centre_i;
centre_j=centre_j;

stack=-ones(255^2,2);
stack_index=1;
stack(stack_index,1)=centre_i;
stack(stack_index,2)=centre_j;
i=centre_i;j=centre_j;
while(stack_index>0 && j+2<size(filteredImg,1) && j-2>0 && i+2<size(filteredImg,2) && i-2>0)
    %oldValue=finalImg(pixel_i,pixel_j);
    i=stack(stack_index,1);
    j=stack(stack_index,2);
    stack_index=stack_index-1;
    finalImg(i,j)=setValue;

    %if( checkNeighbour(filteredImg(i,j),filteredImg(i,j+1))==1 && finalImg(i,j+1)~=setValue )
    if( finalImg(i,j+1)~=setValue )
        
                   if( checkNeighbour(filteredImg(i,j),filteredImg(i,j+1),threshold)==1 )
                         stack_index=stack_index+1;
                         stack(stack_index,1)=i;
                         stack(stack_index,2)=j+1;
                   else
                         finalImg(i,j+1)=setValue;
                   end
     
    end
       
    %if( checkNeighbour(filteredImg(i,j),filteredImg(i,j-1))==1 && finalImg(i,j-1)~=setValue )
    if( finalImg(i,j-1)~=setValue )
        
                    if( checkNeighbour(filteredImg(i,j),filteredImg(i,j-1),threshold)==1 )
                         stack_index=stack_index+1;
                         stack(stack_index,1)=i;
                         stack(stack_index,2)=j-1;
                         %fillAorta(i,j-1);
                    else
                        finalImg(i,j-1)=setValue;
                    end
        
    end

    %if( checkNeighbour(filteredImg(i,j),filteredImg(i+1,j))==1 && finalImg(i+1,j)~=setValue )
    if( finalImg(i+1,j)~=setValue )
        
                    if( checkNeighbour(filteredImg(i,j),filteredImg(i+1,j),threshold)==1 )
                        stack_index=stack_index+1;
                        stack(stack_index,1)=i+1;
                        stack(stack_index,2)=j;
                        %fillAorta(i+1,j)
                    else
                        finalImg(i+1,j)=setValue;
                    end
      
    end

    %if( checkNeighbour(filteredImg(i,j),filteredImg(i-1,j))==1 && finalImg(i-1,j)~=setValue )
    if( finalImg(i-1,j)~=setValue )
        
                if( checkNeighbour(filteredImg(i,j),filteredImg(i-1,j),threshold)==1 )
                       stack_index=stack_index+1;
                       stack(stack_index,1)=i-1;
                       stack(stack_index,2)=j;
                       %fillAorta(i-1,j)
                else
                       finalImg(i-1,j)=setValue;
                end
      
    end
  
end

function isValid=checkNeighbour(curPixel,nextPixel,threshold)
isValid=1;
if (curPixel==threshold) 
    isValid=0;
end
if(abs(nextPixel-curPixel)>0) 
    isValid=0;
end




