function seg=SegResultShow(Img,uu,th,k)
%--------------------------------------------------------------------------
%   This program is showing the segmentaion result of our Two-Stage method.
%    
%   Usage: seg=SegResultShow(Img,uu,th,k);
%
%   Inputs: 
%       - Img: the given image
%       - uu: result of our Stage-one
%       - th: thresholds got at our Stage-two
%       - k: the phases number
%
%   Outputs: 
%       - seg: k phases representation of our segmentation result
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------



% show the contours of each phase
for i=1:k-1
    figure,imshow(Img,'border','tight');
    hold on;
    if 1==i
        temp=uu<th(i);
    elseif k-1==i
        temp=(uu>=th(i-1)).*(uu<th(i));
        temp1=uu>=th(i);       
    else
        temp=(uu>=th(i-1)).*(uu<th(i));   
    end
    contour(temp,'r');  
end

if k>2
    figure,imshow(Img,'border','tight');
    hold on;
    contour(temp1,'r'); 
end


% show all the phases
[xLen,yLen]=size(Img);
seg=zeros(xLen,yLen);
for i=1:k-1
    if 1==i
        temp=uu<th(i);
        seg=seg+temp*sum(sum(temp.*uu))/sum(temp(:)+1e-9);
        if 2==k
            temp=uu>=th(i);
            seg=seg+temp*sum(sum(temp.*uu))/sum(temp(:)+1e-9); 
        end
    elseif k-1==i
        temp=(uu>=th(i-1)).*(uu<th(i));
        seg=seg+temp*sum(sum(temp.*uu))/sum(temp(:)+1e-9);
        temp=uu>=th(i);
        seg=seg+temp*sum(sum(temp.*uu))/sum(temp(:)+1e-9);        
    else
        temp=(uu>=th(i-1)).*(uu<th(i));
        seg=seg+temp*sum(sum(temp.*uu))/sum(temp(:)+1e-9);    
    end
end
figure,imshow(seg,'border','tight');



