function product=dotprod_2dbeams(mode1,mode2,crossDot)

if nargin==2 || (nargin==3 && ~crossDot)
    product=sum(mode1.*mode2,3);
elseif crossDot
    mode1_E=mode1.E; mode2_E=mode2.E;
    mode1_H=mode1.H; mode2_H=mode2.H;
    product=mode1_E(:,:,1).*mode2_H(:,:,2)-mode1_E(:,:,2).*mode2_H(:,:,1);
% elseif nargin==4 && ignoreEz
%     product=mode1_E(:,:,1).*mode2_E(:,:,1)+mode1_E(:,:,2).*mode2_E(:,:,2);
end

