% Dimension of subimage array
subimageDim=325;  % must be larger than 100
% imageFalse: array with all zeros
imageFalse = zeros(subimageDim,subimageDim);
% imageTrue: array with all ones
imageTrue  =  ones(subimageDim,subimageDim);
% Create a black border of width = 1% of subimage 
border = ones(size(imageFalse));
width = fix(0.01*subimageDim);
border([1:width],:)=0; % top edge 
border(:,[1:width])=0; % left edge
border((end-width:end),:)=0; % bottom edge
border(:,(end-width:end))=0; % right edge
% Insert border
imageFalse = imageFalse & border;
imageTrue  = imageTrue  & border;
% Create composite image array inputA:
% Image inputA = [ False, False;
%                  True,  True ]  
inputA=[[imageFalse,imageFalse] ; ...
        [imageTrue ,imageTrue ]];
% Create composite image array inputB:
% Image inputB = [ False, True;
%                  False, True ]  
inputB=[[imageFalse,imageTrue]; ...
        [imageFalse,imageTrue]];
subplot(1,2,1); imshow(inputA); title('A');  % imshow shows image
subplot(1,2,2); imshow(inputB); title('B');  % imshow shows image