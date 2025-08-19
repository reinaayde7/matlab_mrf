function [mask] = Fillin_ProstateMask(mask)
for i=1:size(mask,2)
    X=mask(:,i);
    whites=find(X);
    max_white=max(whites);
    min_white=min(whites);
    for j=min_white:max_white
        mask(j,i)=1;
    end
end
for j=1:size(mask,1)
    X=mask(j,:);
    whites=find(X);
    max_white=max(whites);
    min_white=min(whites);
    for i=min_white:max_white
        mask(j,i)=1;
    end
end

end
%{
How to use this function:

tmpimage = squeeze(abs(recon_images)); % 
tmp_mask = zeros(size(tmpimage));
factor = 1; % this factor can be adjusted based on the intensity of the tempimage
tmp_mask(tmpimage>factor*mean(tmpimage(:))) = 1;
[tmp_mask_filled] = Fillin_ProstateMask(tmp_mask); % calling the function
%}