for j=1:size(t1map,3)
    figure;
    imagesc(abs(squeeze(t1map(:,:,j))));
    axis image;
    rois{j} = get_rois;
end

%%
i=1;
for j=1:size(t1map,3)
    if ~isempty(rois{1,j}) 
        t1 = t1map(:,:,j);
        t1v(i) = 0.5*(mean(t1(rois{1,j}{1,1})) + mean(t1(rois{1,j}{1,2}))); 
        t2 = t2map(:,:,j);
        t2v(i) = 0.5*(mean(t2(rois{1,j}{1,1})) + mean(t2(rois{1,j}{1,2}))); 
        i=i+1;
    end
end
%%

T1avg = mean(t1v);
T1std = std(t1v);
fprintf('Average: %.0f ± %.0f ms \n', T1avg, T1std);
T2avg = mean(t2v);
T2std = std(t2v);
fprintf('Average: %.0f ± %.0f ms \n', T2avg, T2std);