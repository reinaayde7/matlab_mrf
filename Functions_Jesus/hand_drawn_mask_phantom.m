%Provide the T1 / T2 map to drawn ROIs for curve-fitting
%Provide the number of ROIs
function hand_drawn_mask_phantom(quant_map, num_of_ROIs)

%Compute the size of the image
size_image_x = size(quant_map,1);
size_image_y = size(quant_map,2);

%Draw ROIs around the phantom spheres
ROI_list = [];
for phantom_iteration = 0:(num_of_ROIs-1)
    figure, imagesc(1:size_image_x, 1:size_image_y, quant_map)
    set(gca,'YDir','normal')
    colorbar
    %caxis([40 150])
    h = drawfreehand;
    ROI = createMask(h);
    ROI_list = [ROI_list, ROI];
end

%Compute Mean and Standard Deviation of ROI
mag_result = [];
fprintf('Taking the average and standard deviation of 50 pixels ... \n');
for mask_iteration = 0:(num_of_ROIs-1)
    selected_ROI = ROI_list(:, (1+(mask_iteration*size_image_x)):(size_image_x+(mask_iteration*size_image_x)));
    counter = 1;
    for i = 1:size(selected_ROI,1)
        for j = 1:size(selected_ROI,2)
            if selected_ROI(i,j) == 1 && size(mag_result,2) < 50
                if counter == 25
                    fprintf('The 25th pixel for ROI %d is at x location %d \n', mask_iteration+1, i);
                    fprintf('The 25th pixel for ROI %d is at y location %d \n', mask_iteration+1, j);
                end
                counter = counter + 1;
                mag_result = [mag_result, quant_map(i, j)];
            end
        end
    end
    if size(mag_result) < 50
        fprintf("Bigger sphere needed for iteration %d \n", mask_iteration+1)
    end
    mean_of_ROI = mean(mag_result);
    std_of_ROI = std2(mag_result)
    fprintf('The mean for ROI %d is %.2f \n', mask_iteration+1, mean_of_ROI);
    fprintf('The standard deviation for ROI %d is %.2f \n', mask_iteration+1, std_of_ROI);
    mag_result = [];
end
end