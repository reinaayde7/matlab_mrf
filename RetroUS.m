function [usraw, smp_mask] = RetroUS(raw, order, plotpattern)
%Function to retrospectively k-t undersample raw data
%to a given order. Plotpattern is true or false variable
%to plot the current Z k-t undersampling pattern
%Raw data should be 4D, being second index the slice index
%And third index the TimePoint index. Check numRows_orig is the 
%# of slices and numCols_orig the number of timepoints
%Jesus Fajardo 10/18/24

%Take raw data dimensions 
%Define a ones matrix and then multiply the pattern on the edges

if order == 2,
    basePattern = [1 0;
                   0 1];
end
if order == 3,
    basePattern = [1 0 0;
                   0 1 0;
                   0 0 1];
end

% Define the desired matrix size
numRows_orig = size(raw,2); % Number of rows in the final matrix
numCols_orig = size(raw,3); % Number of columns in the final matrix
FullySampled = ones(numRows_orig,numCols_orig);
numRows = round(numRows_orig * 0.3);
% Calculate how many times to repeat the base pattern to cover the desired size
repeatRows = ceil(numRows / size(basePattern, 1)); % Number of times to repeat in row dimension
repeatCols = ceil(numCols_orig / size(basePattern, 2)); % Number of times to repeat in column dimension
% Generate the larger repetitive matrix
repetitiveMatrix = repmat(basePattern, repeatRows, repeatCols);
% Trim the matrix to the exact desired size
repetitiveMatrix = repetitiveMatrix(1:numRows, 1:numCols_orig);

for i = 1:numRows,
        for j = 1:numCols_orig
            FullySampled(i,j) = repetitiveMatrix(i,j);
        end
end
size(repetitiveMatrix)
for i = ceil(size(FullySampled,1)-(numRows-1)):size(FullySampled),
    i;
        for j = 1:numCols_orig
            FullySampled(i,j) = repetitiveMatrix( ceil((size(FullySampled,1))-i+1) , j);
        end
end

% Display the matrix
%disp(repetitiveMatrix);

if plotpattern
    figure('DefaultAxesFontSize',16)
    subplot()
    imagesc(FullySampled(:,1:20));
    colormap('hot');
    set(gca, 'XTick', [], 'YTick', []);
    axis off;
    daspect([1 1 1])
    title('')
    set(gca, 'position', [0 0 1 1], 'units', 'normalized')
end
disp('undersampling raw data...') 
for i = 1:numRows_orig,
        for j = 1:numCols_orig
        usraw(:,i,j,:) = raw(:,i,j,:)*FullySampled(i,j);
        end
end
smp_mask = FullySampled;
end
