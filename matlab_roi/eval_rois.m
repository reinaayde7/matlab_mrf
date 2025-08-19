function eval_rois(varargin)

set(gcf,'KeyPressFcn',@eval_roi);

function eval_roi(~, event)
    switch event.Key
        case 'c'
            roi = drawcircle('Color','red');
            mask = createMask(roi);
            printvals(mask);
        case 'f'
            roi = drawfreehand('Color','red');
            mask = createMask(roi);
            printvals(mask)   
        otherwise
            disp('not a legit input')
    end
end

function printvals(mask)
    for i=1:numel(varargin)
        arr = varargin{i};
        fprintf('%.2f+-%.2f\t', mean(arr(mask)), std(arr(mask)));
    end
    fprintf('\n')
end

end