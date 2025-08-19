function rois = get_rois()

rois = {};
ii = 1;
quitloop = false;

set(gcf,'KeyPressFcn',@add_roi);

while ~quitloop
    pause(0.1);
end

function add_roi(~, event)
    switch event.Key
        case 'c'
            roi = drawcircle('Color','red');
            mask = createMask(roi);
            rois{ii} = mask;
            ii = ii+1;
        case 'f'
            roi = drawfreehand('Color','red');
            mask = createMask(roi);
            rois{ii} = mask;
            ii = ii+1;       
        case 'q'
            quitloop = true;
        otherwise
            disp('not a legit input')
    end
end

end
