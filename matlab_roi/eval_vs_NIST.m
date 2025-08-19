function [eval]=eval_vs_NIST(NIST, t1, t2, rois, b0)

    n = size(rois, 2);
    %-------------------------------------------------------------
    % mean estimates
    eval.means_t1 = zeros(n, 1);
    eval.means_t2 = zeros(n, 1);
    for ii = 1:n
        eval.means_t1(ii) = mean(t1(rois{ii}));
        eval.means_t2(ii) = mean(t2(rois{ii}));
    end
    
    if nargin > 4 %there is b0
        eval.means_b0 = zeros(n, 1);
        for ii = 1:n
            eval.means_b0(ii) = mean(b0(rois{ii}));
        end
    end
    %-------------------------------------------------------------
    % STD estimates
    eval.std_t1 = zeros(n, 1);
    eval.std_t2 = zeros(n, 1);
    
    for ii = 1:n
        eval.std_t1(ii) = std(t1(rois{ii}));
        eval.std_t2(ii) = std(t2(rois{ii}));
    end
    
    if nargin > 4 %there is b0
        eval.std_b0 = zeros(n, 1);
        for ii = 1:n
            eval.std_b0(ii) = std(b0(rois{ii}));    
        end
    end
    %-------------------------------------------------------------
    % COV estimates
    eval.cov_t1 = eval.std_t1./eval.means_t1 *1e2;
    eval.cov_t2 = eval.std_t2./eval.means_t2 *1e2;

    %-------------------------------------------------------------
    % relative errors
    eval.err_t1 = (eval.means_t1 - NIST(1:n, 1))./NIST(1:n, 1) *1e2;
    eval.err_t2 = (eval.means_t2 - NIST(1:n, 2))./NIST(1:n, 2) *1e2;

end