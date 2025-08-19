function r = paramTable(t1series,t2series,varargin)
% Create table of all possible T1 and T2 combinations
% Yun Jiang -- switch from if to Switch...case for more clear presentation

switch (length(varargin))
    case 0
        cnt=0;
        for it1 = 1:length(t1series)
            for it2 = 1:length(t2series)
                if (t2series(it2)<=t1series(it1))
                    cnt=cnt+1;
                end
            end
        end
        r = zeros(cnt,2);
        
        cnt=0;
        for it1 = 1:length(t1series)
            for it2 = 1:length(t2series)
                if (t2series(it2)<=t1series(it1))
                    cnt=cnt+1;
                    r(cnt,1)=t1series(it1);
                    r(cnt,2)=t2series(it2);
                end
            end
        end
        
        
        
    case 1
        dfseries = varargin{1};
        cnt=0;
        for it1 = 1:length(t1series)
            for it2 = 1:length(t2series)
                for idf = 1:length(dfseries)
                    if (t2series(it2)<=t1series(it1))
                        cnt=cnt+1;
                    end
                end
            end
        end
        r = zeros(cnt,3);
        
        cnt = 0;
        for it1 = 1:length(t1series)
            for it2 = 1:length(t2series)
                for idf = 1:length(dfseries)
                    if (t2series(it2)<=t1series(it1))
                        cnt=cnt+1;
                        r(cnt,1)=t1series(it1);
                        r(cnt,2)=t2series(it2);
                        r(cnt,3)=dfseries(idf);
                    end
                end
            end
        end
        
    case 2
        dfseries = varargin{1};
        b1series = varargin{2};
        cnt=0;
        for it1 = 1:length(t1series)
            for it2 = 1:length(t2series)
                for idf = 1:length(dfseries)
                    for idb1 = 1:length(b1series)
                        if (t2series(it2)<=t1series(it1))
                            cnt=cnt+1;
                        end
                    end
                end
            end
        end
        r = zeros(cnt,4);
        
        cnt = 0;
        for it1 = 1:length(t1series)
            for it2 = 1:length(t2series)
                for idf = 1:length(dfseries)
                    for idb1 = 1:length(b1series);
                        if (t2series(it2)<=t1series(it1))
                            cnt=cnt+1;
                            r(cnt,1)=t1series(it1);
                            r(cnt,2)=t2series(it2);
                            r(cnt,3)=dfseries(idf);
                            r(cnt,4) = b1series(idb1);
                        end
                    end
                end
            end
        end
        
    case 3
        
        dfseries = varargin{1};
        b1series = varargin{2};
        diffseries = varargin{3};
        cnt=0;
        for it1 = 1:length(t1series)
            for it2 = 1:length(t2series)
                for idf = 1:length(dfseries)
                    for idb1 = 1:length(b1series)
                        for idiff = 1:length(diffseries)
                            if (t2series(it2)<=t1series(it1))
                                cnt=cnt+1;
                            end
                        end
                    end
                end
            end
        end
        
        r = zeros(cnt,5);
        
        cnt = 0;
        for it1 = 1:length(t1series)
            for it2 = 1:length(t2series)
                for idf = 1:length(dfseries)
                    for idb1 = 1:length(b1series);
                        if (t2series(it2)<=t1series(it1))
                            for idiff = 1:length(diffseries)
                                
                                cnt=cnt+1;
                                r(cnt,1)=t1series(it1);
                                r(cnt,2)=t2series(it2);
                                r(cnt,3)=dfseries(idf);
                                r(cnt,4) = b1series(idb1);
                                r(cnt,5) = diffseries(idiff);
                            end
                        end
                    end
                end
            end
        end
end