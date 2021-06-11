function [scores, riskdeltas, tse, jstats, ht, pt, basescores] = riskscore( x, riskalpha, runs, varargin)
% Calculates risk-sensitive effectiveness of the given runs.
%
% x: run-by-topic score matrix of a performance measure, such as ERR@20.
%
% riskalpha: is the risk-sensitivity level value from 0 to positive infty. 
%            it can be a series of values, [0:5:15], i.e. for alpha=0,5,10,15
%
% runs: a cell array of string with the name of the runs represented by each row of x.

% vargin can be followings:
%
%       'baseline' is a vector of scalar values from 1 to r where r is the number of rows(runs)
%                  in x, e.g. [1,3,5]. Set baseline to [0] for chosing each run as the
%                  baseline in turn.
%       'measure'
%               'urisk'  TREC 2012 official risk-aware measure
%               'trisk'  Risk-aware score based on t-statistics
%               'prisk'  Down-side loss percent
%               'tprisk' Studentized down-side loss percent
%               'zrisk'  Standardized within topic urisk scores
%               'wrisk'  Wilcoxon signed rank statistics
%               
%               'grisk'  Risk score based on Chi-square statistic
%               'zqrisk' ZRisk based on Chi-squared statistic
%               'tgrisk' TRisk using expected scores as the baseline
%               'georisk' GeoRisk: geometric mean of mean-effectiveness and
%                         ZRisk
%
%       'top'   set to 1 if top runs are used as baselines, and 0 otherwise.
%
%       'jB'    the number of resamples used to produce Jackknife estimates.
%               Set to 0 for risk evaluation without Jackknife estimates.
%
% Outputs:
%
% riskdeltas usage: cell2mat( riskdeltas{1,1,2} ) for {run 1, alpha value 1, baseline run 2}
%
% Author: Bekir Taner Dincer.

% Check inputs
narginchk(2,Inf);
if nargin < 3
    runs = [];
end
if nargin < 2
    riskalpha = 0;
end

% number of runs
r = size(x,1);

[~, IX] = sort(mean(x,2), 1, 'descend');

% number of topics
c = size(x,2);

if isempty(runs) || ~iscellstr(runs)
    runs = 1:r;
end

% Get optional parameters
% jB: number of re-samples B for delete-d jackknife estimates
paramNames = {'baseline', 'measure', 'top', 'jB'};
paramDflts = {0, 'urisk', 0, 0};
[baseline, measure, topkeval, jB] = internal.stats.parseArgs(paramNames, paramDflts, varargin{:});

if isempty(measure)
    error('!!! Risk measure option cannot be set to none !!!');
end

if ~isvector(baseline) && (baseline < 0 || baseline > r)
    error('!!! Baseline run number must in between %d and %d !!!', 0, r);
end

% the number of baseline runs in this evaluation
b = size(baseline,2);

if b == 1 && baseline == 0
    baseline = 1:r;
    b = r;
end

% Risk-aware performance scores
scores = zeros(r,size(riskalpha,2),b);

% paired t-test result h and p-value and tstats.sd
tse = zeros(r,size(riskalpha,2),b); % power(tstats.sd,2)/50
ht = zeros(r,size(riskalpha,2),b);
pt = zeros(r,size(riskalpha,2),b);

% Jackknife stats 
jdsize = 2; % one estimation for d=1 and one for d=c/2 (half sample)
jstats.bias = zeros(r,size(riskalpha,2),b,jdsize);
jstats.var = zeros(r,size(riskalpha,2),b,jdsize);
jstats.thetas = zeros(r,size(riskalpha,2),b,jdsize);
jstats.bcthetas = zeros(r,size(riskalpha,2),b,jdsize);
jstats.histdata = cell(r,size(riskalpha,2),b,jdsize);


% Structure of riskdeltas: 
% for (each_run, each_alpha, each_baseline)
%    1 : r x c matrix of risk-deltas i.e. (1+alpha) * (r - b)
%    2 : the name of the run r - the name of the baseline b - alpha value
%  the following is produced only for urisk
%    3 : 5 x 1 matrix of the values for 5 quality metrics
%    where 
%    row 0: mean score based on the effectiveness score.
%    row 1: risk/reward : total risk / total reward
%    row 2: loss/win  : # 
%    row 3: loss : # of topics where r < b
%    row 4: win  : # of topics where r > b
%    row 5: Loss > 20% : # of topics where r/b < 0.80.
%
% usage: 
% cell2mat( riskdeltas{run_i,alpha_i,baseline_i,1} ) 
% The risk-deltas values of run_i at alpha_i using baseline_i across c topics

if (strcmp(measure, 'urisk') == 1)
    riskdeltas = cell(r, size(riskalpha,2), b, 3);
else
    riskdeltas = cell(r, size(riskalpha,2), b, 2);
end

if ((strcmp(measure,'grisk') == 1) || (strcmp(measure,'zgrisk') == 1) || (strcmp(measure,'tgrisk') == 1) || (strcmp(measure,'georisk') == 1))

    % Baseline scores for GRisk and ZGRisk
    basescores = zeros(r,c);
    
    % Calculate sums
    S = sum(x,2);
    T = sum(x,1);
    N = sum(S);
    
    % for each run
    for i=1:r
        % baseline scores for run i
        for j=1:c
            basescores(i,j) = (S(i) * T(j)) / N;
        end
        
        switch(measure)
            case 'grisk'
                [scores(i,1,1),deltas] = grisk( x(i,:), basescores(i,:));
                riskdeltas{i,1,1,1} = num2cell(deltas);
                
                riskdeltas{i,1,1,2} = [runs{i} '-GRisk'];
            case 'zgrisk'
                for j=1:size(riskalpha,2)
                    [scores(i,j,1), deltas] = zgrisk( x(i,:), basescores(i,:), riskalpha(j));
                    riskdeltas{i,j,1,1} = num2cell(deltas);

                    riskdeltas{i,j,1,2} = [runs{i} '-ZGRisk' ' a=' num2str(riskalpha(j))];
                end
            case 'tgrisk'
                for j=1:size(riskalpha,2)
                    [scores(i,j,1), deltas, tse(i,j,1),...
                        ht(i,j,1),pt(i,j,1)] = tgrisk( x(i,:), basescores(i,:), riskalpha(j));
                    
                    riskdeltas{i,j,1,1} = num2cell(deltas);
                    
                    riskdeltas{i,j,1,2} = [runs{i} '-TGRisk' ' a=' num2str(riskalpha(j))];
                end
            case 'georisk'
                for j=1:size(riskalpha,2)
                    [scores(i,j,1), deltas] = georisk( x(i,:), basescores(i,:), riskalpha(j));
                    riskdeltas{i,j,1,1} = num2cell(deltas);

                    riskdeltas{i,j,1,2} = [runs{i} '-GeoRisk' ' a=' num2str(riskalpha(j))];
                end
                
        end % switch
        
    end %for each run
    
else
for baserun=1:b
    if topkeval
       bscores = x(IX(baseline(baserun)),:);
    else
       bscores = x(baseline(baserun),:);
    end
    
    switch(measure)
        case 'urisk'
            % TREC2012 official measure
            for i=1:r
                for j=1:size(riskalpha,2)
                    
                    [scores(i,j,baserun), deltas, quality] = urisk( x(i,:), bscores, riskalpha(j));
                    
                    riskdeltas{i,j,baserun,1} = num2cell(deltas);
                    
                    if topkeval 
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{IX(baseline(baserun))} ' a=' num2str(riskalpha(j))];
                    else
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{baseline(baserun)} ' a=' num2str(riskalpha(j))];
                    end
                    
                    qualityLabels = {'NDCG@10'; 'Risk/Reward'; 'Loss/Win'; 'Loss'; 'Win'; 'Loss > 20%'};
                    riskdeltas{i,j,baserun,3} = [qualityLabels, num2cell(quality)];

                    if jB > 0
                        %Delete-1 Jackknife estimates
                        [jbias, jvar, thetas, bcthetas, histdata, B] = jackknifed(deltas',1,jB);
                        jstats.k(1) = 1;
                        jstats.B(1) = B;
                        jstats.bias(i,j,baserun,1) = jbias; 
                        jstats.var(i,j,baserun,1) = jvar;
                        jstats.thetas(i,j,baserun,1) = thetas;
                        jstats.bcthetas(i,j,baserun,1) = bcthetas;
                        jstats.histdata{i,j,baserun,1} = num2cell(histdata);
                    
                        %Delete-k Jackknife estimates
                        k= c/2;
                        [jbias, jvar, thetas, bcthetas, histdata, B] = jackknifed(deltas',k,jB);
                        jstats.k(2) = k;
                        jstats.B(2) = B;
                        jstats.bias(i,j,baserun,2) = jbias; 
                        jstats.var(i,j,baserun,2) = jvar;
                        jstats.thetas(i,j,baserun,2) = thetas;
                        jstats.bcthetas(i,j,baserun,2) = bcthetas;
                        jstats.histdata{i,j,baserun,2} = num2cell(histdata);
                    end
                end
            end
        case 'trisk'
            % Risk measure based on t-statistic
            for i=1:r
                for j=1:size(riskalpha,2)
                    [scores(i,j,baserun), deltas, tse(i,j,baserun),...
                        ht(i,j,baserun),pt(i,j,baserun)] = trisk( x(i,:), bscores, riskalpha(j));
                    
                    riskdeltas{i,j,baserun,1} = num2cell(deltas);
                    
                    if topkeval 
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{IX(baseline(baserun))} ' a=' num2str(riskalpha(j))];
                    else
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{baseline(baserun)} ' a=' num2str(riskalpha(j))];
                    end
                end
            end
        case 'zrisk'
            % standardized within topic urisk scores
            for j=1:size(riskalpha,2)
                for i=1:r
                    [scores(i,j,baserun),deltas] = urisk( x(i,:), bscores, riskalpha(j));
                    riskdeltas{i,j,baserun,1} = num2cell(deltas);
                    
                    if topkeval 
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{IX(baseline(baserun))}];
                    else
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{baseline(baserun)}];
                    end
                end
                
                % Standardize deltas within topics
                deltas = zeros(r,c);
                for i=1:r
                    deltas(i,:) = cell2mat(riskdeltas{i,j,baserun,1})';
                end
                z = zscore(deltas);
                for i=1:r
                    riskdeltas{i,j,baserun,1} = num2cell(z(i,:));
                end
            end
        case 'prisk'
            % down-side risk percent
            for i=1:r
                for j=1:size(riskalpha,2)
                    [scores(i,j,baserun),deltas] = prisk( x(i,:), bscores, riskalpha(j));
                    riskdeltas{i,j,baserun,1} = num2cell(deltas);
                    
                    if topkeval 
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{IX(baseline(baserun))}];
                    else
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{baseline(baserun)}];
                    end
                end
            end
        case 'tprisk'
            % studentized down-side risk percent
            for i=1:r
                for j=1:size(riskalpha,2)
                    [scores(i,j,baserun),deltas] = tprisk( x(i,:), bscores, riskalpha(j));
                    riskdeltas{i,j,baserun,1} = num2cell(deltas);
                    
                    if topkeval 
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{IX(baseline(baserun))}];
                    else
                        riskdeltas{i,j,baserun,2} = [runs{i} '-' runs{baseline(baserun)}];
                    end
                end
            end
            
    end % switch
    
end % for each baseline

end %outher if

% Jackknife stats: deleted
% for baserun = 1:b
%     if topkeval
%        bscores = x(IX(baseline(baserun)),:);
%     else
%        bscores = x(baseline(baserun),:);
%     end
%     
%     for i=1:r
%         for j=1:size(riskalpha,2)
%             [~, udeltas] = urisk( x(i,:), bscores, riskalpha(j));
%             %Delete-1 Jackknife estimates
%             [jstats.bias(i,j,baserun,1), jstats.var(i,j,baserun,1),... 
%                 jstats.thetas(i,j,baserun,1), jstats.bcthetas(i,j,baserun,1),...
%                 jstats.histdata(i,j,baserun,1)] = jackknifed(udeltas,1,B);
%             %Delete-k Jackknife estimates
%             k = c/2;
%             [jstats.bias(i,j,baserun,2), jstats.var(i,j,baserun,2),... 
%                 jstats.thetas(i,j,baserun,2), jstats.bcthetas(i,j,baserun,2),...
%                 jstats.histdata(i,j,baserun,2)] = jackknifed(udeltas,k,B);
%         end
%     end
% end

% meanalpharisk = zeros(r,size(riskalpha,2));
% for i=1:r
%     for j=1:size(riskalpha,2)
%         scores_ = scores(i,j,:);
%         scores_(scores_(:)==0) = [];
%         if isempty(scores_)
%             meanalpharisk(i,j) = 0;
%         else
%             meanalpharisk(i,j) = mean(scores_);
%         end
%     end
% end
% 
% meanrisk = zeros(r,1);
% for i=1:r
%     meanrisk(i) = mean(meanalpharisk(i,:));
% end

end

% calculates the URisk score for the given model scores "rscores" using
% baseline scores "bscores" at alpha = riskalpha.
function [uscore, deltas, quality]=urisk(rscores, bscores, riskalpha)
    
    c = size(rscores,2); % number of topics
    deltas = zeros(c,1); % delta value for each topic
    risk = 0; % sum of deltas where r < b
    reward = 0; % sum of deltas where r > b
    loss = 0; % # of topics where r < b
    win = 0; % # of topics where r > b
    lossgeq20 = 0;% # of topics where r/b < 0.80
    
    for j=1:c
        sdiff = rscores(j) - bscores(j);
        if sdiff >= 0
            deltas(j) = sdiff;
            reward = reward + deltas(j);
            win = win + 1;
        else
            deltas(j) = ( (1 + riskalpha) * sdiff);
            risk = risk + abs(deltas(j));
            loss = loss + 1;
        end
        
        if (rscores(j) / bscores(j)) < 0.80
            lossgeq20 = lossgeq20 + 1;
        end
    end
    
    uscore = mean(deltas);
    quality = [mean(rscores); risk/reward; loss/win; loss; win; lossgeq20];
end

function [tscore,deltas,tse,h,p]=trisk(rscores, bscores, riskalpha)

    % Standardized deltas by paired standard deviation
    c = size(rscores,2); % number of topics
    deltas = zeros(c,1);
    for j=1:c
        sdiff = rscores(j) - bscores(j);
        if sdiff >= 0
            deltas(j) = sdiff;
        else
            deltas(j) = ( (1 + riskalpha) * sdiff);
        end
    end
    
    [h,p,~,stats] = ttest(deltas, 0, 0.05, 'both');
    if isnan(h) 
        h=0; 
    end
    if isnan(p) 
        p=1; 
    end
    
    if isnan(stats.tstat)
        tscore = 0;
        tse = 0;
    else
        %m = mean(deltas);
        tscore = stats.tstat; %i.e., tscore = sqrt(c) * m / stats.sd;
        tse = stats.sd/sqrt(c);
        % triskdeltas = delete-1 jackknife histdata
        for j=1:c
            deltas(j) = deltas(j) / stats.sd; % d_i / sd / sqrt(c), where c = 1
        end
        
    end
    
end

function [gscore, deltas]=grisk(rscores, bscores)
    % number of topics
    c = size(rscores,2);
    deltas = zeros(c,1);
    for j=1:c
        sdiff = power(rscores(j) - bscores(j),2);
        if bscores(j) > 0
            deltas(j) = sdiff / bscores(j);
        else
            deltas(j) = sdiff / 0.0001;
        end
    end
    gscore = sum(deltas);
end

function [zgscore, deltas]=zgrisk(rscores, bscores, riskalpha)

    % number of topics
    c = size(rscores,2);
    deltas = zeros(c,1);
    for j=1:c
        if bscores(j) > 0
            zg = (rscores(j) - bscores(j)) / sqrt(bscores(j));
        else
            zg = (rscores(j) - bscores(j)) / sqrt(0.0001);
        end
        if zg >= 0
            deltas(j) = zg;
        else
            deltas(j) = ( (1 + riskalpha) * zg);
        end
    end
    
    zgscore = sum(deltas);
    
end

function [geoscore, deltas]=georisk(rscores, bscores, riskalpha)

    % number of topics
    c = size(rscores,2);
    deltas = zeros(c,1);
    for j=1:c
        if bscores(j) > 0
            zg = (rscores(j) - bscores(j)) / sqrt(bscores(j));
        else
            zg = (rscores(j) - bscores(j)) / sqrt(0.0001);
        end
        if zg >= 0
            deltas(j) = zg;
        else
            deltas(j) = ( (1 + riskalpha) * zg);
        end
    end
    
    zgscore = sum(deltas);
    
    geoscore = sqrt(mean(rscores) * normcdf(zgscore/c, 0, 1));
end

function [tscore,deltas,tse,h,p]=tgrisk(rscores, bscores, riskalpha)

    % Standardized deltas by paired standard deviation
    c = size(rscores,2); % number of topics
    deltas = zeros(c,1);
    
    for j=1:c
        sdiff = rscores(j) - bscores(j);
        if sdiff >= 0
            deltas(j) = sdiff;
        else
            deltas(j) = ( (1 + riskalpha) * sdiff);
        end
    end
    
    [h,p,~,stats] = ttest(deltas, 0, 0.05, 'left');
    if isnan(h) 
        h=0; 
    end
    if isnan(p) 
        p=1; 
    end
    
    if isnan(stats.tstat)
        tscore = 0;
        tse = 0;
    else
        %m = mean(deltas);
        tscore = stats.tstat; %i.e., tscore = sqrt(c) * m / stats.sd;
        tse = stats.sd/sqrt(c);
        % triskdeltas = delete-1 jackknife histdata
        for j=1:c
            deltas(j) = deltas(j) / stats.sd; % d_i / sd / sqrt(c), where c = 1
        end
    end
end

function [pscore, deltas]=prisk(rscores, bscores, alpha)

    % number of topics
    c = size(rscores,2);
    deltas = zeros(c,1);
    for j=1:c
        sdiff = rscores(j) - bscores(j);
        if sdiff >= 0
            deltas(j) = sdiff;
        else
            deltas(j) = ( (1 + alpha) * sdiff);
        end
    end
    
    
    deltas(deltas(:) == 0) = [];

    if mean(deltas) == 0 || isempty(deltas)
        pscore = 0;
    else
        deltaplus = 0;
        deltaminus = 0;
        for j=1:size(deltas,1)
            if deltas(j) < 0
                deltaminus = deltaminus + abs(deltas(j));
            else
                deltaplus = deltaplus + deltas(j);
            end
        end

        pscore = deltaminus / (deltaminus + deltaplus);
    end
    
end

function [tpscore, deltas]=tprisk(rscores, bscores, alpha)

    % number of topics
    c = size(rscores,2);
    deltas = zeros(c,1);
    for j=1:c
        sdiff = rscores(j) - bscores(j);
        if sdiff >= 0
            deltas(j) = sdiff;
        else
            deltas(j) = ( (1 + alpha) * sdiff);
        end
    end
    
    deltas_= zeros(c,1);
    deltas_(:) = deltas(:);
    deltas_(deltas(:) == 0) = [];
    
    if mean(deltas_) == 0 || isempty(deltas_)
        tpscore = 0;
    else
        deltastd = std(deltas_);
        deltaplus = 0;
        deltaminus = 0;
        for j=1:size(deltas_,1)
            if deltas_(j) < 0
                deltaminus = deltaminus + (abs(deltas_(j)) / deltastd);
            else
                deltaplus = deltaplus + (deltas_(j) / deltastd);
            end
            deltas(j) = sqrt(c) * deltas(j) / deltastd;
        end

        tpscore = deltaminus / (deltaminus + deltaplus);
    end
    
end
