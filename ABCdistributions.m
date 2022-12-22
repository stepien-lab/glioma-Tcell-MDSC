% PROBABILITY DISTRIBUTIONS AND HISTOGRAMS FOR ABC METHOD
clear

load('ABCfornumsim.mat') %ABCfigures.m output


fs = 11; %font size for axis 

% p - value
    pval = 0.05;

    
% Distribution test used
    disttest = 0;   % 0 - Chi-square goodness-of-fit test
                    % 1 - One-sample Kolmogorov-Smirnov test
                    % 2 - Lilliefors test 
% In each case, 1 implies it does NOT come from that distribution
% so, we're happiest when we see 0

% DISTRIBUTIONS TO TEST:
disttype = {'normal';'lognormal'; 'loglogistic';'exponential';'weibull';'kernel';'ev'};

if disttest ==2
    disttype = {'normal';'exponential';'extreme value';'weibull';'lognormal'};
end

numDist = length(disttype);

% HISTOGRAMS TO GRAPH:  
    chooseparam = [0;     % 1  - lambdaC
               0;     % 2  - Cmax
               0;     % 3  - η 
               0;     % 4  - a_T 
               0;     % 5  - s_T 
               0;     % 6 - ρ        
               0;     % 7 - ε_C 
               0;     % 8  - r 
               0;     % 9  - d_T
               0;     % 10 - s_M 
               0;     % 11 - α
               0;     % 12 - q
               0];    % 13 - d_M   

% PARAMETER NAMES:
pn = ["\lambda_C"; 
      "C_{max}";
      "\eta";
      "a_T";
      "s_T";
      "\rho"; 
      "\epsilon_C";
      "r";
      "d_T"; 
      "s_M"; 
      "\alpha";
      "q";
      "d_M"]; 
  

%% Histograms using all accepted parameter values
param_dist = cell(numDist,13);
kstest_dist = nan(numDist,13);
GoF_dist = nan(numDist,13);
lillie_dist = nan(numDist,13);

for i = 1:13
    if chooseparam(i) ==1
        [y,x] = ksdensity(histoP(:,i));
        figure
        plot(x,y);
        xlabel(pn(i))
        set(gca,'FontSize',fs,'YTick',{})
        ylabel('Frequency')
        title(sprintf('Distribution of %s using all accepted parameter values',pn(i)))
    end
    for j=1:numDist
        param_dist{j,i} = fitdist(histoP(:,i),disttype{j});
        if disttest == 0
            GoF_dist(j,i) = chi2gof(histoP(:,i),'CDF',param_dist{j,i},'Alpha',pval,'NBins',10,'EMin',0);
        end
        if disttest == 1
            kstest_dist(j,i) = kstest(histoP(:,i),'CDF',param_dist{j,i},'Alpha',pval);
        end
        if disttest == 2
             if strcmpi(disttype{j},'weibull')
                lillie_dist(j,i) = lillietest(log(histoP(:,i)),'Distribution','extreme value','Alpha',pval);
             elseif strcmpi(disttype{j},'lognormal')
                lillie_dist(j,i) = lillietest(log(histoP(:,i)),'Distribution','normal','Alpha',pval);
             else
                lillie_dist(j,i) = lillietest(histoP(:,i),'Distribution',disttype{j},'Alpha',pval);
             end
        end
    end 
end


%% Histograms using analysis paramater values
analysis_param_dist = cell(numDist,13);
analysis_GoF_dist = nan(numDist,13); 
analysis_kstest_dist = nan(numDist,13); 
analysis_lillie_dist = nan(5,13); 

for i = 1:13
    if chooseparam(i) ==1
        [y,x] = ksdensity(histoanalysisP(:,i));
        figure
        plot(x,y);
        xlabel(pn(i))
        set(gca,'FontSize',fs)
        ylabel('Frequency')
        title(sprintf('Distribution of %s using analysis parameter values',pn(i)))
    end
    for j=1:numDist
        analysis_param_dist{j,i} = fitdist(histoanalysisP(:,i),disttype{j});
        if disttest == 0
            analysis_GoF_dist(j,i) = chi2gof(histoanalysisP(:,i),'CDF',analysis_param_dist{j,i}, 'Alpha',pval,'NBins',10,'EMin',0);
        end
        if disttest == 1
            analysis_kstest_dist(j,i) = kstest(histoanalysisP(:,i),'CDF',analysis_param_dist{j,i},'Alpha',pval);
        end
        if disttest == 2
            if strcmpi(disttype{j},'weibull')
                analysis_lillie_dist(j,i) = lillietest(log(histoanalysisP(:,i)),'Distribution','extreme value','Alpha',pval);
            elseif strcmpi(disttype{j},'lognormal')
                analysis_lillie_dist(j,i) = lillietest(log(histoanalysisP(:,i)),'Distribution','normal','Alpha',pval);
            else
                analysis_lillie_dist(j,i) = lillietest(histoanalysisP(:,i),'Distribution',disttype{j},'Alpha',pval);
            end
        end
    end    
end

if disttest == 0
    GoF_dist
    analysis_GoF_dist
end

if disttest == 1
    kstest_dist
    analysis_kstest_dist
end

if disttest == 2
    lillie_dist
    analysis_lillie_dist
end


%% Checking for uniform distribution of accepted and analysis parameters

accepted_uniform = zeros(2,13);
n1 = length(histoP(:,1));
n2 = length(histoanalysisP(:,1));
uniformv1 = rand(n1,1);
uniformv2 = rand(n2,1);

for i = 1:13
    maxwell1 = max(histoP(:,i));
    minnie1 = min(histoP(:,i));
    testme1 = (histoP(:,i)-minnie1)./(maxwell1-minnie1);
    accepted_uniform(1,i) = kstest2(uniformv1,testme1,'Alpha',pval);
    
    maxwell2 = max(histoanalysisP(:,i));
    minnie2 = min(histoP(:,i));
    testme2 = (histoanalysisP(:,i)-minnie1)./(maxwell2-minnie2);
    accepted_uniform(2,i) = kstest2(uniformv2,testme2,'Alpha',pval);
end

accepted_uniform

% first row for accepted parameters (from ABC)
% second row for best specified % of the accepted parameters

% 1 implies that it rejects the null hypothesis (the null hypothesis is
% that these two vectors are of the same distributions)
% therefore, 0 implies that the parameter has a uniform distribution
% (p-value 0.05)