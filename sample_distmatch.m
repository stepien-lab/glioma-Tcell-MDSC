%clear variables global;
%clc;

%load('ABCfornumsim.mat') %ABCfigures.m output (obtained after using average data for ABCrejection.m)


% Plot modes
    plotmode = 3;   % 0 - average of data points at a certain day
                    % 1 - all data points
                    % 2 - uses both 0 and 1
                    % 3 - errorbar where center is mean and outer bounds
                    % are min and max data points
                    
    datamode = 1;   % 0 - use accepted values for std and mean
                    % 1 - use analysis values for std and mean
                    
    width = 1.5; % line thickness for mean plot
    tickwidth = 0.8; % line thickness otherwise
    scattersize = 25; % scatter plot size
    fs = 17; % font size for plots

num_param = 13; % number of parameters
num_hold = 13; 
    

if datamode==0
    param_sort_hold = histoP; % matrix of accepted parameters that have been sorted (by error) and held onto after ABC rejection
end

if datamode==1
    param_sort_hold = histoanalysisP; % matrix of top accepted parameters that have been sorted (by error) and held onto after ABC rejection
end


bound = [0            0.5;       % 1  - lambdaC
        1e6          5e7;        % 2  - Cmax
        0            1e-6;       % 3  - η 
        5e1          5e6;        % 4  - a_T 
        1e2          1e7;        % 5  - s_T 
        0            0.5;        % 6  - ρ       
        1            100;        % 7  - ε_C 
        0            1e-4;       % 8  - r 
        0            0.75;       % 9  - d_T
        0            0.1;        % 10 - s_M 
        1e7          5e8;        % 11 - α
        1e9          1e11;       % 12 - q
        0            0.5];       % 13 - d_M
% matrix that has lower bounds of parameters in the first column and upper bounds in the second column

%% Fit the data to different probabiltiy distributions

dist_type = {'Normal';'Lognormal';'Gamma';'Exponential';'Weibull';...
    'Logistic';'Uniform'};
%%% didn't use these distributions:
%%% 'beta';'birnbaumsaunders';'burr';'negative binomial';'extreme value';'kernel';
%%% 'generalized extreme value';'generalized pareto';'inversegaussian';
%%% 'nakagami';'loglogistic';'poisson';'rayleigh';'rician';'tlocationscale';
num_dist = length(dist_type);

dist_param = cell(num_dist,num_param);

% all distributions but uniform
for i=1:num_dist-1
    for j=1:num_param
        [j]
        dist_param{i,j} = fitdist(param_sort_hold(:,j),dist_type{i});
    end
end

% uniform
for j=1:num_param
    dist_param{num_dist,j}.Lower = bound(j,1);
    dist_param{num_dist,j}.Upper = bound(j,2);
end

%% Create synthetic data based on the fitted distributions
rng(100,'twister')

dist_synth = cell(num_dist,num_param);

for i=1:num_dist
    for j=1:num_param
        if strcmp(dist_type{i},'Normal')==1 || strcmp(dist_type{i},'Lognormal')==1 ...
                || strcmp(dist_type{i},'Logistic')==1
            dist_synth{i,j} = random(dist_type{i},dist_param{i,j}.mu,...
                dist_param{i,j}.sigma,num_hold,1);

         elseif strcmp(dist_type{i},'Gamma')==1
            dist_synth{i,j} = random(dist_type{i},dist_param{i,j}.a,...
                dist_param{i,j}.b,num_hold,1);

        elseif strcmp(dist_type{i},'Exponential')==1
            dist_synth{i,j} = random(dist_type{i},dist_param{i,j}.mu,...
                num_hold,1);

        elseif strcmp(dist_type{i},'Weibull')==1
            dist_synth{i,j} = random(dist_type{i},dist_param{i,j}.A,...
                dist_param{i,j}.B,num_hold,1);

        elseif strcmp(dist_type{i},'Uniform')==1
            dist_synth{i,j} = dist_param{i,j}.Lower ...
                + (dist_param{i,j}.Upper - dist_param{i,j}.Lower).*rand(num_hold,1);
        end
    end
end


%% Calculate difference between synthetic and data probability distributions 
%%% using Weisserstein metric / Earth mover's distance

wsd1 = zeros(num_dist,num_param);
wsd2 = zeros(num_dist,num_param);

for i=1:num_dist
    for j=1:num_param
        wsd1(i,j) = ws_distance(dist_synth{i,j},param_sort_hold(:,j),1);
        wsd2(i,j) = ws_distance(dist_synth{i,j},param_sort_hold(:,j),2);
    end
end

bestfitdist = cell(1,num_param);
bestfitdist_param = cell(1,num_param);

%% Determine best fitting distribution

numsamp = 10000;
rng(numsamp,'twister')

% collect parameter samples from dist
randomsamples = zeros(13,numsamp);
sz = size(randomsamples(1,:));

for j=1:num_param
    ind = min(wsd1(:,j))==wsd1(:,j);
    bestfitdist{j} = dist_type{ind};
    bestfitdist_param(j) = dist_param(ind,j);
    if ind(7) == 1 % ie param is uniform
        randomsamples(j,:) = (bound(j,2) - bound(j,1)).*rand(numsamp,1) + bound(j,1);
    else
        randomsamples(j,:) = random(bestfitdist_param{j}, sz);
    end
    for i = 1:numsamp
        while randomsamples(j,i)<0
            randomsamples(j,i) = random(bestfitdist_param{j}); 
        end
    end
end
disp(bestfitdist);

%% Plot histograms

tiledlayout(4,4,'TileSpacing','compact','Padding','compact')
for i=1:num_param
    nexttile(i)

    if strcmp(bestfitdist{i},'Uniform')==0
        h = histfit(param_sort_hold(:,i),[],bestfitdist{i});
        h(1).FaceColor = 'none';
        h(2).Color = 'k';
        box on
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',round(yt/num_hold,2));
    else
        hold on
        histogram(param_sort_hold(:,i),'Normalization','probability',...
            'BinMethod','sqrt','FaceColor','none');
        hh = get(gca,'YLim');
        plot(linspace(bound(i,1),bound(i,2),100),hh(2)/2*ones(1,100),'k',...
            'LineWidth',2.5)
        box on
        hold off
    end

    xlabel(pn{i},'Interpreter','latex')
    ylabel('Percentage','Interpreter','latex')
    xlim([0,bound(i,2)])
end

%% Determine mean and std of output

% solve ODE for each of parameter sample 
% collect numerical data points at time points
timerange =  linspace(0,40,40*24); 
tsize = length(timerange);
initialcondition = [35000;100;0]; 
Prandom = zeros(13,1);
TUMORnumsim = zeros(tsize,numsamp); 
MDSCnumsim = zeros(tsize,numsamp);
TCELLnumsim = zeros(tsize,numsamp);


for i=1:numsamp
    [i]
    
    [T1,X1] = ode23s(@(t,x) GBMFunc(t,x,randomsamples(:,i)), timerange, initialcondition);
    
    TUMORnumsim(:,i) = X1(:,1);

    TCELLnumsim(:,i) = X1(:,2); 
    
    MDSCnumsim(:,i) = X1(:,3);
end


% take std and mean of data points at time points

meanrandom = zeros(tsize,3);
clear mean

meanrandom(:,1) = mean(TUMORnumsim,2);
meanrandom(:,2) = mean(TCELLnumsim,2);
meanrandom(:,3) = mean(MDSCnumsim,2);


stdrandom = zeros(tsize,3);

stdrandom(:,1) = std(TUMORnumsim,1,2); %population std
stdrandom(:,2) = std(TCELLnumsim,1,2);
stdrandom(:,3) = std(MDSCnumsim,1,2);

plusstdrandom = meanrandom + stdrandom; % tsize x 3 (tsize times, 3 cell pops)
plusstdrandomquarter = meanrandom + stdrandom/4; 
plusstdrandomhalf = meanrandom + stdrandom/2;
plusstdrandom3quarter = meanrandom + 3*stdrandom/4; 
minusstdrandom = meanrandom - stdrandom; 
minusstdrandomquarter = meanrandom - stdrandom/4; 
minusstdrandomhalf = meanrandom - stdrandom/2; 
minusstdrandom3quarter = meanrandom - 3*stdrandom/4; 

%% Restructure data for plot mode 3

numoftimedata = length(Averagetumordata(:,1));
 
maxtumordata = zeros(numoftimedata,1);
mintumordata = zeros(numoftimedata,1);
    
maxTcelldata = zeros(numoftimedata,1);
minTcelldata = zeros(numoftimedata,1);
    
maxMDSCdata = zeros(numoftimedata,1);
minMDSCdata = zeros(numoftimedata,1);
  
  
if plotmode == 3
    for i = 1:numoftimedata
        elm = Averagetumordata(i,1);
        
        
        maxtumordata(i) = max(tumordata(:,2).*(tumordata(:,1) == elm));
        mintumordata(i) = min(nonzeros(tumordata(:,2).*(tumordata(:,1) == elm)));
    
        maxTcelldata(i) = max(Tcelldata(:,2).*(Tcelldata(:,1) == elm));
        minTcelldata(i) = maxTcelldata(i);
        
        if elm == 24
            minTcelldata(i) = min(Tcelldata(2:4,2)); 
        end
        
    
        maxMDSCdata(i) = max(MDSCdata(:,2).*(MDSCdata(:,1) == elm));
        minMDSCdata(i) = min(nonzeros(MDSCdata(:,2).*(MDSCdata(:,1) == elm)));
        
    end
    
    
    postumor = maxtumordata - Averagetumordata(:,2);
    negtumor = Averagetumordata(:,2) - mintumordata;
    
    posTcell = nonzeros(maxTcelldata) - AverageTcelldata(:,2);
    negTcell = AverageTcelldata(:,2) - nonzeros(minTcelldata);
    
    posMDSC = maxMDSCdata - AverageMDSCdata(:,2);
    negMDSC = AverageMDSCdata(:,2) - minMDSCdata;
end

%% Plot std and mean with shading

x2 = [timerange, fliplr(timerange)];
    
inBetweentumor = [plusstdrandom(:,1)', fliplr(minusstdrandom(:,1)')];
inBetweenTcell = [plusstdrandom(:,2)', fliplr(minusstdrandom(:,2)')];
inBetweenMDSC = [plusstdrandom(:,3)', fliplr(minusstdrandom(:,3)')];

inBetweentumorquarter = [plusstdrandomquarter(:,1)', fliplr(minusstdrandomquarter(:,1)')];
inBetweenTcellquarter = [plusstdrandomquarter(:,2)', fliplr(minusstdrandomquarter(:,2)')];
inBetweenMDSCquarter = [plusstdrandomquarter(:,3)', fliplr(minusstdrandomquarter(:,3)')];

inBetweentumorhalf = [plusstdrandomhalf(:,1)', fliplr(minusstdrandomhalf(:,1)')];
inBetweenTcellhalf = [plusstdrandomhalf(:,2)', fliplr(minusstdrandomhalf(:,2)')];
inBetweenMDSChalf = [plusstdrandomhalf(:,3)', fliplr(minusstdrandomhalf(:,3)')];

inBetweentumor3quarter = [plusstdrandom3quarter(:,1)', fliplr(minusstdrandom3quarter(:,1)')];
inBetweenTcell3quarter = [plusstdrandom3quarter(:,2)', fliplr(minusstdrandom3quarter(:,2)')];
inBetweenMDSC3quarter = [plusstdrandom3quarter(:,3)', fliplr(minusstdrandom3quarter(:,3)')];


figure
tiledlayout(1,3,'TileSpacing','compact');
nexttile(1)
fill(x2, inBetweentumor, [0.25 0.25 0.25]);
hold on;
fill(x2, inBetweentumor3quarter, [0.5 0.5 0.5]);
fill(x2, inBetweentumorhalf, [0.75 0.75 0.75]);
fill(x2, inBetweentumorquarter, [1 1 1]);
plot(timerange, meanrandom(:,1),'--','Color','r','LineWidth',width);
if plotmode == 0
    scatter(Averagetumordata(:,1), Averagetumordata(:,2),scattersize,'k','filled')
end
if plotmode == 1
    scatter(tumordata(:,1), tumordata(:,2),scattersize,'k','filled')
end
if plotmode == 2
    scatter(tumordata(:,1), tumordata(:,2),scattersize,'d','filled',"MarkerFaceColor","#77AC30")
    scatter(Averagetumordata(:,1), Averagetumordata(:,2),scattersize,'k','filled')
end
if plotmode == 3
    scatter(tumordata(:,1), tumordata(:,2),scattersize,'d','filled',"MarkerFaceColor",'b')
    errorbar(Averagetumordata(:,1), Averagetumordata(:,2), negtumor, postumor,'.','color','b')
end

hold off 
xlabel('t (days)')
ylabel('C (cells)')
set(gca,'FontSize',fs,'LineWidth',tickwidth)


%---T-Cell--Plot-----------------------------------------------------------

nexttile(2)
fill(x2, inBetweenTcell, [0.25 0.25 0.25]);
hold on;
fill(x2, inBetweenTcell3quarter, [0.5 0.5 0.5]);
fill(x2, inBetweenTcellhalf, [0.75 0.75 0.75]);
fill(x2, inBetweenTcellquarter, [1 1 1]);
plot(timerange, meanrandom(:,2),'--','Color','r','LineWidth',width);
if plotmode == 0
   scatter(AverageTcelldata(:,1), AverageTcelldata(:,2),scattersize,'k','filled')
end
if plotmode == 1
   scatter(Tcelldata(:,1), Tcelldata(:,2),scattersize,'k','filled')
end
if plotmode == 2
   scatter(Tcelldata(:,1), Tcelldata(:,2),scattersize,'d','filled',"MarkerFaceColor","#77AC30")
   scatter(AverageTcelldata(:,1), AverageTcelldata(:,2),scattersize,'k','filled')
end
if plotmode == 3
    scatter(Tcelldata(:,1), Tcelldata(:,2),scattersize,'d','filled',"MarkerFaceColor",'b')
    errorbar(AverageTcelldata(:,1), AverageTcelldata(:,2), negTcell, posTcell,'.','color','b')
end

hold off
xlabel('t (days)')
ylabel('T (cells)')
set(gca,'FontSize',fs,'LineWidth',tickwidth)


%---MDSC--Plot-----------------------------------------------------------

ax1 = nexttile(3);
s100 = fill(x2, inBetweenMDSC, [0.25 0.25 0.25]);
hold on;
s75 = fill(x2, inBetweenMDSC3quarter, [0.5 0.5 0.5]);
s50 = fill(x2, inBetweenMDSChalf, [0.75 0.75 0.75]);
s25 = fill(x2, inBetweenMDSCquarter, [1 1 1]);
meanplot = plot(timerange, meanrandom(:,3),'--','Color','r','LineWidth',width);
if plotmode == 0
    scatter(AverageMDSCdata(:,1), AverageMDSCdata(:,2),scattersize,'k','filled')
end
if plotmode == 1
    scatter(MDSCdata(:,1), MDSCdata(:,2),scattersize,'k','filled')
end
if plotmode == 2
    p2 = scatter(MDSCdata(:,1), MDSCdata(:,2),scattersize,'d','filled',"MarkerFaceColor","#77AC30");
    p1 = scatter(AverageMDSCdata(:,1), AverageMDSCdata(:,2),scattersize,'k','filled');
end
if plotmode == 3
    databar = scatter(MDSCdata(:,1), MDSCdata(:,2),scattersize,'d','filled',"MarkerFaceColor",'b');
    errorbar(AverageMDSCdata(:,1), AverageMDSCdata(:,2), negMDSC, posMDSC,'.','Color','b');
    leg = legend(ax1,[databar, meanplot, s25,s50,s75,s100],["data", "mean", "\sigma/4","\sigma/2","3\sigma/4", "\sigma"],'FontSize',fs);
    leg.Layout.Tile = 'east';
end
hold off
xlabel('t (days)')
ylabel('M (cells)')
set(gca,'FontSize',fs,'LineWidth',tickwidth)


set(gcf,'Position',[0 300 900 400],'PaperPositionMode','auto');
