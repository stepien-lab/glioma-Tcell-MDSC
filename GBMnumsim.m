% GBM Numerical Simulations

clear

load('ABCminfornumsim.mat') %ABCfigures.m output (obtained after using minimum data for ABCrejection.m)
minminP = minP;
errormin =  min(error(1:n));

load('ABCmaxfornumsim.mat') %ABCfigures.m output (obtained after using maximum data for ABCrejection.m)
maxminP = minP;
errormax = min(error(1:n));

load('ABCfornumsim.mat') %ABCfigures.m output (obtained after using average data for ABCrejection.m)
F = min(error(1:n));


width = 1.5;
tickwidth =0.8;
scattersize = 25; 

fs = 17;

% Plot modes
    plotmode = 0;   % 0 - average of data points at a certain day
                    % 1 - all data points
                    % 2 - uses both 0 and 1
                    % 3 - errorbar where center is mean and outer bounds
                    % are min and max data points
                    
    datamode = 1;   % 0 - use accepted values for std and mean
                    % 1 - use analysis values for std and mean
               
                    
                    
Legendlabel = ["\lambda_C";
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
      "d_M";
      "";
      "";
      "";
      "";
      "";
      "";
      "";
      "";
      "";
      "";
      "";
      "";
      "";
      "Mean"];

%% RESTRUCTURE DATA FOR PLOT MODE 3

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
                    
%% DETERMINE STANDARD DEVIATIONS AND MEAN 

if datamode==0
    paramdata = histoP;
end

if datamode==1
    paramdata = histoanalysisP;
end

clear mean
meanparam = mean(paramdata)'; %13x1 column vector containing the mean of each parameter

stdparam = std(paramdata)'; %13x1 column vector containing the standard deviation of each parameter

meanminusstdparam = zeros(13,13);
% Mean minus a standard deviation
for i =1:13
    meanminusstdparam(:,i) = meanparam;
    meanminusstdparam(i,i) = meanparam(i) - stdparam(i);
end

meanplusstdparam = zeros(13,13);
% Mean plus a standard deviation
for i =1:13
    meanplusstdparam(:,i) = meanparam;
    meanplusstdparam(i,i) = meanparam(i) + stdparam(i);
end



%% NUMERICAL SIMULATION FOR MEAN AND STANDARD DEVIATIONS 

% Initial Values
C0=35000;
T0=0;
M0=0; 

% Time Range
t0 = 0;
tf = 40;

% Plot

[TlambdaCn,XlambdaCn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,1)), [t0 tf],[C0;T0;M0]);
[TCmaxn,XCmaxn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,2)), [t0 tf],[C0;T0;M0]);
[Tetan,Xetan] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,3)), [t0 tf],[C0;T0;M0]);
[TaTn,XaTn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,4)), [t0 tf],[C0;T0;M0]);
[TsTn,XsTn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,5)), [t0 tf],[C0;T0;M0]);
[Trhon,Xrhon] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,6)), [t0 tf],[C0;T0;M0]);
[TeCn,XeCn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,7)), [t0 tf],[C0;T0;M0]);
[Trn,Xrn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,8)), [t0 tf],[C0;T0;M0]);
[TdTn,XdTn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,9)), [t0 tf],[C0;T0;M0]);
[TsMn,XsMn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,10)), [t0 tf],[C0;T0;M0]);
[Talphan,Xalphan] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,11)), [t0 tf],[C0;T0;M0]);
[Tqn,Xqn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,12)), [t0 tf],[C0;T0;M0]);
[TdMn,XdMn] = ode23s(@(t,x) GBMFunc(t,x,meanminusstdparam(1:13,13)), [t0 tf],[C0;T0;M0]);

[T,X] = ode23s(@(t,x) GBMFunc(t,x,meanparam(1:13)), [t0 tf],[C0;T0;M0]); % Compute Fit Solution

[TlambdaCp,XlambdaCp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,1)), [t0 tf],[C0;T0;M0]);
[TCmaxp,XCmaxp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,2)), [t0 tf],[C0;T0;M0]);
[Tetap,Xetap] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,3)), [t0 tf],[C0;T0;M0]);
[TaTp,XaTp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,4)), [t0 tf],[C0;T0;M0]);
[TsTp,XsTp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,5)), [t0 tf],[C0;T0;M0]);
[Trhop,Xrhop] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,6)), [t0 tf],[C0;T0;M0]);
[TeCp,XeCp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,7)), [t0 tf],[C0;T0;M0]);
[Trp,Xrp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,8)), [t0 tf],[C0;T0;M0]);
[TdTp,XdTp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,9)), [t0 tf],[C0;T0;M0]);
[TsMp,XsMp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,10)), [t0 tf],[C0;T0;M0]);
[Talphap,Xalphap] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,11)), [t0 tf],[C0;T0;M0]);
[Tqp,Xqp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,12)), [t0 tf],[C0;T0;M0]);
[TdMp,XdMp] = ode23s(@(t,x) GBMFunc(t,x,meanplusstdparam(1:13,13)), [t0 tf],[C0;T0;M0]);


%---Tumor--Plot-----------------------------------------------------------

figure
meannumsim = tiledlayout(1,3,'TileSpacing','compact');
nexttile(1)
plot(TlambdaCn,XlambdaCn(:,1),'--','LineWidth',width)
hold on
plot(TCmaxn,XCmaxn(:,1),'--','LineWidth',width)
plot(Tetan,Xetan(:,1),'--','LineWidth',width)
plot(TaTn,XaTn(:,1),'--','LineWidth',width)
plot(TsTn,XsTn(:,1),'--','LineWidth',width)
plot(Trhon,Xrhon(:,1),'k','LineWidth',width)
plot(TeCn,XeCn(:,1),'Color','#0072BD','LineWidth',width)
plot(Trn,Xrn(:,1),'Color','#D95319','LineWidth',width)
plot(TdTn,XdTn(:,1),'Color','#EDB120','LineWidth',width)
plot(TsMn,XsMn(:,1),'Color','#7E2F8E','LineWidth',width)
plot(Talphan,Xalphan(:,1),'Color','#77AC30','LineWidth',width)
plot(Tqn,Xqn(:,1),'Color','#4DBEEE','LineWidth',width)
plot(TdMn,XdMn(:,1),'Color','#A2142F','LineWidth',width)

plot(TlambdaCp,XlambdaCp(:,1),'--','Color','#0072BD','LineWidth',width)
plot(TCmaxp,XCmaxp(:,1),'--','Color','#D95319','LineWidth',width)
plot(Tetap,Xetap(:,1),'--','Color','#EDB120','LineWidth',width)
plot(TaTp,XaTp(:,1),'--','Color','#7E2F8E','LineWidth',width)
plot(TsTp,XsTp(:,1),'--','Color','#77AC30','LineWidth',width)
plot(Trhop,Xrhop(:,1),'k','LineWidth',width)
plot(TeCp,XeCp(:,1),'Color','#0072BD','LineWidth',width)
plot(Trp,Xrp(:,1),'Color','#D95319','LineWidth',width)
plot(TdTp,XdTp(:,1),'Color','#EDB120','LineWidth',width)
plot(TsMp,XsMp(:,1),'Color','#7E2F8E','LineWidth',width)
plot(Talphap,Xalphap(:,1),'Color','#77AC30','LineWidth',width)
plot(Tqp,Xqp(:,1),'Color','#4DBEEE','LineWidth',width)
plot(TdMp,XdMp(:,1),'Color','#A2142F','LineWidth',width)

plot(T,X(:,1),'r','LineWidth',width)

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
plot(TlambdaCn,XlambdaCn(:,2),'--','LineWidth',width)
hold on
plot(TCmaxn,XCmaxn(:,2),'--','LineWidth',width)
plot(Tetan,Xetan(:,2),'--','LineWidth',width)
plot(TaTn,XaTn(:,2),'--','LineWidth',width)
plot(TsTn,XsTn(:,2),'--','LineWidth',width)
plot(Trhon,Xrhon(:,2),'k','LineWidth',width)
plot(TeCn,XeCn(:,2),'Color','#0072BD','LineWidth',width)
plot(Trn,Xrn(:,2),'Color','#D95319','LineWidth',width)
plot(TdTn,XdTn(:,2),'Color','#EDB120','LineWidth',width)
plot(TsMn,XsMn(:,2),'Color','#7E2F8E','LineWidth',width)
plot(Talphan,Xalphan(:,2),'Color','#77AC30','LineWidth',width)
plot(Tqn,Xqn(:,2),'Color','#4DBEEE','LineWidth',width)
plot(TdMn,XdMn(:,2),'Color','#A2142F','LineWidth',width)

plot(TlambdaCp,XlambdaCp(:,2),'--','Color','#0072BD','LineWidth',width)
plot(TCmaxp,XCmaxp(:,2),'--','Color','#D95319','LineWidth',width)
plot(Tetap,Xetap(:,2),'--','Color','#EDB120','LineWidth',width)
plot(TaTp,XaTp(:,2),'--','Color','#7E2F8E','LineWidth',width)
plot(TsTp,XsTp(:,2),'--','Color','#77AC30','LineWidth',width)
plot(Trhop,Xrhop(:,2),'k','LineWidth',width)
plot(TeCp,XeCp(:,2),'Color','#0072BD','LineWidth',width)
plot(Trp,Xrp(:,2),'Color','#D95319','LineWidth',width)
plot(TdTp,XdTp(:,2),'Color','#EDB120','LineWidth',width)
plot(TsMp,XsMp(:,2),'Color','#7E2F8E','LineWidth',width)
plot(Talphap,Xalphap(:,2),'Color','#77AC30','LineWidth',width)
plot(Tqp,Xqp(:,2),'Color','#4DBEEE','LineWidth',width)
plot(TdMp,XdMp(:,2),'Color','#A2142F','LineWidth',width)

plot(T,X(:,2),'r','LineWidth',width)

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

ax = nexttile(3);
lambdaC = plot(TlambdaCn,XlambdaCn(:,3),'--','Color','#0072BD','LineWidth',width);
hold on
Cmax = plot(TCmaxn,XCmaxn(:,3),'--','Color','#D95319','LineWidth',width);
eta =plot(Tetan,Xetan(:,3),'--','Color','#EDB120','LineWidth',width);
aT = plot(TaTn,XaTn(:,3),'--','Color','#7E2F8E','LineWidth',width);
sT = plot(TsTn,XsTn(:,3),'--','Color','#77AC30','LineWidth',width);
rho = plot(Trhon,Xrhon(:,3),'k','LineWidth',width);
eC = plot(TeCn,XeCn(:,3),'Color','#0072BD','LineWidth',width);
r = plot(Trn,Xrn(:,3),'Color','#D95319','LineWidth',width);
dT = plot(TdTn,XdTn(:,3),'Color','#EDB120','LineWidth',width);
sM = plot(TsMn,XsMn(:,3),'Color','#7E2F8E','LineWidth',width);
alpha = plot(Talphan,Xalphan(:,3),'Color','#77AC30','LineWidth',width);
q = plot(Tqn,Xqn(:,3),'Color','#4DBEEE','LineWidth',width);
dM = plot(TdMn,XdMn(:,3),'Color','#A2142F','LineWidth',width);

plot(TlambdaCp,XlambdaCp(:,3),'--','Color','#0072BD','LineWidth',width)
plot(TCmaxp,XCmaxp(:,3),'--','Color','#D95319','LineWidth',width)
plot(Tetap,Xetap(:,3),'--','Color','#EDB120','LineWidth',width)
plot(TaTp,XaTp(:,3),'--','Color','#7E2F8E','LineWidth',width)
plot(TsTp,XsTp(:,3),'--','Color','#77AC30','LineWidth',width)
plot(Trhop,Xrhop(:,3),'k','LineWidth',width)
plot(TeCp,XeCp(:,3),'Color','#0072BD','LineWidth',width)
plot(Trp,Xrp(:,3),'Color','#D95319','LineWidth',width)
plot(TdTp,XdTp(:,3),'Color','#EDB120','LineWidth',width)
plot(TsMp,XsMp(:,3),'Color','#7E2F8E','LineWidth',width)
plot(Talphap,Xalphap(:,3),'Color','#77AC30','LineWidth',width)
plot(Tqp,Xqp(:,3),'Color','#4DBEEE','LineWidth',width)
plot(TdMp,XdMp(:,3),'Color','#A2142F','LineWidth',width)

mean = plot(T,X(:,3),'r','LineWidth',width);

if plotmode == 0
    scatter(AverageMDSCdata(:,1), AverageMDSCdata(:,2),scattersize,'k','filled')
end
if plotmode == 1
    scatter(MDSCdata(:,1), MDSCdata(:,2),scattersize,'k','filled')
end
if plotmode == 2
    scatter(MDSCdata(:,1), MDSCdata(:,2),scattersize,'d','filled',"MarkerFaceColor","#77AC30")
    scatter(AverageMDSCdata(:,1), AverageMDSCdata(:,2),scattersize,'k','filled')
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
lg = legend(ax,Legendlabel,'FontSize',fs);
lg.Layout.Tile = 'east';



%% Histogram of error

histotitle = ["Glioma Error"; "T cell Error"; "MDSC Error"];
medianerror= zeros(1,3);

for i =1:3
    medianerror(i)= median(individualerror(:,i));
    relevantindices = ones(n,1).*(individualerror(:,i)<1);
    relevanterror = individualerror(:,i).*relevantindices;
    nonzerorelevanterror = nonzeros(relevanterror);
    figure
    histogram(nonzerorelevanterror)
    title(histotitle(i))
end

medianerror



%% NUMERICAL SIMULATION FOR MIN, AVERAGE, and MAX
% Plot of parameter combination with the lowest error (out of sampled)

% Initial Values
C0=35000;
T0=0;
M0=0; 

% Time Range
t0 = 0;
tf = 40;

% Plot
[Tmin,Xmin] = ode23s(@(t,x) GBMFunc(t,x,minminP(1:13)), [t0 tf],[C0;T0;M0]);
[T,X] = ode23s(@(t,x) GBMFunc(t,x,minP(1:13)), [t0 tf],[C0;T0;M0]); % Compute Fit Solution
[Tmax,Xmax] = ode23s(@(t,x) GBMFunc(t,x,maxminP(1:13)), [t0 tf],[C0;T0;M0]);


%---Tumor--Plot-----------------------------------------------------------

figure
maxminnumsim = tiledlayout(1,3,'TileSpacing','compact');
nexttile(1)
plot(T,X(:,1),'LineWidth',width)
hold on
plot(Tmin,Xmin(:,1),'LineWidth',width)
plot(Tmax,Xmax(:,1),'LineWidth',width)
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
plot(T,X(:,2),'LineWidth',width)
hold on
plot(Tmin,Xmin(:,2),'LineWidth',width)
plot(Tmax,Xmax(:,2),'LineWidth',width)
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
m2 = plot(T,X(:,3),'LineWidth',width);
hold on
m1 = plot(Tmin,Xmin(:,3),'LineWidth',width);
m3 = plot(Tmax,Xmax(:,3),'LineWidth',width);
if plotmode == 0
    scatter(AverageMDSCdata(:,1), AverageMDSCdata(:,2),scattersize,'k','filled')
    leg = legend(ax1,[m3,m2,m1],["Maximum","Average","Minimum"],'FontSize',fs);
    leg.Layout.Tile = 'east';
end
if plotmode == 1
    scatter(MDSCdata(:,1), MDSCdata(:,2),scattersize,'k','filled')
    leg = legend(ax1,[m3,m2,m1],["Maximum","Average","Minimum"],'FontSize',fs);
    leg.Layout.Tile = 'east';
end
if plotmode == 2
    p2 = scatter(MDSCdata(:,1), MDSCdata(:,2),scattersize,'d','filled',"MarkerFaceColor","#77AC30");
    p1 = scatter(AverageMDSCdata(:,1), AverageMDSCdata(:,2),scattersize,'k','filled');
    leg = legend(ax1,[m3,m2,m1,p1,p2],["Maximum","Average","Minimum", "Average Data","All Data"],'FontSize',fs);
    leg.Layout.Tile = 'east';
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
