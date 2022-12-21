% PLOTS FOR ABC METHOD

clear

load('ABC.mat') %load output from ABCrejection.m


fs = 17; % font size
width = 1.15;
scattersize = 100; 

errorthreshold = [0.75; 0.72; 0.78]; %Tumor, T cell, and MDSC threshold


% Plot modes
    plotmode = 1;   % 0 - average of data points at a certain day
                    % 1 - all data points   
% Plot modes                 
    percenttoanalyze = 0.25;                 

% Choose which parameters to graph   
    chooseparam = [1;     % 1  - lambdaC
               1;     % 2  - Cmax
               1;     % 3  - η 
               0;     % 4  - a_T 
               0;     % 5  - s_T 
               1;     % 6 - ρ        
               1;     % 7 - ε_C 
               1;     % 8  - r 
               0;     % 9  - d_T
               0;     % 10 - s_M 
               0;     % 11 - α
               0;     % 12 - q
               0];    % 13 - d_M        

axislabels = {'0',            '0.5';            % 1  - lambdaC
              '10^6',         '5 \times 10^7';  % 2  - Cmax
              '0',            '10^{-6}';        % 3  - η 
              '50',           '5 \times 10^6';  % 4  - a_T 
              '100',          '10^7';           % 5  - s_T
              '0',            '0.5';            % 6  - ρ        
              '1',            '100';            % 7  - ε_C 
              '0',            '10^{-4}';        % 8  - r 
              '0',            '0.75';           % 9  - d_T
              '0',            '0.1';            % 10 - s_M 
              '10^7',         '5 \times 10^8';  % 11 - α
              '10^9',         '10^{11}';        % 12 - q
              '0',            '0.5'};           % 13 - d_M

%% Step 4: Keep parameter sets with error below a determined threshhold

error = sum(individualerror,2);

COnes = ones(n,1).*(individualerror(:,1)<errorthreshold(1));
CTOnes = COnes.*(individualerror(:,2)<errorthreshold(2));
CTMOnes = CTOnes.*(individualerror(:,3)<errorthreshold(3));

acceptedindividualerror = individualerror.*CTMOnes;
acceptederror = error.*CTMOnes;
P = Psample.*CTMOnes;


%% Step 5a: Keep all of these parameter sets (in Step 4) 
% these are posterior samples to represent the distribution of parameter 
% values that are optimal

[sortacceptederror, sortindices]=sort(acceptederror,'descend');

nonzeroacceptederror = nonzeros(sortacceptederror);
l = length(nonzeroacceptederror);
percentofacceptedparameters = l/n;

sortP = P(sortindices,:);
nonzeroP = nonzeros(sortP); 
histoP = reshape(nonzeroP,l,13);


%% Step 5b: Do analyses with smallest (specified)% (in terms of accepted error) of these parameter sets

analysiserror = zeros(l,1);
analysisP = zeros(l,13);
analysisthresholdindex = ceil(l*(1-percenttoanalyze)); 
analysiserrorthreshold = nonzeroacceptederror(analysisthresholdindex,1); 
for i = 1:l
     if nonzeroacceptederror(i) < analysiserrorthreshold
            analysiserror(i) = nonzeroacceptederror(i);
            analysisP(i,1:13) = histoP(i,1:13);    
     end
end
nonzeroanalysiserror = nonzeros(analysiserror);
    
nonzeroanalysisP = nonzeros(analysisP); 
m = length(nonzeroanalysisP)/13;
histoanalysisP = reshape(nonzeroanalysisP,m,13);
      
%% Parameter combination with the lowest error (out of sampled)

% Parameter set with minimum error
minP = zeros(13,1);  % parameter set of total minimum error
minPofC = zeros(13,1); % " of tumor minimum error
minPofT = zeros(13,1); % " T cell "
minPofM = zeros(13,1); % " MDSC "
for i = 1:13
    minP(i,1) = max(Psample(1:n,i).*(error(1:n,1) == min(nonzeroacceptederror)));
    minPofC(i,1) = max(Psample(1:n,i).*(individualerror(1:n,1) == min(individualerror(1:n,1))));
    minPofT(i,1) = max(Psample(1:n,i).*(individualerror(1:n,2) == min(individualerror(1:n,2))));
    minPofM(i,1) = max(Psample(1:n,i).*(individualerror(1:n,3) == min(individualerror(1:n,3))));
end

minerror = min(error(1:n));
minindividualerror = zeros(1,3); % the minimum of each individual error assessed
individualofminerror = zeros(1,3); % the individual errors associated with the TOTAL minimum error
for i = 1:3
    individualofminerror(i) = max(individualerror(:,i).*(error(1:n) == min(nonzeroacceptederror)));
    minindividualerror(i) = min(individualerror(:,i));
end


%% 13x13 SCATTER PLOT OF PARAMETERS AND ERRORS
% histograms for parameters are on the diagonals

chosenparam = find(chooseparam); % vector of the indices with 1 in choose param
noc = length(chosenparam); % number of chosen (parameters)

% SCATTER PLOT USING ALL ACCEPTED PARAMETER VALUES
figure
% i is row, j is column
acceptedscatter = tiledlayout(noc,noc,'TileSpacing','compact');
for i = 1:noc
    for j = 1:noc
        if i >j
            ax4 = nexttile(noc*(i-1) + j); %position in 169
            scatter(histoP(:,chosenparam(j)),histoP(:,chosenparam(i)),[],nonzeroacceptederror)
            colormap(hot(100)) 
            hold on
            scatter(minP(chosenparam(j)), minP(chosenparam(i)),scattersize,'k','d','filled')
            scatter(minPofC(chosenparam(j)), minPofC(chosenparam(i)),scattersize,'b','p','filled')
            scatter(minPofT(chosenparam(j)), minPofT(chosenparam(i)),scattersize,'h','filled', "MarkerFaceColor","#77AC30")
            scatter(minPofM(chosenparam(j)), minPofM(chosenparam(i)),scattersize,'s','filled',"MarkerFaceColor","#7E2F8E")
            hold off
            if i ~= noc
                set(gca,'XTick',[])
            end
            if i == noc 
                xlabel(sprintf('%s',pn(chosenparam(j))))
                xticks(mpr(chosenparam(j),:))
                xticklabels(axislabels(chosenparam(j),:))
                set(gca,'FontSize',fs)
            end
            if j ~= 1
                set(gca,'YTick',[])
            end
            if j == 1
                ylabel(sprintf('%s',pn(chosenparam(i))))
                yticks(mpr(chosenparam(i),:))
                yticklabels(axislabels(chosenparam(i),:))
                set(gca,'FontSize',fs)
            end
        
        end
        if i == j
            [y,x] = ksdensity(histoanalysisP(:,chosenparam(i)));
            nexttile(noc*(i-1)+i)
            plot(x,y,'k','LineWidth',width);
            if i==noc
                set(gca,'YTick',[],'FontSize',fs)
            else
                set(gca,'XTick',[],'YTick',[])
            end
            if i == noc 
                xlabel(sprintf('%s',pn(chosenparam(j))))
                xticks(mpr(chosenparam(j),:))
                xticklabels(axislabels(chosenparam(j),:))
            end
        end
    end
end
set(gcf,'Position',[0 300 1100 800],'PaperPositionMode','auto');
C=colorbar(ax4(end));
C.Layout.Tile='west';
title(C,'error','Interpreter','Latex','FontSize',fs);
lg = legend(ax4,{'','total','glioma','T cell','MDSC'},'FontSize',fs); %, 'Orientation', 'Vertical');
lg.Layout.Tile = 2*noc;
lg.Title.String = 'Minimum Relative Error';     
title(acceptedscatter,'Scatter plot of accepted parameters','FontSize',fs+2)


% SCATTER PLOT USING ONLY ANALYSIS PARAMETER VALUES (for instance, the 25% of
% accepted parameter values with the lowest error)
figure
analysisscatter = tiledlayout(noc,noc,'TileSpacing','compact');
for i = 1:noc
    for j = 1:noc
        if j<i
            ax5 = nexttile(noc*(i-1) + j); %position in 169 (13x13 block)
            scatter(histoanalysisP(:,chosenparam(j)),histoanalysisP(:,chosenparam(i)), [],nonzeroanalysiserror)
            colormap(hot(100)) 
            hold on
            scatter(minP(chosenparam(j)), minP(chosenparam(i)),scattersize,'k','d','filled')
            scatter(minPofC(chosenparam(j)), minPofC(chosenparam(i)),scattersize,'b','p','filled')
            scatter(minPofT(chosenparam(j)), minPofT(chosenparam(i)),scattersize,'h','filled',"MarkerFaceColor","#77AC30")
            scatter(minPofM(chosenparam(j)), minPofM(chosenparam(i)),scattersize,'s','filled',"MarkerFaceColor","#7E2F8E")
            hold off
            if i ~= noc
                set(gca,'XTick',[])
            end
            if i == noc 
                xlabel(sprintf('%s',pn(chosenparam(j))))
                xticks(mpr(chosenparam(j),:))
                xticklabels(axislabels(chosenparam(j),:))
                set(gca,'FontSize',fs)
            end
            if j ~= 1
                set(gca,'YTick',[])
            end
            if j == 1
                ylabel(sprintf('%s',pn(chosenparam(i))))
                yticks(mpr(chosenparam(i),:))
                yticklabels(axislabels(chosenparam(i),:))
                set(gca,'FontSize',fs)
            end
        
        end
        if i == j
            [y,x] = ksdensity(histoanalysisP(:,chosenparam(i)));
            nexttile(noc*(i-1)+i)
            plot(x,y,'k','LineWidth',width);
            if i==noc
                set(gca,'YTick',[],'FontSize',fs)
            else
                set(gca,'XTick',[],'YTick',[])
            end
            if i == noc 
                xlabel(sprintf('%s',pn(chosenparam(j))))
                xticks(mpr(chosenparam(j),:))
                xticklabels(axislabels(chosenparam(j),:))
            end
        end
    end
end
set(gcf,'Position',[0 300 1100 800],'PaperPositionMode','auto');
C=colorbar(ax5(end));
C.Layout.Tile='west';
title(C,'error','Interpreter','Latex','FontSize',fs);
lg = legend(ax5,{'','total','glioma','T cell','MDSC'},'FontSize',fs); 
lg.Layout.Tile = 2*noc;
lg.Title.String = 'Minimum Relative Error';
title(analysisscatter,'Scatter plot of analysis parameters','FontSize',fs+2)
       

%% 13x13 CONTOUR PLOT OF PARAMETERS

% CONTOUR PLOT USING ONLY ACCEPTED PARAMETER VALUES
figure
acceptedcontour = tiledlayout(noc,noc,'TileSpacing','compact');
for i = 1:noc
    for j = 1:noc
        if chosenparam(i) > chosenparam(j)
            ax2 = nexttile(noc*(i-1) + j); %position in noc x noc block
            x = histoP(:,chosenparam(j));
            y = histoP(:,chosenparam(i));
            [pdfx,xi] = ksdensity(x);
            [pdfy,yi] = ksdensity(y);
            [xxi,yyi] = meshgrid(xi,yi);
            [pdfxx,pdfyy] = meshgrid(pdfx,pdfy);
            pdfxy = pdfxx.*pdfyy;
            contourf(xxi,yyi,pdfxy,'LineColor','none')
            colormap(flipud(hot))
            hold on
            scatter(minP(chosenparam(j)), minP(chosenparam(i)),scattersize,'k','d','filled')
            scatter(minPofC(chosenparam(j)), minPofC(chosenparam(i)),scattersize,'b','p','filled')
            scatter(minPofT(chosenparam(j)), minPofT(chosenparam(i)),scattersize,'h','filled',"MarkerFaceColor","#77AC30")
            scatter(minPofM(chosenparam(j)), minPofM(chosenparam(i)),scattersize,'s','filled',"MarkerFaceColor","#7E2F8E")
            hold off
            if i ~= noc
                set(gca,'XTick',[])
            end
            if i == noc 
                xlabel(sprintf('%s',pn(chosenparam(j))))
                xticks(mpr(chosenparam(j),:))
                xticklabels(axislabels(chosenparam(j),:))
                set(gca,'FontSize',fs)
            end
            if j ~= 1
                set(gca,'YTick',[])
            end
            if j == 1
                ylabel(sprintf('%s',pn(chosenparam(i))))
                yticks(mpr(chosenparam(i),:))
                yticklabels(axislabels(chosenparam(i),:))
                set(gca,'FontSize',fs)
            end
        end
        if i == j
            [y,x] = ksdensity(histoP(:,chosenparam(i)));
            nexttile(noc*(i-1)+i)
            plot(x,y,'k','LineWidth',width);
            if i==noc
                set(gca,'YTick',[],'FontSize',fs)
                xticks(mpr(chosenparam(i),:))
                xticklabels(axislabels(chosenparam(i),:))
            else
                set(gca,'XTick',[],'YTick',[])
            end
            set(gcf,'Position',[0 300 1000 800],'PaperPositionMode','auto');
            if i == noc 
                xlabel(sprintf('%s',pn(chosenparam(j))))
                xticks(mpr(chosenparam(j),:))
                xticklabels(axislabels(chosenparam(j),:))
            end
        end
    end
end
lgd = legend(ax2,{'','total','glioma','T cell','MDSC'},'FontSize',fs); 
lgd.Layout.Tile = 2*noc;
lgd.Title.String = 'Minimum Relative Error';
title(acceptedcontour,'Contour plot for accepted parameters','FontSize',fs+2)


% CONTOUR PLOT USING ONLY ANALYSIS PARAMETERS 
figure
analysiscontour = tiledlayout(noc,noc,'TileSpacing','compact');
for i = 1:noc
    for j = 1:noc
        if chosenparam(i) > chosenparam(j)
            ax3 = nexttile(noc*(i-1) + j); %position in noc x noc block
            x = histoanalysisP(:,chosenparam(j));
            y = histoanalysisP(:,chosenparam(i));
            [pdfx,xi] = ksdensity(x);
            [pdfy,yi] = ksdensity(y);
            [xxi,yyi] = meshgrid(xi,yi);
            [pdfxx,pdfyy] = meshgrid(pdfx,pdfy);
            pdfxy = pdfxx.*pdfyy;
            contourf(xxi,yyi,pdfxy,'LineColor','none')
            colormap(flipud(hot))
            hold on
            scatter(minP(chosenparam(j)), minP(chosenparam(i)),scattersize,'k','d','filled')
            scatter(minPofC(chosenparam(j)), minPofC(chosenparam(i)),scattersize,'b','p','filled')
            scatter(minPofT(chosenparam(j)), minPofT(chosenparam(i)),scattersize,'h','filled',"MarkerFaceColor","#77AC30")
            scatter(minPofM(chosenparam(j)), minPofM(chosenparam(i)),scattersize,'s','filled',"MarkerFaceColor","#7E2F8E")
            hold off
            if i ~= noc
                set(gca,'XTick',[])
            end
            if i == noc 
                xlabel(sprintf('%s',pn(chosenparam(j))))
                xticks(mpr(chosenparam(j),:))
                xticklabels(axislabels(chosenparam(j),:))
                set(gca,'FontSize',fs)
            end
            if j ~= 1
                set(gca,'YTick',[])
            end
            if j == 1
                ylabel(sprintf('%s',pn(chosenparam(i))))
                yticks(mpr(chosenparam(i),:))
                yticklabels(axislabels(chosenparam(i),:))
                set(gca,'FontSize',fs)
            end
        end
        if i == j
            [y,x] = ksdensity(histoanalysisP(:,chosenparam(i)));
            nexttile(noc*(i-1)+i)
            plot(x,y,'k','LineWidth',width);
            if i==noc
                set(gca,'YTick',[],'FontSize',fs)
                xticks(mpr(chosenparam(i),:))
                xticklabels(axislabels(chosenparam(i),:))
            else
                set(gca,'XTick',[],'YTick',[])
            end
            set(gcf,'Position',[0 300 1000 800],'PaperPositionMode','auto');
            if i == noc 
                xlabel(sprintf('%s',pn(chosenparam(j))))
                xticks(mpr(chosenparam(j),:))
                xticklabels(axislabels(chosenparam(j),:))
            end
        end
    end
end
leg = legend(ax3,{'','total','glioma','T cell','MDSC'},'FontSize',fs); %, 'Orientation', 'Vertical');
leg.Layout.Tile = 2*noc;
leg.Title.String = 'Minimum Relative Error';   
title(analysiscontour,'Contour plot of analysis parameters','FontSize',fs+2)