% ABC rejection method
clear 

tic
n = 1000000; % number of sample parameters from range
percenttoaccept = 0.2;


% Modes
    mode = 0;   % 0 - average of data points at a certain day
                % 1 - finds the closest data point
                % 2 - min of data points at a certain day
                % 3 - max of data points at a certain day
% Mouse
    mouse = 0; % mouse 0 is the average of mice 1-4
               % mouse 5 is the min of the mice initial conditions
               % mouse 6 is the max "
               
% Save file?
    savefile = 1; % 0 is no, 1 is yes
    filename = 'ABC.mat';

    
%% Parameters
% the following vector excludes (0) or includes (1) parameters from
% the ABC rejection program

chooseparam = [1;     % 1  - lambdaC
               1;     % 2  - Cmax
               1;     % 3  - η 
               1;     % 4  - a_T 
               1;     % 5  - s_T 
               1;     % 6 - ρ        
               1;     % 7 - ε_C 
               1;     % 8  - r 
               1;     % 9  - d_T
               1;     % 10 - s_M 
               1;     % 11 - α
               1;     % 12 - q
               1];    % 13 - d_M

param = [0.431;      % 1  - lambdaC
         3.04e6;     % 2  - Cmax
         2.57e-8;    % 3  - η 
         3.26e6;     % 4  - a_T 
         8.56e6;     % 5  - s_T 
         0.107;      % 6 - ρ        
         16.2;       % 7 - ε_C 
         6.92e-6;    % 8  - r 
         0.0221;     % 9  - d_T
         0.0372;     % 10 - s_M 
         3.6e8;      % 11 - α
         3.51e10;    % 12 - q
         0.251];     % 13 - d_M
     

%% Data in terms of number of cells

Timedata = [7;13; 20; 24; 27; 34];

TcellTimedata = [7; 24; 27; 34];

MDSCdata = table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',3,'Range','J3:K17'));
% two columns of data containing the average number of MDSCs (column 2)
% on days 7(x4), 13(x2), 20(x4), 24(x2), 27(x2), 34 (in that order) (column 1)

Tcelldata = table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',2,'Range','J3:K8'));
% two columns of data containing the average number of T cells (column 2)
% on days 7, 24(x3), 27, 34 (in that order)

tumordata = table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',1,'Range','J3:K17'));
% two columns of data containing the average number of tumor cells (col 2)
% on days 7(x4), 13(x2), 20(x4), 24(x2), 27(x2), 34 (in that order) (column 1)

AverageMDSCdata = table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',3,'Range','E3:F8'));
% two columns of data containing the average number of MDSCs (column 2)
% on days 7, 13, 20, 24, 27, 34 (in that order) (column 1)

AverageTcelldata = table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',2,'Range','E3:F6'));
% two columns of data containing the average number of T cells (column 2)
% on days 7, 24, 27, 34 (in that order)

Averagetumordata = table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',1,'Range','E3:F8'));
% two columns of data containing the average number of tumor cells (col 2)
% on days 7,13, 20, 24, 27, 34 (in that order) (col 1)

data = zeros(15,3);
data(:,1) = tumordata(:,2);
data(:,2) = MDSCdata(:,2);
data(1:6,3) = Tcelldata(:,2);

AverageData = zeros(6,3);
AverageData(:,1) = Averagetumordata(:,2);
AverageData(:,2) = AverageMDSCdata(:,2);
AverageData(1:4,3) = AverageTcelldata(:,2);
AverageData(5,3) = 0;
AverageData(6,3) = 0;

%%  Data based on mode

if mode == 0 % uses average data
    
% Average Data
MDSCdatainput = AverageData(:,2);

Tcelldatainput = AverageData(1:4,3);

tumordatainput = AverageData(:,1);
    
end

if mode == 1 % for finding nearest data point
    % T cell data
    PRETcelldata = data(1:6,3);
    
    Tcellerrorperday = abs(PRETcelldata-TCELLnumsim);
    
    Tcelldata = zeros(4,1);
    Tcelldata(1) = PRETcelldata(1); %day 7
    Tcelldata(2) = max(PRETcelldata(2:4).*(Tcellerrorperday(2:4) == min(Tcellerrorperday(2:4)))); %day 24
    Tcelldata(3) = PRETcelldata(5); %day 27
    Tcelldata(4) = PRETcelldata(6); %day 34
    
   
    % TUMOR
    Tumorerrorperday = abs(data(:,1)-TUMORnumsim);
   
    tumordatainput =zeros(6,1);
    tumordatainput(1) = max(data(1:4,1).*(Tumorerrorperday(1:4) == min(Tumorerrorperday(1:4))));   %day 7
    tumordatainput(2) = max(data(5:6,1).*(Tumorerrorperday(5:6) == min(Tumorerrorperday(5:6))));  %day 13
    tumordatainput(3) = max(data(7:10,1).*(Tumorerrorperday(7:10) == min(Tumorerrorperday(7:10)))); %day 20
    tumordatainput(4) = max(data(11:12,1).*(Tumorerrorperday(11:12) == min(Tumorerrorperday(11:12))));%day 24
    tumordatainput(5) = max(data(13:14,1).*(Tumorerrorperday(13:14) == min(Tumorerrorperday(13:14))));%day 27
    tumordatainput(6) = data(15,1);   %day 34
    
    
    % MDSC
    MDSCerrorperday = abs(data(:,2)-MDSCnumsim);

    MDSCdatainput =zeros(6,1);
    MDSCdatainput(1) = max(data(1:4,2).*(MDSCerrorperday(1:4) == min(MDSCerrorperday(1:4))));   %day 7
    MDSCdatainput(2) = max(data(5:6,2).*(MDSCerrorperday(5:6) == min(MDSCerrorperday(5:6))));  %day 13
    MDSCdatainput(3) = max(data(7:10,2).*(MDSCerrorperday(7:10) == min(MDSCerrorperday(7:10)))); %day 20
    MDSCdatainput(4) = max(data(11:12,2).*(MDSCerrorperday(11:12) == min(MDSCerrorperday(11:12))));%day 24
    MDSCdatainput(5) = max(data(13:14,2).*(MDSCerrorperday(13:14) == min(MDSCerrorperday(13:14))));%day 27
    MDSCdatainput(6) = data(15,2);   %day 34
end

if mode == 2 % Uses min data points
    
% Minimum Data
    
% MDSC
MDSCdatainput = zeros(6,1);
MDSCdatainput(1) = min(data(1:4,2));   %day 7
MDSCdatainput(2) = min(data(5:6,2));   %day 13
MDSCdatainput(3) = min(data(7:10,2));  %day 20
MDSCdatainput(4) = min(data(11:12,2)); %day 24
MDSCdatainput(5) = min(data(13:14,2)); %day 27
MDSCdatainput(6) = data(15,2);         %day 34
    
    
% T cell data
PRETcelldata = data(1:6,3);
    
Tcelldatainput = zeros(4,1);
Tcelldatainput(1) = PRETcelldata(1); %day 7
Tcelldatainput(2) = min(PRETcelldata(2:4)); %day 24
Tcelldatainput(3) = PRETcelldata(5); %day 27
Tcelldatainput(4) = PRETcelldata(6); %day 34
    
   
% TUMOR
tumordatainput = zeros(6,1);
tumordatainput(1) = min(data(1:4,1));   %day 7
tumordatainput(2) = min(data(5:6,1));   %day 13
tumordatainput(3) = min(data(7:10,1));  %day 20
tumordatainput(4) = min(data(11:12,1)); %day 24
tumordatainput(5) = min(data(13:14,1)); %day 27
tumordatainput(6) = data(15,1);         %day 34

    
end
    
    
if mode == 3 % Uses max data points

% Maximum Data
    
% MDSC
MDSCdatainput = zeros(6,1);
MDSCdatainput(1) = max(data(1:4,2));   %day 7
MDSCdatainput(2) = max(data(5:6,2));   %day 13
MDSCdatainput(3) = max(data(7:10,2));  %day 20
MDSCdatainput(4) = max(data(11:12,2)); %day 24
MDSCdatainput(5) = max(data(13:14,2)); %day 27
MDSCdatainput(6) = data(15,2);         %day 34
    
    
% T cell data
PRETcelldata = data(1:6,3);
    
Tcelldatainput = zeros(4,1);
Tcelldatainput(1) = PRETcelldata(1); %day 7
Tcelldatainput(2) = max(PRETcelldata(2:4)); %day 24
Tcelldatainput(3) = PRETcelldata(5); %day 27
Tcelldatainput(4) = PRETcelldata(6); %day 34
    
   
% TUMOR
tumordatainput = zeros(6,1);
tumordatainput(1) = max(data(1:4,1));   %day 7
tumordatainput(2) = max(data(5:6,1));   %day 13
tumordatainput(3) = max(data(7:10,1));  %day 20
tumordatainput(4) = max(data(11:12,1)); %day 24
tumordatainput(5) = max(data(13:14,1)); %day 27
tumordatainput(6) = data(15,1);         %day 34

    
end


%% Initial Conditions

t0=7; %actually day 7
tf=34; %actually day 34
    
timerange = [t0 13 20 24 27 tf];
    
%Initial (average) Conditions at day 7
T7= 66131.7773; % we only have one T cell data point for day 7
    
if mouse == 0 %average
    C7= 135472.1199;
    M7= 4912.9204;
end
    
if mouse == 1
    C7 = data(1,1);
    M7 = data(1,2);
end
    
if mouse == 2
    C7 = data(2,1);
    M7 = data(2,2);
end
    
if mouse == 3
    C7 = data(3,1);
    M7 = data(3,2);
end
    
if mouse == 4
    C7 = data(4,1);
    M7 = data(4,2);
end
    
if mouse == 5 %min
    C7 = min(data(1:4,1));
    M7 = min(data(1:4,2));
end
    
if mouse == 6 %max
    C7 = max(data(1:4,1));
    M7 = max(data(1:4,2));
end
    
initialcondition = [C7;T7;M7];

%% Step 1: Specify prior distribution
% specify range for parameter (a,b)
% we have no previous inclination, so performed uniform sampling of
% parameter range

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

% MODEL PARAMETER RANGES: (Min/Max) 

mpr = [0            0.5;        % 1  - lambdaC 
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
   
rng(0, 'twister'); % to uniformly generate random numbers 

%% Step 2: Sample randomly from the distribution of each of the parameters
% for each sampling, you get a single set of parameters

Psample = zeros(n,13);

for i = 1:13
    if chooseparam(i) == 1
        Psample(:,i) = (mpr(i,2) - mpr(i,1)).*rand(n,1) + mpr(i,1);
    end
    if chooseparam(i) == 0
        Psample(:,i) = param(i);
    end
end

toc
%% Step 3: Calculate error using the ABC error function

individualerror = zeros(n,3);

for i = 1:n
    [i]
    individualerror(i,:) = e(Psample(i,1:13), MDSCdatainput, Tcelldatainput, tumordatainput, mode, initialcondition, timerange);
end

error = sum(individualerror,2);

findthreshold= sort(individualerror);
errorthreshold=findthreshold(ceil(n*percenttoaccept),:);

toc 

%% SAVE FILE

if savefile == 1
    save(filename)
end

%% COMPLETE

f = msgbox("The operation has been successfully completed","DONE");

%% OBJECTIVE FUNCTION:
function E = e(Psample, MDSCdata, Tcelldata, tumordata, mode, initialcondition, timerange)
    
    % EVALUATE FUNCTION WITH CURRENT GUESS:
    [T1,X1] = ode23s(@(t,x) GBMFunc(t,x,Psample), timerange, initialcondition);
   
    
    if mode ~= 1 % for modes 0 (average), 2 (min), 3 (max)
        % Update num sim values
    TUMORnumsim = [X1(1,1);  % day 7
        X1(2,1);  % day 13
        X1(3,1);  % day 20 
        X1(4,1);  % day 24
        X1(5,1);  % day 27  
        X1(6,1)]; % day 34

    MDSCnumsim = [X1(1,3); % day 7
        X1(2,3);  % day 13  
        X1(3,3);  % day 20 
        X1(4,3);  % day 24  
        X1(5,3);  % day 27
        X1(6,3)]; % day 34

    TCELLnumsim = [X1(1,2);  % day 7
        X1(4,2); % day 24 
        X1(5,2); % day 27
        X1(6,2)];% day 34
    
    end
    
     
    if mode == 1 % find nearest data point
        % Update num sim values
    TUMORnumsim = [X1(1,1);  % day 7 (x4)
        X1(1,1);  
        X1(1,1);  
        X1(1,1);  
        X1(2,1);  % day 13 (x2)
        X1(2,1);  
        X1(3,1);  % day 20 (x4)
        X1(3,1);  
        X1(3,1);  
        X1(3,1);  
        X1(4,1);  % day 24 (x2)
        X1(4,1);  
        X1(5,1);  % day 27 (x2)
        X1(5,1);  
        X1(6,1)]; % day 34

    MDSCnumsim = [X1(1,3); % day 7 (x4)
        X1(1,3); 
        X1(1,3); 
        X1(1,3); 
        X1(2,3);  % day 13 (x2)
        X1(2,3);  
        X1(3,3);  % day 20 (x4)
        X1(3,3);  
        X1(3,3);  
        X1(3,3);  
        X1(4,3);  % day 24 (x2)
        X1(4,3);  
        X1(5,3);  % day 27 (x2)
        X1(5,3);  
        X1(6,3)]; % day 34

    TCELLnumsim = [X1(1,2);  % day 7
        X1(4,2); % day 24 (x3)
        X1(4,2); 
        X1(4,2); 
        X1(5,2); % day 27
        X1(6,2)];% day 34
    
    end
    
    % COMPUTE ERROR:
    tumorerror = relativeError(tumordata,TUMORnumsim,6);
    Tcellerror = relativeError(Tcelldata,TCELLnumsim,4);
    MDSCerror = relativeError(MDSCdata,MDSCnumsim,6);
    
    E = [tumorerror Tcellerror MDSCerror];    
   
end

%% Mean Square Error/Quadratic Loss/L2 Loss
function g = MSE(data,numsim,n)

x=0;

for i = 1:n
    xx=(data(i)-numsim(i))^2;
    x = x + xx;
end

g = x/n;

end


%% Mean Absolute Error/L1 Loss

function g = MAE(data,numsim,n)

y=0;

for i = 1:n
    yy=abs(data(i)-numsim(i));
    y = y + yy;
end

g = y/n;
end


%% Mean Bias Error

function g = MBE(data,numsim,n)

z=0;

for i = 1:n
    zz=(data(i)-numsim(i));
    z = z + zz;
end

g = z/n;
end

%% Relative Error

function g = relativeError(data,numsim,n)

z=0;

for i = 1:n
    zz=abs((data(i)-numsim(i))/data(i));
    z = z + zz;
end


g = z/n; 
end