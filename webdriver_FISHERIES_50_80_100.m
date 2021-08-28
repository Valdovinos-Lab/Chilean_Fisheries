%%No harvesting
load Data\ChileanIntertidal_ECIMBiomass.mat 
t_init = 0;
tspan = t_init:3650; 
Data.K.type = 'Constant';
Data.K.mean = 176299; 
Data.K.standard_deviation = [];
Data.K.autocorrelation = [];
Data = updateGuildInfo(Data);
Data = addDerivedQuantities(Data);
Data.nYearsFwd = length(tspan)- 1;
Data = compileOdeData(Data);
GuildInfo = Data.GuildInfo;
Binds = 1:GuildInfo.nGuilds;
Binit = vertcat(Data.Guilds.binit);
%Di = Binit(38,:);
Pl = Binit(95,:);
Pl = Pl - ((Pl*12)/100);
Binit(95,:) = Pl;

cpuTime = zeros(1,Data.nYearsFwd);
    tic
    cycle = 1;
    avgCpuTime = mean(cpuTime(1:cycle-1));
    timeLeftEst = round((Data.nYearsFwd-cycle+1)*avgCpuTime); 
    
    fprintf('Cycle %i/%i starting. BF. Estimated time left %i seconds\n', ...
        cycle, Data.nYearsFwd, timeLeftEst);
      
    Data.year = cycle;
        
    options = odeset('AbsTol',1.000000000000000e-10,'RelTol',1.000000000000000e-08);
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
Before = BC;    
save('Before','Before');
clear all
%%
%%50
load Data\ChileanIntertidal_ECIMBiomass.mat 
t_init = 0;
tspan = t_init:3650; 
Data.K.type = 'Constant';
Data.K.mean = 176299; 
Data.K.standard_deviation = [];
Data.K.autocorrelation = [];
Data = updateGuildInfo(Data);
Data = addDerivedQuantities(Data);
Data.nYearsFwd = length(tspan)- 1;
Data = compileOdeData(Data);
GuildInfo = Data.GuildInfo;
Binds = 1:GuildInfo.nGuilds;
load('Before.mat')
Binit = Before(end,:)';

cpuTime = zeros(1,Data.nYearsFwd);
    tic
    cycle = 1;
    avgCpuTime = mean(cpuTime(1:cycle-1));
    timeLeftEst = round((Data.nYearsFwd-cycle+1)*avgCpuTime); 
    
    fprintf('Cycle %i/%i starting. BF. Estimated time left %i seconds\n', ...
        cycle, Data.nYearsFwd, timeLeftEst);
      
    Data.year = cycle;
        
    options = odeset('AbsTol',1.000000000000000e-10,'RelTol',1.000000000000000e-08);
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_H50(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H50 = BC;    
save('H50','H50');
clear all

%%
%80
load Data\ChileanIntertidal_ECIMBiomass.mat 
t_init = 0;
tspan = t_init:3650; 
Data.K.type = 'Constant';
Data.K.mean = 176299; 
Data.K.standard_deviation = [];
Data.K.autocorrelation = [];
Data = updateGuildInfo(Data);
Data = addDerivedQuantities(Data);
Data.nYearsFwd = length(tspan)- 1;
Data = compileOdeData(Data);
GuildInfo = Data.GuildInfo;
Binds = 1:GuildInfo.nGuilds;
load('Before.mat')
Binit = Before(end,:)';

cpuTime = zeros(1,Data.nYearsFwd);
    tic
    cycle = 1;
    avgCpuTime = mean(cpuTime(1:cycle-1));
    timeLeftEst = round((Data.nYearsFwd-cycle+1)*avgCpuTime); 
    
    fprintf('Cycle %i/%i starting. BF. Estimated time left %i seconds\n', ...
        cycle, Data.nYearsFwd, timeLeftEst);
      
    Data.year = cycle;
        
    options = odeset('AbsTol',1.000000000000000e-10,'RelTol',1.000000000000000e-08);
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_H80(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H80 = BC;    
save('H80','H80');
clear all
%%
%100
load Data\ChileanIntertidal_ECIMBiomass.mat 
t_init = 0;
tspan = t_init:3650; 
Data.K.type = 'Constant';
Data.K.mean = 176299; 
Data.K.standard_deviation = [];
Data.K.autocorrelation = [];
Data = updateGuildInfo(Data);
Data = addDerivedQuantities(Data);
Data.nYearsFwd = length(tspan)- 1;
Data = compileOdeData(Data);
GuildInfo = Data.GuildInfo;
Binds = 1:GuildInfo.nGuilds;
load('Before.mat')
Binit = Before(end,:)';

cpuTime = zeros(1,Data.nYearsFwd);
    tic
    cycle = 1;
    avgCpuTime = mean(cpuTime(1:cycle-1));
    timeLeftEst = round((Data.nYearsFwd-cycle+1)*avgCpuTime); 
    
    fprintf('Cycle %i/%i starting. BF. Estimated time left %i seconds\n', ...
        cycle, Data.nYearsFwd, timeLeftEst);
      
    Data.year = cycle;
        
    options = odeset('AbsTol',1.000000000000000e-10,'RelTol',1.000000000000000e-08);
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_H100(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H100 = BC;  
save('H100','H100');
clear all
