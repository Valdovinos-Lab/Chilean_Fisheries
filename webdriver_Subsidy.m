%%
%Before
%
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
%Deplet subsidy 50%
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_D50(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
D50 = BC;    
save('D50','D50');

clear all

%%
%%Deplet subsidy 80%
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_D80(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
D80 = BC;    
save('D80','D80');
clear all
%%
%%Deplet subsidy 100%
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_D100 (t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
D100 = BC;    
save('D100','D100');
clear all
%%
%Enrich subsidy 50%
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_E50(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
E50 = BC;    
save('E50','E50');
clear all
%%
%Enrich subsidy 80%
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_E80(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
E80 = BC;    
save('E80','E80');
clear all
%%
%Enrich subsidy 100%
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_E100(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
E100 = BC;    
save('E100','E100');
clear all

