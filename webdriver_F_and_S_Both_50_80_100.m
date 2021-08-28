%%HSd50
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_HSd50(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H50Sd50 = BC;    
save('H50Sd50','H50Sd50');
clear all
%%
%%HSd80
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_HSd80(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H50Sd80 = BC;    
save('H50Sd80','H50Sd80');
clear all
%%
%%HSd100
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_HSd100(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H50Sd100 = BC;    
save('H50Sd100','H50Sd100');
clear all
%%
%%HSe50
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_HSe50(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H50Se50 = BC;    
save('H50Se50','H50Se50');
clear all
%%
%%HSe80
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_HSe80(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H50Se80 = BC;    
save('H50Se80','H50Se80');
clear all
%%
%%HSe100
load Data\ChileanIntertidal_ECIMBiomass.mat %biomasas estan en unidades de g/m2
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
    [~,BC] = ode45(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast_HSe100(t,BC,Data), ...
        tspan,[Binit], options);
% ----------------------------------------------------------------------- %
% extinction threshold
BC(BC< 1e-06) = 0;
% ----------------------------------------------------------------------- %
    
    if ~isreal(BC) || any(any(isnan(BC)))
        error('Differential eqution solver returned non real biomasses. Check model parameters.');
    end
	
	cpuTime(cycle) = toc;
    
H50Se100 = BC;    
save('H50Se100','H50Se100');
clear all
