%% Before Fsheries KL

topmin = round(150/4);
topmax = 150*4;
Anmin = round(500/4);
Anmax = 500*4;
Ppmin = round(15000/4);
Ppmax = 15000*4;

Opt = cell(4, 2);%(10 es el num total de iteraciones que seran mostradas en los raws y el 2 dos columnas)
Bio = cell(4, 2);
BTL = cell(4, 2);
Err = cell(4, 2);
EvB0 = cell(4, 2);

for i = 1:4
B0top = randi([topmin,topmax],1,1);
B0an = randi([Anmin,Anmax],1,1);
B0pp = randi([Ppmin,Ppmax],1,1);
OptBo = horzcat(B0top,B0an,B0pp);
EvaluationB0 = OptBo(1,1) < OptBo(1,2) && OptBo(1,2) < OptBo(1,3);

Opt(i, :) = {i OptBo}; %para hacer el cell que guarda los datos
RandomB0s(i,:) = cell2mat(Opt(i,2)); %para pasar de cell a matrix, las columnas son las especies y las filas las iteraciones

EvB0(i, :) = {i EvaluationB0}; %para hacer el cell que guarda los datos
BoEvaluation(i,:) = cell2mat(EvB0(i,2));

load Data\ChileanIntertidal_ECIMBiomass.mat %biomasas estan en unidades de g/m2
t_init = 0;
tspan = t_init:700; %significa que va a hacer 3
Data.K.type = 'Constant';
Data.K.mean = 176299; %Using Pablo criteria were I calculated by each group of algae
Data.K.standard_deviation = [];
Data.K.autocorrelation = [];

B0_opt = Data.B0_6;
B0_opt (B0_opt == 150) = B0top;
B0_opt (B0_opt == 500) = B0an;
B0_opt (B0_opt == 15000) = B0pp;
Data.B0_opt = B0_opt;

Data = updateGuildInfo(Data);
Data = addDerivedQuantities(Data);
Data.nYearsFwd = length(tspan)- 1;
Data = compileOdeData(Data);
GuildInfo = Data.GuildInfo;
Binds = 1:GuildInfo.nGuilds;
Binit = vertcat(Data.Guilds.binit);
Di = Binit(38,:);
Pl = Binit(95,:);
Di = Di - ((Di*12)/100);
Pl = Pl - ((Pl*12)/100);
Binit(38,:) = Di;
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
 
    gTL = {Data.Guilds.TL};
    gTL = gTL';
    iTop = find(strcmpi(gTL,'''Top'''));
    iPredator = find(strcmpi(gTL,'''Predator'''));
    iOmnivore = find(strcmpi(gTL,'''Omnivore'''));
    iHerbivore = find(strcmpi(gTL,'''Herbivore''')); 
    iFilter = find(strcmpi(gTL,'''Filter'''));
    iProducers = find(strcmpi(gTL,'''Producers'''));
    
    SimulatedBiomass = BC(end,:)';   
    EmpiricalBiomass = vertcat(Data.Guilds.binit);
    
    %COmmunity
    ObsCom = sum(SimulatedBiomass);
    EmpCom = sum(EmpiricalBiomass); 
    %TOP
    OBTop = SimulatedBiomass(iTop);
    OBTop = sum(OBTop);
    SBTop = EmpiricalBiomass(iTop);
    SBTop = sum(SBTop);
    %CARNIVORES
    OBPred = SimulatedBiomass(iPredator);
    OBPred = sum(OBPred);
    SBPred = EmpiricalBiomass(iPredator);
    SBPred = sum(SBPred);
    %Omnivores
    OBOmn = SimulatedBiomass(iOmnivore);
    OBOmn = sum(OBOmn);
    SBOmn = EmpiricalBiomass(iOmnivore);
    SBOmn = sum(SBOmn);
    %Hervibores
    OBHerb = SimulatedBiomass(iHerbivore);
    OBHerb = sum(OBHerb);
    SBHerb = EmpiricalBiomass(iHerbivore);
    SBHerb = sum(SBHerb);
    %Filters
    OBFil = SimulatedBiomass(iFilter);
    OBFil = sum(OBFil);
    SBFil = EmpiricalBiomass(iFilter);
    SBFil = sum(SBFil);
    %Producers 
    OBPp = SimulatedBiomass(iProducers);
    OBPp = sum(OBPp);
    SBPp = EmpiricalBiomass(iProducers);
    SBPp = sum(SBPp);
    
    Simulated = vertcat(ObsCom, OBTop, OBPred, OBOmn, OBHerb, OBFil, OBPp);
    Empirical = vertcat(EmpCom, SBTop, SBPred, SBOmn, SBHerb, SBFil, OBPp);
    Error = abs((Simulated-Empirical)./Empirical);
    
    %Error = abs((SimulatedBiomass-EmpiricalBiomass)./EmpiricalBiomass);
    
    %save biomass of t=700
    Bio(i, :) = {i SimulatedBiomass}; %para hacer el cell que guarda los datos
    Biomass(i,:) = cell2mat(Bio(i,2));
    %save biomass by TL
    BTL(i, :) = {i Simulated}; %para hacer el cell que guarda los datos
    BiomassTL(i,:) = cell2mat(BTL(i,2));
    
    %save erros
    Err(i, :) = {i Error}; %para hacer el cell que guarda los datos
    Errors(i,:) = cell2mat(Err(i,2));
        
	cpuTime(cycle) = toc;
    
end
save('Optimization\Biomass','Biomass')
save('Optimization\BiomassTL','BiomassTL')
save('Optimization\Errors','Errors')
save('Optimization\RandomB0s','RandomB0s')
save('Optimization\BoEvaluation','BoEvaluation')

%Biomass = Biomass';
%BiomassTL = BiomassTL'
%Errors = Errors';
%RandomB0s;
