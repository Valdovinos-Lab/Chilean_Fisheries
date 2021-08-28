function dXdt = atn_ode_biomass_catch_gainfish_gain_loss_fast_E50(~,X,Data)
% ----------------------------------------------------------------------- %
% ----------------------------- PARAMETERS ------------------------------ %
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% FEEDING LINKS
% ----------------------------------------------------------------------- %

% Feeding links between feeders and food
communityMatrix = Data.communityMatrix;
% Functional response exponent parameters for each feeding link
q = Data.q;
%q(isnan(q))=0;
% Feeding interference coefficients for each feeding link
d = Data.d;
%d(isnan(d))=0;
% Relative food preferences of feeders for each feeding link
omega = Data.omega;
% Assimilation efficiencies for each feeding link
e = Data.e;
%e(isnan(e))=0;
% Maximum consumption rates for each feeding link
y = Data.y;
%y(isnan(y))=0;
% Half saturation constants to the power of q for each feeding link
B0_pow_q = Data.B0_pow_q;

% ----------------------------------------------------------------------- %
% PRODUCERS
% ----------------------------------------------------------------------- %

% Carrying capacity of autotrophs for current year
K = Data.K.values(Data.year);
% Intrinsic growth rate of producers
r = Data.r;
% Competition coefficients of producers
c = Data.c;

% ----------------------------------------------------------------------- %
% CONSUMERS AND FISHES
% ----------------------------------------------------------------------- %
% Metabolic rate of consumers and fishes
x = Data.x;
% Fractions of assimilated C used for biomass gains for feeders
ax = Data.ax;
% Fractions of assimilated C respired by maintenance for feeders
bx = Data.bx;


% ----------------------------------------------------------------------- %
% MISC
% ----------------------------------------------------------------------- %

% Precalculated index vectors etc.
GuildInfo = Data.GuildInfo;
% New variables for speed up by avoiding repetitive computations or indexing
B = X(1:GuildInfo.nGuilds);  

% ----------------------------------------------------------------------- %
% extinction threshold
B(B< 1.00000000000000e-06) = 0;
% ----------------------------------------------------------------------- %

% Producer biomasses
bProducerGuilds = B(GuildInfo.iProducerGuilds);
% Consumer biomasses
bConsumerGuilds = B(GuildInfo.iConsumerGuilds);
% Feeder biomasses
bFeederGuilds = B(GuildInfo.iFeederGuilds);
% Food biomasses
bFoodGuilds = B(GuildInfo.iFoodGuilds);
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% OMEGA VARIABLE
% ----------------------------------------------------------------------- %
A = bFoodGuilds'.*Data.communityMatrix;
A(A > 0) = 1;
ng = size(A,2);
nps = 1./sum(A,2);
nps(isinf(nps)) = 0;
w = repmat(nps,1,ng);
omega = A.*w; 

% ----------------------------------------------------------------------- %
% Fisheries
% ----------------------------------------------------------------------- %
Catchable = vertcat(Data.Guilds.catchablenew);
Catch = Catchable(1:GuildInfo.nGuilds);
ProducerHarvest = Catch(GuildInfo.iProducerGuilds);
ConsumerHarvest = Catch(GuildInfo.iConsumerGuilds);

ProducerFisheries = ((0*bProducerGuilds)/100).*ProducerHarvest;
ConsumerFisheries = ((0*bConsumerGuilds)/100).*ConsumerHarvest;

% ----------------------------------------------------------------------- %
% -------------------------------- MODEL -------------------------------- %
% ----------------------------------------------------------------------- %

% Preallocate biomass derivative vector
dBdt = zeros(GuildInfo.nGuilds,1);

% ----------------------------------------------------------------------- %
% FUNCTIONAL RESPONSE MATRIX
% ----------------------------------------------------------------------- %

Bmx_pow_q = repmat(bFoodGuilds',GuildInfo.nFeederGuilds,1).^q;
Ftop = omega.*repmat(bFoodGuilds',GuildInfo.nFeederGuilds,1).^q;
Fbot1 = B0_pow_q;
Fbot2 = diag(sum(omega.*Bmx_pow_q,2))*communityMatrix; 
Fbot3 = diag(bFeederGuilds)*d.*B0_pow_q;
Fbot = Fbot1 + Fbot2 + Fbot3;
F = Ftop./Fbot;
F(isnan(F)) = 0;

% ----------------------------------------------------------------------- %
% GAINS AND LOSSES
% ----------------------------------------------------------------------- %

ytimesxtimesBtimesF = y.*repmat(x.*bFeederGuilds,1,GuildInfo.nFoodGuilds).*F;

% Gain from resources
gainMTX = diag(ax)*ytimesxtimesBtimesF;
gainVEC = sum(gainMTX,2);


% Loss to consumers
lossMTX = (ytimesxtimesBtimesF./e)';
lossMTX(isnan(lossMTX)) = 0;
lossVEC = sum(lossMTX,2);

% Loss to maintenance
lConsumerMaintenance = bx.*x.*bConsumerGuilds;

% ----------------------------------------------------------------------- %
% PRODUCER BIOMASS DERIVATIVES
% ----------------------------------------------------------------------- %

% PRODUCER gain from logistic growth
G = 1 - (c'*bProducerGuilds + bProducerGuilds)/K;
gProducerGrowth = r.*bProducerGuilds.*G;

% PRODUCER loss to consumption
lProducerConsumption = lossVEC(GuildInfo.iProducerGuildsInFood_position_vector);

lproducerHarvested = ProducerFisheries;

%%%%%%%%%%%%%%%the original below%%%%%%%%%%%%%
dBdt(GuildInfo.iProducerGuilds) = + gProducerGrowth ...
                                  - lProducerConsumption...
                                  - lproducerHarvested;
                                  
plankton = dBdt(95,:);
plankton = plankton + ((61291*18)/100);
dBdt(95,:) = plankton;

% ----------------------------------------------------------------------- %
% CONSUMER BIOMASS DERIVATIVES
% ----------------------------------------------------------------------- %

% CONSUMER maintenance loss
lConsumerMaintenance = lConsumerMaintenance;

% CONSUMER gain from consumption
gConsumerConsumption = gainVEC;


% CONSUMER loss to consumption
NewlConsumerConsumption = lossVEC(GuildInfo.iConsumerGuildsInFood_position_vector); 
xxA = 1:GuildInfo.nGuilds; 
xxB = NewlConsumerConsumption'; 
xxBi = GuildInfo.iConsumerGuildsInFood; 
xxAnew = zeros(1,length(xxA)); 
xxAnew(xxBi) = xxB; 
xxAnew1 = xxAnew'; 
lConsumerConsumption2 = xxAnew1(GuildInfo.iConsumerGuilds); 

%------------------------------------------------------------------------%
lconsumerHarvested = ConsumerFisheries;

dBdt(GuildInfo.iConsumerGuilds) = - lConsumerMaintenance ...
                                  + gConsumerConsumption ...
                                  - lConsumerConsumption2 ...
                                  -lconsumerHarvested;

% ----------------------------------------------------------------------- %
% CONCATENATE DERIVATIVES OF BIOMASSES, CATCHES, GAIN FISH, GAIN AND LOSS
% ----------------------------------------------------------------------- %
dXdt = dBdt;



