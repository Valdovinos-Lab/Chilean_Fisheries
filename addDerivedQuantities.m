function Data = addDerivedQuantities(Data)
% calculate relative prey preference matrix
A = Data.communityMatrix;
ng = size(A,1);
nps = 1./sum(A,2);
nps(isinf(nps)) = 0;
w = repmat(nps,1,ng);
w = A.*w;
Data.omega = w;

% %% THESE HAVE BEEN MOVED ELSEWHERE
% % calculate intrinsic growth rate and metabolic rate vectors
% GuildInfo = Data.GuildInfo;
% r = zeros(GuildInfo.nGuilds,1);
% r(GuildInfo.iProducerGuilds) = [Data.Guilds(GuildInfo.iProducerGuilds).igr];
% Data.r = r;
% 
% % metabolic rate calculated from avgl and l-w params if available
% x = zeros(GuildInfo.nGuilds,1);
% x(GuildInfo.iConsumerGuilds) = [Data.Guilds(GuildInfo.iConsumerGuilds).mbr];
% LW_A = [Data.Guilds(GuildInfo.iFishGuilds).lw_a];
% LW_B = [Data.Guilds(GuildInfo.iFishGuilds).lw_b];
% AVGL = [Data.Guilds(GuildInfo.iFishGuilds).avgl];
% % WCC = LW_A.*AVGL.^LW_B*0.2*0.53*1e6;
% % x(GuildInfo.iFishGuilds) = 0.88*(6.4e-5./WCC).^0.11;
% x(GuildInfo.iFishGuilds) = metabolicRate(AVGL,LW_A,LW_B);
% Data.x = x;