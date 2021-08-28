function Results = initializeResultStruct(GuildInfo)
Results.B = zeros(GuildInfo.nGuilds,0);
%%%Results.C = zeros(GuildInfo.nAdultFishGuilds,0); %%%this comment was ommited
%%%Results.GF = zeros(GuildInfo.nFishGuilds,0);
%%%Results.G = zeros(GuildInfo.nGuilds,GuildInfo.nGuilds,0);
%%%%Results.L = zeros(GuildInfo.nGuilds,GuildInfo.nGuilds,0);
Results.allbiomasses = zeros(GuildInfo.nGuilds,0);
%%%Results.ARE = zeros(GuildInfo.nAdultFishGuilds,0);
