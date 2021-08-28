%%%function S = selective_fishery(Data, gearStr)


%%%nAFG = Data.GuildInfo.nAdultFishGuilds;
%%%S = zeros(nAFG,1);

%%%for i = 1:nAFG
    %%%g = Data.Guilds(Data.GuildInfo.iAdultFishGuilds(i));
    
    %%%if g.catchable
       %%% switch lower(gearStr)
             %%%  case 'small-mesh net'
             %%%   mu = 20;
                %%%sigma = 3;
                %%%s = exp(-(g.avgl-mu)^2/(2*sigma^2));
            %%%case 'large-mesh net'
               %%% mu = 30;
                %%%sigma = 4;
                %%%s = exp(-(g.avgl-mu)^2/(2*sigma^2));
            %%%case 'fyke net'
                %a = 10;
                %b = 0.5;
                %s = exp(a+b*g.avgl)/(1+exp(a+b*g.avgl));
               %%% m= - 0.000012; Boit
                %%%s = m; Boit
            %%%case 'hatchery-test'
               %%% Stest = [0.25 0.5 0.75]/0.75;
                %%%s = Stest(g.age-1);
            %%%case 'stochasticity-simulations'
               %%% a50 = 3;
               %%% if g.age < 4
                  %%%  s = 1/(1 + exp(-2*(g.age - a50)));
                %%%else
                   %%% s = 1;
                %%%end
            %%%otherwise
               %%% error('Unknown gear type')
        %%%end
        %%%S(i) = s*Data.hmax;
    %%%else
       %%% S(i) = 0;
    %%%end
%%%end

%%% OLD STUFF
% cycle = 0;
% nAFG = Data.GuildInfo.nAdultFishGuilds;
% S = zeros(nAFG,1);
% for i = 1:nAFG
%     g = Data.Guilds(Data.GuildInfo.iAdultFishGuilds(i));
%     switch lower(gearStr)
%         case 'small-mesh net'
%             a50 = 2*(1+0.0063)^cycle;
%             S(i) = Data.hmax/(1 + exp(-2*(g.age - a50)));
%         case 'large-mesh net'
%             a50 = 3*(1+0.0063)^cycle;
%             S(i) = Data.hmax/(1 + exp(-2*(g.age - a50)));
%         case 'fyke net'
%             a50 = 1*(1+0.0063)^cycle;
%             S(i) = Data.hmax/(1 + exp(-2*(g.age - a50)));
%         otherwise
%             error('Unknown gear type')
%     end
%     if cycle == 0 && g.age == 4
%         S(i) = Data.hmax;
%     end
% end