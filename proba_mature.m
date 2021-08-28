function P = proba_mature(cycle)
a50 = 3*(1- 0.005)^cycle;
a = [2 3 4];
sumL = 1 + exp(-3*(a-a50));
P = [0 0 1./sumL]';              % added zero for larvae and juveniles