function nLLY = objectivefn_Stable_Caushy(x0)
global data
sig=x0^2;
pd= makedist('Stable','alpha',1,'beta',0,'gam',sig,'delta',0);
y = pdf(pd,data);
logY=log(y);
nLLY=-sum(logY);
