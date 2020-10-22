function nLLY = objectivefn_Stable_Levy_xmin(x0)
global data
sig=x0^2;
pd= makedist('Stable','alpha',0.5,'beta',1,'gam',sig,'delta',0);
y = pdf(pd,data);
logY=log(y);
nLLY=-sum(logY);
