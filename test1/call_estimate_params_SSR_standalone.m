clear all
load gba_test1.mat

global o snpLength mm rr

o.plotSSR = 0;
o.usePrior = false;
o.useRconstraint = false;

rr = gba_test1.rr;
mm = gba_test1.mm;

zTarget = gba_test1.zTarget;
train   = gba_test1.train;
wti0    = ones(size(train.m));
mleTest = gba_test1.mle;
    
snippet   = gba_test1.snippet; 
snpLength = gba_test1.snpLength;
nsim = 30; 

[mle] = estimate_params_SSR_standalone(zTarget,[],train,wti0,snippet,1:9,{'Z'},nsim,[],[]);