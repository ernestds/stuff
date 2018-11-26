function [xpred xsd] = differentclass(xselected, vselected, tnow, nsteps)
% based on the data given predicts the trajectory starting at point tnow + 1
% for the next nsteps 
%
%
load('data/sonigXdiff.mat')
load('data/sonigVdiff.mat')


vPDist = createDistribution([vselected(tnow+1)], eye(1)*1e-13);
xPDist = createDistribution([xselected(tnow+1)], eye(1)*1e-13);

j = 1;
xpred = [];
vpred = [];
xsd = [];
vsd = [];

for i=tnow+1:(size(xselected,2)-1)
    
    st = 1e-8;
    
    vDist = createDistribution([i; vPDist.mean], [st 0; 0 vPDist.cov]);
    vPDist = makeSonigStochasticPrediction(sonigV, vDist);
    vpred(j) = vPDist.mean;
    vsd(j) = vPDist.cov; %predicts velocity
    
    xDist = createDistribution([vpred(j); i; xPDist.mean],...
        [0.03 0 0; 0 st 0 ; 0 0  xPDist.cov]);
    %predicts the position based on last prediction of position and
    %velocity
    xPDist = makeSonigStochasticPrediction(sonigX, xDist);
    xpred(j) = xPDist.mean;
    xsd(j) = sqrt(xPDist.cov);
    
    
    j = j+1;
    if j > nsteps
        break
    end
end


