function [xpred xsd] = similiartrajectories(xselected, vselected, tnow, nsteps)
% based on the data given predicts the trajectory starting at point tnow + 1
% for the next nsteps 
%
%
load('data/sonigVsimoptimized.mat')
load('data/sonigXsimoptimized.mat')

% i = 1;
% fid = fopen('trajx2','wt');
% fprintf(fid,'%g\t',size(x,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(x,2));
% fprintf(fid,'\n');
% for ii = 1:size(x,1)
%     fprintf(fid,'%g\t',x(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% fid = fopen('sonigX.hyp.lx','wt');
% fprintf(fid,'%g\t',size(sonigX.hyp.lx,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigX.hyp.lx,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigX.hyp.lx,1)
%     fprintf(fid,'%g\t',sonigX.hyp.lx(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('sonigX.hyp.ly','wt');
% fprintf(fid,'%g\t',size(sonigX.hyp.ly,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigX.hyp.ly,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigX.hyp.ly,1)
%     fprintf(fid,'%g\t',sonigX.hyp.ly(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('Xfu_mean','wt');
% fprintf(fid,'%g\t',size(sonigX.fu{i}.mean,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigX.fu{i}.mean,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigX.fu{i}.mean,1)
%     fprintf(fid,'%g\t',sonigX.fu{i}.mean(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('Xfu_cov','wt');
% fprintf(fid,'%g\t',size(sonigX.fu{i}.cov,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigX.fu{i}.cov,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigX.fu{i}.cov,1)
%     fprintf(fid,'%g\t',sonigX.fu{i}.cov(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('XKuu','wt');
% fprintf(fid,'%g\t',size(sonigX.Kuu{i},1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigX.Kuu{i},2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigX.Kuu{i},1)
%     fprintf(fid,'%g\t',sonigX.Kuu{i}(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('sonigV.Xu.','wt');
% fprintf(fid,'%g\t',size(sonigV.Xu,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigV.Xu,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigV.Xu,1)
%     fprintf(fid,'%g\t',sonigV.Xu(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% fid = fopen('sonigV.hyp.lx','wt');
% fprintf(fid,'%g\t',size(sonigV.hyp.lx,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigV.hyp.lx,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigV.hyp.lx,1)
%     fprintf(fid,'%g\t',sonigV.hyp.lx(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('sonigV.hyp.ly','wt');
% fprintf(fid,'%g\t',size(sonigV.hyp.ly,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigV.hyp.ly,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigV.hyp.ly,1)
%     fprintf(fid,'%g\t',sonigV.hyp.ly(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('Vfu_mean','wt');
% fprintf(fid,'%g\t',size(sonigV.fu{i}.mean,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigV.fu{i}.mean,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigV.fu{i}.mean,1)
%     fprintf(fid,'%g\t',sonigV.fu{i}.mean(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('Vfu_cov','wt');
% fprintf(fid,'%g\t',size(sonigV.fu{i}.cov,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigV.fu{i}.cov,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigV.fu{i}.cov,1)
%     fprintf(fid,'%g\t',sonigV.fu{i}.cov(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('VKuu','wt');
% fprintf(fid,'%g\t',size(sonigV.Kuu{i},1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonigV.Kuu{i},2));
% fprintf(fid,'\n');
% for ii = 1:size(sonigV.Kuu{i},1)
%     fprintf(fid,'%g\t',sonigV.Kuu{i}(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% 
% fid = fopen('trajectoryX','wt');
% fprintf(fid,'%g\t',size(xselected,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(xselected,2));
% fprintf(fid,'\n');
% for ii = 1:size(xselected,1)
%     fprintf(fid,'%g\t',xselected(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('trajectoryV','wt');
% fprintf(fid,'%g\t',size(vselected,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(vselected,2));
% fprintf(fid,'\n');
% for ii = 1:size(vselected,1)
%     fprintf(fid,'%g\t',vselected(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)


vPDist = createDistribution([vselected(tnow+1)], eye(1)*1e-13);
xPDist = createDistribution([xselected(tnow+1)], eye(1)*1e-13);

j = 1;
xpred = [];
vpred = [];
xsd = [];
vsd = [];

total = toc;
for k =tnow+1:(size(xselected,2)-1)
    
    st = 1e-8;
    
    vDist = createDistribution([k; vPDist.mean], [st 0; 0 vPDist.cov]);
    temp = toc;
    vPDist = makeSonigStochasticPrediction(sonigV, vDist);
    predV = toc - temp
    vpred(j) = vPDist.mean;
    vsd(j) = vPDist.cov; %predicts velocity
    
    xDist = createDistribution([vpred(j); k; xPDist.mean],...
        [0.03 0 0; 0 st 0 ; 0 0  xPDist.cov]);
    %predicts the position based on last prediction of position and
    %velocity
    temp = toc;
    xPDist = makeSonigStochasticPrediction(sonigX, xDist);
    predX = toc - temp
    
    xpred(j) = xPDist.mean;
    xsd(j) = sqrt(xPDist.cov);
    
    
    j = j+1;
    if j > nsteps
        break
    end
end
tempotodaspred = toc - total
