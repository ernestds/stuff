function [fDist] = makeSonigStochasticPrediction(sonig, xDist)
%makeSonigStochasticPrediction will predict the output value for a given stochastic trial input point. (Only one point is allowed.)
% This function will apply GP regression, predicting the function output value for given trial input points. The input should be the following.
%	sonig: a SONIG object with (hopefully) already measurements implemented into it.
%	xDist: a trial input point distribution for which to predict the output. It should just be a distribution of size sonig.dx.
%
% The output of the function is subsequently given by a single parameter.
%	fDist: a distribution object giving the posterior distribution of the output. This distribution is of size sonig.dy.

% We check the trial input point that we have been given.
dx = getDistributionSize(xDist);
if dx ~= sonig.dx
    error(['The makeSonigStochasticPrediction function was called with a point of size ',num2str(dx),', while the given SONIG object has points of size ',num2str(sonig.dx),'.']);
end
% i = 1;
% fid = fopen('xDist.mean','wt');
% fprintf(fid,'%g\t',size(xDist.mean,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(xDist.mean,2));
% fprintf(fid,'\n');
% for ii = 1:size(xDist.mean,1)
%     fprintf(fid,'%g\t',xDist.mean(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% fid = fopen('xDist.cov','wt');
% fprintf(fid,'%g\t',size(xDist.cov,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(xDist.cov,2));
% fprintf(fid,'\n');
% for ii = 1:size(xDist.cov,1)
%     fprintf(fid,'%g\t',xDist.cov(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('sonig.Xu.','wt');
% fprintf(fid,'%g\t',size(sonig.Xu,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonig.Xu,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonig.Xu,1)
%     fprintf(fid,'%g\t',sonig.Xu(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% fid = fopen('sonig.hyp.lx','wt');
% fprintf(fid,'%g\t',size(sonig.hyp.lx,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonig.hyp.lx,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonig.hyp.lx,1)
%     fprintf(fid,'%g\t',sonig.hyp.lx(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('sonig.hyp.ly','wt');
% fprintf(fid,'%g\t',size(sonig.hyp.ly,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonig.hyp.ly,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonig.hyp.ly,1)
%     fprintf(fid,'%g\t',sonig.hyp.ly(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('fu_mean','wt');
% fprintf(fid,'%g\t',size(sonig.fu{i}.mean,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonig.fu{i}.mean,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonig.fu{i}.mean,1)
%     fprintf(fid,'%g\t',sonig.fu{i}.mean(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('fu_cov','wt');
% fprintf(fid,'%g\t',size(sonig.fu{i}.cov,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonig.fu{i}.cov,2));
% fprintf(fid,'\n');
% for ii = 1:size(sonig.fu{i}.cov,1)
%     fprintf(fid,'%g\t',sonig.fu{i}.cov(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('Kuu','wt');
% fprintf(fid,'%g\t',size(sonig.Kuu{i},1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(sonig.Kuu{i},2));
% fprintf(fid,'\n');
% for ii = 1:size(sonig.Kuu{i},1)
%     fprintf(fid,'%g\t',sonig.Kuu{i}(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 

% fMean(i,:) = q(:,i)'/sonig.Kuu{i}*sonig.fu{i}.mean;
%    b= trace(sonig.Kuu{i}\(sonig.Kuu{i} - sonig.fu{i}.cov)/sonig.Kuu{i}*Q(:,:,i,i))
% We set up some helpful vectors/matrices.
q = zeros(sonig.nu,sonig.dy);
% Q = zeros(sonig.nu,sonig.nu,sonig.dy,sonig.dy);
% 
% q = gpuArray(q);
% Q = gpuArray(Q);
% % a = sonig.Xu
% % b = xDist.mean
% % c =sonig.nu
% sonig.Xu = gpuArray(sonig.Xu);
% xDist.mean = gpuArray(xDist.mean);
% sonig.nu = gpuArray(sonig.nu);
% sonig.dy = gpuArray(sonig.dy);
% sonig.hyp.lx = gpuArray(sonig.hyp.lx);
% sonig.hyp.ly = gpuArray(sonig.hyp.ly);
% xDist.cov = gpuArray(xDist.cov);
% 
% 
% sonig.Kuu{1} = gpuArray(sonig.Kuu{1});
% sonig.fu{1}.mean = gpuArray(sonig.fu{1}.mean);
% sonig.fu{1}.cov = gpuArray(sonig.fu{1}.cov);
diff = sonig.Xu - repmat(xDist.mean, [1,sonig.nu]); % This is the difference matrix for the trialInput mean with Xu.
diffNormalized = repmat(diff,[1,1,sonig.dy])./repmat(permute(sonig.hyp.lx.^2, [1,3,2]),[1,sonig.nu,1]); %
for i = 1:sonig.dy
    q(:,i) = sonig.hyp.ly(i)^2/sqrt(det(xDist.cov)*det(inv(xDist.cov) + diag(1./sonig.hyp.lx(:).^2)))*exp(-1/2*sum((diff'/(diag(sonig.hyp.lx(:).^2) + xDist.cov)).*diff', 2));
    for j = 1:sonig.dy
        xn = repmat(permute(diffNormalized(:,:,i),[1,3,2,4]),[1,1,1,sonig.nu]) + repmat(permute(diffNormalized(:,:,j),[1,3,4,2]),[1,1,sonig.nu,1]);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
        Q(:,:,i,j) = sonig.hyp.ly(i)^2*sonig.hyp.ly(j)^2/sqrt(det(xDist.cov)*det(inv(xDist.cov) + diag(1./sonig.hyp.lx(:).^2) + diag(1./sonig.hyp.lx(:).^2)))*(exp(-1/2*sum(diffNormalized(:,:,i).*diff,1))'*exp(-1/2*sum(diffNormalized(:,:,j).*diff,1))).*permute(exp(1/2*mmat(mmat(permute(xn,[2,1,3,4]), repmat(inv(inv(xDist.cov) + diag(1./sonig.hyp.lx(:).^2) + diag(1./sonig.hyp.lx(:).^2)), [1,1,sonig.nu,sonig.nu])), xn)),[3,4,1,2]);
        
        mmat(permute(xn,[2,1,3,4]), repmat(inv(inv(xDist.cov) + diag(1./sonig.hyp.lx(:).^2) + diag(1./sonig.hyp.lx(:).^2)), [1,1,sonig.nu,sonig.nu]));
                                                                                                                                                                                                                                                                                                                                                       
    end
end
% A = 1./sonig.hyp.lx(:).^2;
% B = 2*diag(A) + inv(xDist.cov);
% C = det(B);
% D = det(xDist.cov) * C;
% E = sqrt(D);
% F = sonig.hyp.ly(i)^4 * inv(E);
% %sf = size(F)
% G = sum (diffNormalized(:,:,i).*diff,1)*-1/2;
% %sg = size(G)
% H = G'*G;
% 
% 		I = 1./sonig.hyp.lx(:).^2 ;
% 		J = 2*diag(I);
% 		K = J + inv(xDist.cov);
% 		L = inv(K);
%         %size(L)
% 		M = repmat(L,[1,1,sonig.nu,sonig.nu]);
%         %teste(1:size(M,3),1:size(M,4)) = M(2,2,:,:)
%         
%         %size(M)
% 		N = mmat( permute(xn,[2,1,3,4]) , M);
%         
%         %size(N)
% 		O = mmat(N,xn);
%         %teste(1:size(O,3),1:size(O,4)) = O(1,1,:,:)
%         %size(O)
% 		P = exp(O/2);
%         teste(1:size(P,3),1:size(P,4)) = P(1,1,:,:);
% 		Qq = permute(O,[3,4,1,2]);
%         %size(Q)
% 
% 		R = F*H.*Qq;
% %         size(R)
% We calculate the posterior distribution of the output.
fMean = zeros(sonig.dy,1);
fCov = zeros(sonig.dy,sonig.dy);
for i = 1:sonig.dy
    fMean(i,:) = (q(:,i)'/sonig.Kuu{i}*sonig.fu{i}.mean);
    fCov(i,i) = (sonig.hyp.ly(i)^2 - trace(sonig.Kuu{i}\(sonig.Kuu{i} - sonig.fu{i}.cov)/sonig.Kuu{i}*Q(:,:,i,i)));

end
for i = 1:sonig.dy
    for j = 1:sonig.dy
        fCov(i,j) = fCov(i,j) + (sonig.fu{i}.mean'/sonig.Kuu{i}*Q(:,:,i,j)/sonig.Kuu{j}*sonig.fu{j}.mean - fMean(i,:)*fMean(j,:));
    end
end
% We walk through the diagonal elements of the covariance matrix to make sure they're not negative. This is an extra check, because sometimes Matlab has numerical issues which cause diagonal
% elements of the fCov matrix to be negative. Yes, it's a crude fix, but it helps prevent a few of the problems.
for i = 1:sonig.dy
    fCov(i,i) = max(fCov(i,i),1e-16);
end
fDist = createDistribution(fMean, fCov);

% fid = fopen('fMean','wt');
% fprintf(fid,'%g\t',size(fMean,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(fMean,2));
% fprintf(fid,'\n');
% for ii = 1:size(fMean(i,:),1)
%     fprintf(fid,'%g\t',fMean(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)
% 
% fid = fopen('fCov','wt');
% fprintf(fid,'%g\t',size(fCov,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(fCov,2));
% fprintf(fid,'\n');
% for ii = 1:size(fCov(i,:),1)
%     fprintf(fid,'%g\t',fCov(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)

end