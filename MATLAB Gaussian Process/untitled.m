% clear b11 b12 b21 b22
% a(1:2,1:2,1:2,1:2) = 0;
% a(1,1,:,:) = [1 2; 3 4];
% a(1,2,:,:) = [1 2; 3 4];
% a(2,1,:,:) = [1 2; 3 4];
% a(2,2,:,:) = [1 2; 3 4];
% %a = permute(a,[4 3 2 1]) 
% b = mmat(a,a);
% %a = permute(a,[4 3 2 1]);
% b11(1:2,1:2) = b(1,1,:,:)
% b12(1:2,1:2) = b(1,2,:,:)
% b21(1:2,1:2) = b(2,1,:,:)
% b22(1:2,1:2) = b(2,2,:,:)
hyp.lx = [1 2 3];
xDist.cov = [4 5 7; 4 5 6; 7 8 9];
xDist.mean = [1;2;2]
sonig.hyp.ly = 5;
nu = 3;
Xu = [3 5 7; 9 1 1; 1 1 1];
% % tic
% % for i = 1:100
A = diag(1./hyp.lx(:).^2)
B = A + inv(xDist.cov)
C =  det(xDist.cov)*det(B)
D = sqrt(C)
E = sonig.hyp.ly^2 * inv(D)
% 
diff = Xu - repmat(xDist.mean, [1,nu])
F = hyp.lx(:).^ 2
G = diag(F) + xDist.cov
H = diff'*inv(G)
I = H.*diff'
J = sum(I,2)*-1/2
K = exp(J)
L = E*K
% % % end
% % % toc
diffNormalized = repmat(diff,[1,1,1])./repmat(permute(hyp.lx.^2, [1,3,2]),[1,nu,1])
xn = repmat(permute(diffNormalized(:,:,1),[1,3,2,4]),[1,1,1,nu]) + repmat(permute(diffNormalized(:,:,1),[1,3,4,2]),[1,1,nu,1])
A = 1./hyp.lx(:).^2
B = 2*diag(A) + inv(xDist.cov)
C = det(B)
D = det(xDist.cov) * C
E = sqrt(D)
F = sonig.hyp.ly^4 * inv(E)

G = sum (diffNormalized(:,:,1).*diff,1)*-1/2
H = G'*G

I = 1./hyp.lx(:).^2 
		J = 2*diag(I)
		K = J + inv(xDist.cov);
		L = inv(K)
		M = repmat(L,[1,1,nu,nu])
        
		N = mmat( permute(xn,[2,1,3,4]) , M)
		O = mmat(N,xn)
		P = exp(O/2)
		Q = permute(P,[3,4,1,2])
%         
R = F*H.*Q

fMean = L'/Xu*L
fCov = sonig.hyp.ly^2 - trace(Xu\((Xu - Xu/2)/Xu*R))
