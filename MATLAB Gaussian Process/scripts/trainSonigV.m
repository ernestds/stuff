load('../data/traindataallclass1.mat')
stepsback=1;stepsahead=1;
X = Xtrain;

xtrain =[]; ttrain  =[];vtrain =[];ytrain =[]; endtraj = [];vtrainshift = [];
ff = 2;
for i = 1:2:size(X,2)
    ttrain = [ttrain; [1:size(X{i},2)]'];
    xtrain = [xtrain; [X{i}(ff,:)]'];
    vtrain = [vtrain; [X{i}(ff+3,:)]'];
    vtrainshift = [vtrainshift; [X{i}(ff+3,2:end) X{i}(ff+3,end)]'];
    endtraj = [endtraj; zeros(size([X{i}(ff,:)]'))];
    endtraj((end-1)) = 1;
end
endtraj((end)) = 0;

%try to tune hyperparameters using a NIGP model
xm = [ ttrain'; vtrain';];
ym = vtrainshift';
evalc('[model, nigp] = trainNIGP(permute(xm(:,1:end),[2,1]),permute(ym(:,1:end),[2,1]),-500,1);')
hyp = NIGPModelToHyperparameters(model)

%%
load('../data/traindataallclass1.mat')
stepsback=1;stepsahead=1;
X = Xtrain;

xtrain =[]; ttrain  =[];vtrain =[];ytrain =[]; endtraj = [];vtrainshift = [];
ff = 3;
for i = 1:2:size(X,2)
    ttrain = [ttrain; [1:size(X{i},2)]'];
    xtrain = [xtrain; [X{i}(ff,:)]'];
    vtrain = [vtrain; [X{i}(ff+3,:)]'];
    vtrainshift = [vtrainshift; [X{i}(ff+3,2:end) X{i}(ff+3,end)]'];
    endtraj = [endtraj; zeros(size([X{i}(ff,:)]'))];
    endtraj((end-1)) = 1;
end
endtraj((end)) = 0;
distance = 0.4;
ly = 0.6;
lv = ly;
la = 6;
lt = 60;
sy = ly * .03;
sv = lv *.001;
sa = la *.01;
st = lt *.01;
features = [];
feat = '';
hyp.lx = [];
hyp.sx = [];

prevTSize = 0;
stepsbackT = 1;
for i = (stepsbackT-1):-1:0
    features =  [features; ttrain(stepsbackT - i);];
    thing = "(i - " + num2str(i) + ")";
    feat = strcat(feat,"ttrain" + thing + ';');
    hyp.lx = [hyp.lx; lt];
    hyp.sx = [hyp.sx; st];
    prevTSize = prevTSize+1;
end
stepsbackA = 0;
prevASize = 0;
for i = (stepsbackA-1):-1:0
    features =  [features; atrain(i + 1);];
    thing = "(i - " + num2str(i) + ")";
    feat = strcat(feat,"atrain" + thing + ';');
    hyp.lx = [hyp.lx; la];
    hyp.sx = [hyp.sx; sa];
    prevASize = prevASize+1;
end

prevXSize = 0;

stepsbackX = stepsback;
for i = 1:stepsbackX
    features =  [features; vtrain(i)];
    thing = "(size(features,1) - stepsback + 1 +" + num2str(i) + ")";
    feat = strcat(feat,"jointMean" + thing + ';');
    hyp.lx = [hyp.lx; ly];
    hyp.sx = [hyp.sx; sy];
    prevXSize = prevXSize+1;
end


feat = "features = [" + feat + '];';

s = hyp.sx.^2;

hyp.ly = [];
hyp.sy = [];
outputs = [];
so = [];

for j = (stepsback):(stepsahead + stepsback - 1)
    outputs = [outputs; vtrain(j + 1)];
    so = [so sy^2];
    hyp.ly = [hyp.ly; ly];
    hyp.sy  = [hyp.sy ; sy];
end

hyp.sy = sy;
% selxu = randperm(size(xtrain,1),40);
% xu = [vtrain(selxu)'; ttrain(selxu)'; xtrain(selxu)'];


hyp.lx = [25;0.8]*1;
hyp.sx = [4.57630004290503;0.0052451756094707]*01;
hyp.ly = 0.560919610368498*1;
hyp.sy = 0.00155019635324873;
sonigV = createSONIG(hyp);
sonigV.addIIPDistance = distance;
jointMean = [features; outputs];jointCov = [eye(size(features,1)).*s zeros( size(features,1) , size(outputs,1));...
    zeros( size(outputs,1) , size(features,1)) eye(size(outputs,1)) .* so];

jointDist = createDistribution(jointMean, jointCov); % We create a distribution from the mean and the covariance. This is the joint distribution for the SONIG input and SONIG output.
inputDist = getSubDistribution(jointDist, [1:size(features,1)]); % We extract the SONIG input, which is evidently required by the SONIG algorithm.
outputDist = getSubDistribution(jointDist, (size(features,1)+1):(size(features,1)+size(outputs,1)) ); % We extract the SONIG output, which is also required.

[sonigV, inputPost, outputPost, jointPost] = implementMeasurement(sonigV, inputDist, outputDist, jointDist);

jointMean = jointPost.mean; % We will need this one for the next iteration.
jointCov = jointPost.cov;


a = [];
i = stepsback + 1;
for count = (stepsback + 1):(size(vtrain,1)-stepsahead-1)
    tic
    
    if endtraj(i)
        i = i + stepsahead + 1;
    else
        i = i + 1;
    end
    if (i+stepsahead) > size(xtrain,1)
        break;
    end
    [size(sonigV.Xu,2) i]
    
    outputs = [];
    for j = 1:(stepsahead)
        outputs = [outputs; vtrain(j+i) ];
    end

    eval(feat);
    a(count,:) = [vtrain(i) outputs(end)];
    
    jointMean = [features; outputs];

    
    %
    jointCov = [jointCov(2:end,2:end) zeros(size(jointCov(2:end,2:end),1),1);...
        zeros(size(jointCov(2:end,2:end),1),1)' sy^2];
    jointCov(1:size(features,1)-prevXSize,:) = 0;
    jointCov(:,1:size(features,1)-prevXSize) = 0;
    jointCov(1:size(features,1)-prevXSize, 1:size(features,1)-prevXSize) = eye(size(features,1)-prevXSize).*s(1:end-prevXSize);
    
    
  
    jointDist = createDistribution(jointMean, jointCov); % We create a distribution from the mean and the covariance. This is the joint distribution for the SONIG input and SONIG output.
    inputDist = getSubDistribution(jointDist, [1:size(features,1)]); % We extract the SONIG input, which is evidently required by the SONIG algorithm.
    outputDist = getSubDistribution(jointDist, (size(features,1)+1):(size(features,1)+size(outputs,1)) ); % We extract the SONIG output, which is also required.
    [sonigV, inputPost, outputPost, jointPost] = implementMeasurement(sonigV, inputDist, outputDist, jointDist);
    %     We update the distributions of all our points.
    jointMean = jointPost.mean; % We will need this one for the next iteration.
    jointCov = jointPost.cov;
    toc
    %a(i,:) = outputDist.mean';
    %
end

which = 2;
load('../data/testdataclass1.mat')
xselected = Xpred{which}(ff,:);
xselectedshift = [Xpred{which}(ff,1) Xpred{which}(ff,1:(end-1))];
vselectedshift = [Xpred{which}(ff+3,1) Xpred{which}(ff+3,1:(end-1))];
vselected = Xpred{which}(ff+3,:);
aselected = Xpred{which}(ff+6,:);
tselected = 1:size(xselected,2);

figure(2);

h6 = plot(vselected);

xmarked = [];xshiftmarked = [];
tmarked = [];
vmarked = [];
amarked = [];
for iter = (stepsback+1):5:(size(xselected,2)-1)
    tic
    xmarked = [xmarked; xselected(iter)];
    tmarked = [tmarked; tselected(iter)];
    xshiftmarked = [xshiftmarked; xselectedshift(iter)];
    vmarked = [vmarked; vselected(iter)];
    amarked = [amarked; aselected(iter)];
    
    aPDist = createDistribution([aselected(iter)], eye(1)*1e-13);
    vPDist = createDistribution([vselected(iter)], eye(1)*1e-13);
    xPDist = createDistribution([xselected(iter)], eye(1)*1e-13);
    
    j = 1;
    xshift = xselected(iter);
    xpred = 0;
    vpred = [];
    apred = [];
    xsd = [];
    vsd = [];
    asd = [];
    

    
    for i=iter+1:(size(xselected,2)-1)
        
        
        
        t(j) = i;
        
        st = 1e-8;
        
        vDist = createDistribution([t(j); vPDist.mean], [st 0; 0 vPDist.cov]);
        [vPDist] = makeSonigStochasticPrediction(sonigV, vDist);
        vpred(j) = vPDist.mean;
        vsd(j) = vPDist.cov;
        
        xDist = createDistribution([vpred(j); tselected(i); xPDist.mean],...
            [0.03 0 0; 0 st 0 ; 0 0  xPDist.cov]);
        
        j = j+1;
        if j > 30
            break
        end
    end
    toc
    %
    
    
    if iter>stepsback+1 %delete old plots
        delete(h7)
        delete(h22)
        delete(h33)
    end
    
    figure(2);
    grid on;
    hold on;
    h7 = plot(t,vpred);
    hold on;
    h22=plot(t, vpred + 2*vsd, 'g-');
    hold on;
    h33=plot(t, vpred - 2*vsd, 'g-');

    
    drawnow;
    waitforbuttonpress
    
    
end





