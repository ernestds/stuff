load('../data/traindataallclass1.mat')
stepsback=1;stepsahead=1;
X = Xtrain;

xtrain =[]; ttrain  =[];vtrain =[];ytrain =[]; endtraj = [];xtrainshift = [];
for i = 1:size(X,2)
    ttrain = [ttrain; [1:size(X{i},2)]'];
    xtrain = [xtrain; [X{i}(1,:)]'];
    vtrain = [vtrain; [X{i}(4,:)]'];
    xtrainshift = [xtrainshift; [X{i}(1,2:end) X{i}(1,end)]'];
    endtraj = [endtraj; zeros(size([X{i}(1,:)]'))];
    endtraj((end-1)) = 1;
end
endtraj((end)) = 0;

%try to tune hyperparameters using a NIGP model
xm = [vtrain'; ttrain'; xtrain'];
ym = xtrainshift';
evalc('[model, nigp] = trainNIGP(permute(xm(:,1:end),[2,1]),permute(ym(:,1:end),[2,1]),-500,1);')
hyp = NIGPModelToHyperparameters(model)

%%
distance = 0.18;


%parameters that I found through trial and error
lyy = 0.1;
lx = lyy;
lv = 4;
la = 40;
lt = 100;
syy = 0.005;
sv = 0.03
st = 1;
sx = 0.001;
features = [];
feat = '';
lxx = [];
sxx = [];

%this is just the way I did so that I could automatically increase the
%amount of steps included as a feature, as well as the number of predicted
%steps in the case of the batch way
prevVSize = 0;
stepsbackV = 1;
for i = (stepsbackV-1):-1:0
    features =  [features; vtrain(stepsbackV - i);];
    thing = "(i - " + num2str(i) + ")";
    feat = strcat(feat,"vtrain" + thing + ';');
    lxx = [lxx; lv];
    sxx = [sxx; sv];
    prevVSize = prevVSize+1;
end

prevTSize = 0;
stepsbackT = 1;
for i = (stepsbackT-1):-1:0
    features =  [features; ttrain(stepsbackT - i);];
    thing = "(i - " + num2str(i) + ")";
    feat = strcat(feat,"ttrain" + thing + ';');
    lxx = [lxx; lt];
    sxx = [sxx; st];
    prevTSize = prevTSize+1;
end
prevASize = 0;
% for i = (stepsback-1):-1:0
%     features =  [features; atrain(i + 1);];
%     thing = "(i - " + num2str(i) + ")";
%     feat = strcat(feat,"atrain" + thing + ';');
%     lxx = [lxx; la];
%     sxx = [sxx; sa];
%     prevASize = prevASize+1;
% end
prevXSize = 0;
%
stepsbackX = stepsback;
for i = 1:stepsbackX
    features =  [features; xtrain(i)];
    thing = "(size(features,1) - stepsback + 1 +" + num2str(i) + ")";
    feat = strcat(feat,"jointMean" + thing + ';');
    lxx = [lxx; lx];
    sxx = [sxx; sx];
    prevXSize = prevXSize+1;
end

% for i = (stepsback-1):-1:0
%     features =  [features; xtrain(stepsback - i);];
%     thing = "(i - " + num2str(i) + ")";
%     feat = strcat(feat,"xtrain" + thing + ';');
%     lxx = [lxx; lyy];
%     sxx = [sxx; syy];
%     prevXSize = prevXSize+1;
% end
%

feat = "features = [" + feat + '];';

s = sxx.^2;
sxx = sxx*1;

ly = [];
sy = [];
outputs = [];
so = [];

sx = [];
lx = [];

for j = (stepsback):(stepsahead + stepsback - 1)
    outputs = [outputs; xtrain(j + 1)];
    so = [so syy^2];
    ly = [ly; lyy];
    sy = [sy; syy];
end


hyp.lx = lxx;
hyp.sx = sxx;
hyp.ly = ly;
hyp.sy = sy;
hyp.lx = [22.13; 269.5; 0.27];
hyp.lx = [.54; 269.5; 0.27];
hyp.sx = [0.000656451484235225; 15.0827137839194; 5.89558238910542e-05];

hyp.ly = 0.086225582343152;
hyp.sy = 4.72025688279776e-05; %those parameters were found using the optimizer that uses
%the NIGP, as is done in the examples from the creator of SONIG
sonigX = createSONIG(hyp);
sonigX.addIIPDistance = distance;
jointMean = [features; outputs]; % This will be the joint vector of inputs and outputs which we're currently applying. It will basically "shift through time" upon each iteration. The first three entries are for inputs (being u(k-2), u(k-1), u(k)) and the other two entries are for outputs (being y(k) and y(k+1)).
jointCov = [eye(size(features,1)).*s zeros( size(features,1) , size(outputs,1));...
    zeros( size(outputs,1) , size(features,1)) eye(size(outputs,1)) .* so];

jointDist = createDistribution(jointMean, jointCov); % We create a distribution from the mean and the covariance. This is the joint distribution for the SONIG input and SONIG output.
inputDist = getSubDistribution(jointDist, [1:size(features,1)]); % We extract the SONIG input, which is evidently required by the SONIG algorithm.
outputDist = getSubDistribution(jointDist, (size(features,1)+1):(size(features,1)+size(outputs,1)) ); % We extract the SONIG output, which is also required.

[sonigX, inputPost, outputPost, jointPost] = implementMeasurement(sonigX, inputDist, outputDist, jointDist);

jointMean = jointPost.mean; % We will need this one for the next iteration.
jointCov = jointPost.cov;



a = [];
i = stepsback + 1;
for count = (stepsback + 1):(size(xtrain,1)-stepsahead-1)
    tic
    
    if endtraj(i)
        i = i + 2 + 1;
    else
        i = i + 1;
    end
    if (i+stepsahead) > size(xtrain,1)
        break;
    end
    [size(sonigX.Xu,2) i]
    
    outputs = [];
    for j = 1:(stepsahead)
        outputs = [outputs; xtrain(j+i) ];
    end
    
 
    eval(feat); %this is basically
    %features = [vtrain(i - 0);ttrain(i - 0);jointMean(4);]
    
%     a(count,:) = [xtrain(i)+0.1 outputs(end) -endtraj(i)*0.5];
    

    jointMean = [features; outputs];
    
    jointCov = [jointCov(2:end,2:end) zeros(size(jointCov(2:end,2:end),1),1);...
        zeros(size(jointCov(2:end,2:end),1),1)' syy^2];
    jointCov(1:size(features,1)-prevXSize,:) = 0;
    jointCov(:,1:size(features,1)-prevXSize) = 0;
    jointCov(1:size(features,1)-prevXSize, 1:size(features,1)-prevXSize) = eye(size(features,1)-prevXSize).*s(1:end-prevXSize);
    %shifts the covariance matrix by one, that means the last predictions
    %covariance now is the position covariance.. the other inputs
    %disregards the post and replaces it with the specified covariances
    
    
    jointDist = createDistribution(jointMean, jointCov); % We create a distribution from the mean and the covariance. This is the joint distribution for the SONIG input and SONIG output.
    inputDist = getSubDistribution(jointDist, [1:size(features,1)]); % We extract the SONIG input, which is evidently required by the SONIG algorithm.
    outputDist = getSubDistribution(jointDist, (size(features,1)+1):(size(features,1)+size(outputs,1)) ); % We extract the SONIG output, which is also required.
    [sonigX, inputPost, outputPost, jointPost] = implementMeasurement(sonigX, inputDist, outputDist, jointDist);
    %     We update the distributions of all our points.
    jointMean = jointPost.mean; % We will need this one for the next iteration.
    jointCov = jointPost.cov;
    toc
    %a(i,:) = outputDist.mean';
    %
end
%% 
%validate

which = 1;
load('../data/testdataclass1.mat')

xselected = Xpred{which}(1,:);
xselectedshift = [Xpred{which}(1,1) Xpred{which}(1,1:(end-1))];
vselectedshift = [Xpred{which}(4,1) Xpred{which}(4,1:(end-1))];
vselected = Xpred{which}(4,:);
aselected = Xpred{which}(7,:);
tselected = 1:size(xselected,2);

figure(1)
hold all
h5 = plot(tselected, xselected, 'k');
clear xlim
xlim([0  size(xselected,2)]);
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
    %pretends that the first i-1 is a prediction with high certainty
    
    
    j = 1;
    xshift = xselected(iter);
    xpred = 0;
    vpred = [];
    apred = [];
    xsd = [];
    vsd = [];
    asd = [];
    

    %predicts step i based on predicted i-1
    for i=iter+1:(size(xselected,2)-1)
        
        
        
        t(j) = i;
        
        st = 1e-8;
        
        vDist = createDistribution([t(j); vPDist.mean], [st 0; 0 vPDist.cov]);
        [vPDist] = makeSonigStochasticPrediction(sonigV, vDist);
        vpred(j) = vPDist.mean;
        vsd(j) = sqrt(vPDist.cov);
        
        xDist = createDistribution([vpred(j); tselected(i); xPDist.mean],...
            [vsd(j) 0 0; 0 st 0 ; 0 0  xPDist.cov]);
%         xDist = createDistribution([vpred(j); xPDist.mean],...
%             [0.03 0 ; 0 xPDist.cov]);
        [xPDist] = makeSonigStochasticPrediction(sonigX, xDist);
        xpred(j) = xPDist.mean;
        xsd(j) = sqrt(xPDist.cov);
        

        j = j+1;
        if j > 30
            break
        end
    end
    toc
    %
    
    %all below is just plot stuff
    if iter>stepsback+1 %delete old plots
        delete(h1)
        delete(h2)
        delete(h3)
        delete(h4)
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

    
    figure(1);
    tplot = (iter+1):(iter+size(xpred,2));
    grid on;
    hold on;
    
    h1=plot(tplot, xpred, 'b-');
    h2=plot(tplot, xpred + 2*xsd, 'g-');
    h3=plot(tplot, xpred - 2*xsd, 'g-');
    h4=plot(tmarked, xmarked, 'ko');
    xlim = [0  size(xselected,2)];
    
    
    
    drawnow;
    waitforbuttonpress
    
    
end

