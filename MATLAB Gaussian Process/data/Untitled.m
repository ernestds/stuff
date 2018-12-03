Xtrain = {};
j=1;
for m = 1:5
    
    for i = 1:size(Traj{m,1},2)
        
        Xtrain{j} = [Traj{m,1}{i}';TrajVel{m,1}{i}'];
        j = j+1;
    end
end

hold all
for i=1:size(Xtrain,2)
    plot(Xtrain{i}(3,:))
end
