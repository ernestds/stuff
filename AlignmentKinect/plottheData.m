load('../TrajectoriesDataSameCoordinateSystem.mat');
%figure;
counter(1:6) = 0;
avg(1:6,1:3) = 0;
for m = 1:5
    for n = [1 3 5 7 9 11 13 15 17]
        for i = 1:length(Traj{m,n})
            for k = 1:1 ... length(Traj{m,n}{i})
                    %                 M = Rotation_ZYX(TrajRot{m,n}{i}(k,3),TrajRot{m,n}{i}(k,2),TrajRot{m,n}{i}(k,1));
                %                 x_vec = M*[0.05;0;0];
                %                 y_vec = M*[0;0.05;0];
                %                 z_vec = M*[0;0;0.05];
                %                 quiver3(Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3),x_vec(1),x_vec(2),x_vec(3),'b');
                %                 hold on
                %                 quiver3(Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3),y_vec(1),y_vec(2),y_vec(3),'r');
                %                 hold on
                %                 quiver3(Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3),z_vec(1),z_vec(2),z_vec(3),'g');
                %                 hold on
                %                 scatter3(Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3),'m')
                %                 hold on
                %plot3(Traj{m,n}{i}([1 end],1),Traj{m,n}{i}([1 end],2),Traj{m,n}{i}([1 end],3), '.')
                %hold on
                if sum(n==[1 3 5])>0
                    counter(1) = counter(1) + 1;
                    avg(1,:) = avg(1,:) + [Traj{m,n}{i}(k,1) Traj{m,n}{i}(k,2) Traj{m,n}{i}(k,3)];
                end
                if sum(n==[7 9 11])>0
                    counter(2) = counter(2) + 1;
                    avg(2,:) = avg(2,:) + [Traj{m,n}{i}(k,1) Traj{m,n}{i}(k,2) Traj{m,n}{i}(k,3)];
                end
                if sum(n==[13 15 17])>0
                    counter(3) = counter(3) + 1;
                    avg(3,:) = avg(3,:) + [Traj{m,n}{i}(k,1) Traj{m,n}{i}(k,2) Traj{m,n}{i}(k,3)];
                end
                if sum(n==[1 7 13])>0 %ending point 1
                    counter(4) = counter(4) + 1;
                    avg(4,:) = avg(4,:) + [Traj{m,n}{i}(end,1) Traj{m,n}{i}(end,2) Traj{m,n}{i}(end,3)];
                end
                if sum(n==[3 9 15])>0 %ending point 1
                    counter(5) = counter(5) + 1;
                    avg(5,:) = avg(5,:) + [Traj{m,n}{i}(end,1) Traj{m,n}{i}(end,2) Traj{m,n}{i}(end,3)];
                end
                if sum(n==[5 11 17])>0 %ending point 1
                    counter(6) = counter(6) + 1;
                    avg(6,:) = avg(6,:) + [Traj{m,n}{i}(end,1) Traj{m,n}{i}(end,2) Traj{m,n}{i}(end,3)];
                end
            end
        end
    end
end

avg(1,:) = avg(1,:)/counter(1);
avg(2,:) = avg(2,:)/counter(2);
avg(3,:) = avg(3,:)/counter(3);
avg(4,:) = avg(4,:)/counter(4);
avg(5,:) = avg(5,:)/counter(5);
avg(6,:) = avg(6,:)/counter(6);

distto1 = [];
distto2 = [];
distto3 = [];
distto4 = [];
distto5 = [];
distto6 = [];

a = 0; b = 0; c = 0; d = 0; e = 0; f = 0;
for m = 1:5
    for n = [1 3 5 7 9 11 13 15 17]
        for i = 1:length(Traj{m,n})
            for k = 1:1 ... length(Traj{m,n}{i})
                    %                 M = Rotation_ZYX(TrajRot{m,n}{i}(k,3),TrajRot{m,n}{i}(k,2),TrajRot{m,n}{i}(k,1));
                %                 x_vec = M*[0.05;0;0];
                %                 y_vec = M*[0;0.05;0];
                %                 z_vec = M*[0;0;0.05];
                %                 quiver3(Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3),x_vec(1),x_vec(2),x_vec(3),'b');
                %                 hold on
                %                 quiver3(Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3),y_vec(1),y_vec(2),y_vec(3),'r');
                %                 hold on
                %                 quiver3(Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3),z_vec(1),z_vec(2),z_vec(3),'g');
                %                 hold on
                %                 scatter3(Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3),'m')
                %                 hold on
                %plot3(Traj{m,n}{i}([1 end],1),Traj{m,n}{i}([1 end],2),Traj{m,n}{i}([1 end],3), '.')
                %hold on
                if sum(n==[1 3 5])>0
                    a = a+ 1;
                    distto1(a,:) = [norm( avg(1,:) - [ Traj{m,n}{i}(k,1) ,Traj{m,n}{i}(k,2) ,Traj{m,n}{i}(k,3) ] ) Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3) TrajRot{m,n}{i}(k,1),TrajRot{m,n}{i}(k,2),TrajRot{m,n}{i}(k,3)];
                end
                if sum(n==[7 9 11])>0
                    b = b+ 1;
                    distto2(b,:) = [norm( avg(2,:) - [ Traj{m,n}{i}(k,1) ,Traj{m,n}{i}(k,2) ,Traj{m,n}{i}(k,3) ] ) Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3) TrajRot{m,n}{i}(k,1),TrajRot{m,n}{i}(k,2),TrajRot{m,n}{i}(k,3)];
                end
                if sum(n==[13 15 17])>0
                    c = c+ 1;
                    distto3(c,:) = [norm( avg(3,:) - [ Traj{m,n}{i}(k,1) ,Traj{m,n}{i}(k,2) ,Traj{m,n}{i}(k,3) ] ) Traj{m,n}{i}(k,1),Traj{m,n}{i}(k,2),Traj{m,n}{i}(k,3) TrajRot{m,n}{i}(k,1),TrajRot{m,n}{i}(k,2),TrajRot{m,n}{i}(k,3)];
                end
                if sum(n==[1 7 13])>0
                    d = d+ 1;
                    distto4(d,:) = [norm( avg(4,:) - [ Traj{m,n}{i}(end,1) ,Traj{m,n}{i}(end,2) ,Traj{m,n}{i}(end,3) ] ) Traj{m,n}{i}(end,1),Traj{m,n}{i}(end,2),Traj{m,n}{i}(end,3) TrajRot{m,n}{i}(end,1),TrajRot{m,n}{i}(end,2),TrajRot{m,n}{i}(end,3)];
                end
                if sum(n==[3 9 15])>0
                    e = e+ 1;
                    distto5(e,:) = [norm( avg(5,:) - [ Traj{m,n}{i}(end,1) ,Traj{m,n}{i}(end,2) ,Traj{m,n}{i}(end,3) ] ) Traj{m,n}{i}(end,1),Traj{m,n}{i}(end,2),Traj{m,n}{i}(end,3) TrajRot{m,n}{i}(end,1),TrajRot{m,n}{i}(end,2),TrajRot{m,n}{i}(end,3)];
                end
                if sum(n==[5 11 17])>0
                    f = f+ 1;
                    distto6(f,:) = [norm( avg(6,:) - [ Traj{m,n}{i}(end,1) ,Traj{m,n}{i}(end,2) ,Traj{m,n}{i}(end,3) ] ) Traj{m,n}{i}(end,1),Traj{m,n}{i}(end,2),Traj{m,n}{i}(end,3) TrajRot{m,n}{i}(end,1),TrajRot{m,n}{i}(end,2),TrajRot{m,n}{i}(end,3)];
                end
            end
        end
        
        
    end
end
[B I] = sort(distto1);
distto1 = distto1(I(:,1),:);
[B I] = sort(distto2);
distto2 = distto2(I(:,1),:);
[B I] = sort(distto3);
distto3 = distto3(I(:,1),:);
[B I] = sort(distto4);
distto4 = distto4(I(:,1),:);
[B I] = sort(distto5);
distto5 = distto5(I(:,1),:);
[B I] = sort(distto6);
distto6 = distto6(I(:,1),:);

points{1} = distto1(1:3,:);
points{2} = distto2(1:3,:);
points{3} = distto3(1:3,:);
points{4} = distto4(1:3,:);
points{5} = distto5(1:3,:);
points{6} = distto6(1:3,:);
%plot3([distto1(:,2);distto2(:,2);distto3(:,2)],[distto1(:,3);distto2(:,3);distto3(:,3)],[distto1(:,4);distto2(:,4);distto3(:,4)],'b.')
hold on
grid
plot3([avg(:,1);points{1}(:,2);points{2}(:,2);points{3}(:,2)],[avg(:,2);points{1}(:,3);points{2}(:,3);points{3}(:,3)],[avg(:,3);points{1}(:,4);points{2}(:,4);points{3}(:,4)],'r.')
pointsQuart = [];
for i=1:6
    for j=1:3
        R_ZYX = Rotation_ZYX(points{i}(j,7), points{i}(j,6), points{i}(j,5));
        
        x_vec = R_ZYX*[0.02;0;0];
        y_vec = R_ZYX*[0;0.02;0];
        z_vec = R_ZYX*[0;0;0.2];
        
        quiver3(points{i}(j,2),points{i}(j,3),points{i}(j,4),x_vec(1),x_vec(2),x_vec(3),'b');
        hold on
        quiver3(points{i}(j,2),points{i}(j,3),points{i}(j,4),y_vec(1),y_vec(2),y_vec(3),'r');
        hold on
        quiver3(points{i}(j,2),points{i}(j,3),points{i}(j,4),z_vec(1),z_vec(2),z_vec(3),'g');
        hold on
        pointsQuart{i}(j,:) = [points{i}(j,2),points{i}(j,3),points{i}(j,4),eulerToQuartenion([points{i}(j,5), points{i}(j,6), points{i}(j,7)])];
        
    end
end

legend('Point','X axis', 'Y axis', 'Z axis')