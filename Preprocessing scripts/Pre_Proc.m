clear all
close all

%% Load Data
load('MarkusRHandOut.mat');
load('MinRHandOut.mat');
load('SwenRHandOut.mat');
load('YakunRHandOut.mat');
load('YanhaoRHandOut.mat');

load('YanhaoRot.mat');
load('YakunRot.mat');
load('SwenRot.mat');
load('MinRot.mat');
load('MarkusRot.mat');

Data_Teilnehmer_1_X = Skeleton2RHandOut;
Data_Teilnehmer_1_Y = -Skeleton2RHandOut2;
Data_Teilnehmer_1_Z = Skeleton2RHandOut1;

Data_Teilnehmer_R1_X = MarkusRot; %rotation data, same thing as position
% MarkusRot = X, MarkusRot2 = Z, MarkusRot1 = Y
Data_Teilnehmer_R1_Y = -MarkusRot2;
Data_Teilnehmer_R1_Z = MarkusRot1;

clear Skeleton2RHandOut
clear Skeleton2RHandOut1
clear Skeleton2RHandOut2
clear MarkusRot
clear MarkusRot1
clear MarkusRot2

Data_Teilnehmer_2_X = MinRHandOut;
Data_Teilnehmer_2_Y = -MinRHandOut2;
Data_Teilnehmer_2_Z = MinRHandOut1;

Data_Teilnehmer_R2_X = MinRot;
Data_Teilnehmer_R2_Y = -MinRot2;
Data_Teilnehmer_R2_Z = MinRot1;

clear MinRHandOut
clear MinRHandOut1
clear MinRHandOut2
clear MinRot
clear MinRot1
clear MinRot2

Data_Teilnehmer_3_X = SwenRHandOut;
Data_Teilnehmer_3_Y = -SwenRHandOut2;
Data_Teilnehmer_3_Z = SwenRHandOut1;

Data_Teilnehmer_R3_X = SwenRot;
Data_Teilnehmer_R3_Y = -SwenRot2;
Data_Teilnehmer_R3_Z = SwenRot1;

clear SwenRHandOut
clear SwenRHandOut1
clear SwenRHandOut2
clear SwenRot
clear SwenRot1
clear SwenRot2

Data_Teilnehmer_4_X = YakunRHandOut;
Data_Teilnehmer_4_Y = -YakunRHandOut2;
Data_Teilnehmer_4_Z = YakunRHandOut1;

Data_Teilnehmer_R4_X = YakunRot;
Data_Teilnehmer_R4_Y = -YakunRot2;
Data_Teilnehmer_R4_Z = YakunRot1;

clear YakunRHandOut
clear YakunRHandOut1
clear YakunRHandOut2
clear YakunRot
clear YakunRot1
clear YakunRot2

Data_Teilnehmer_5_X = YanhaoRHandOut;
Data_Teilnehmer_5_Y = -YanhaoRHandOut2;
Data_Teilnehmer_5_Z = YanhaoRHandOut1;

Data_Teilnehmer_R5_X = YanhaoRot;
Data_Teilnehmer_R5_Y = -YanhaoRot2;
Data_Teilnehmer_R5_Z = YanhaoRot1;

clear YanhaoRHandOut
clear YanhaoRHandOut1
clear YanhaoRHandOut2
clear YanhaoRot
clear YanhaoRot1
clear YanhaoRot2


%% Trajectory segmentation
%Mark the starting- and end points
start_position_circ{1} = [-0.22,-0.15;-0.15,-0.02;-0.02,-0.16];
start_position_circ{2} = [-0.23,-0.14;-0.13,-0.03;-0.03,-0.13];
start_position_circ{3} = [-0.23,-0.12;-0.12,-0.03;-0.03,-0.13];
start_position_circ{4} = [-0.25,-0.15;-0.15,-0.02;-0.05,-0.13];
start_position_circ{5} = [-0.25,-0.13;-0.15,-0.02;-0.05,-0.10];

end_position_circ{1} = [-0.58,0.11;-0.44,0.18;-0.25,0.17];
end_position_circ{2} = [-0.59,0.12;-0.44,0.18;-0.25,0.17];
end_position_circ{3} = [-0.58,0.13;-0.44,0.20;-0.22,0.17];
end_position_circ{4} = [-0.60,0.12;-0.45,0.20;-0.22,0.20];
end_position_circ{5} = [-0.59,0.13;-0.44,0.20;-0.22,0.20];

Traj = cell(5,21);
TrajRot = cell(5,21);
anzahl = ones(5,21);

for i=1:1:5
    
    rename_temp_X=['Data_Teilnehmer_', num2str(i), '_X'];
    Data_X_temp = eval(rename_temp_X);
    rename_temp_Y=['Data_Teilnehmer_', num2str(i), '_Y'];
    Data_Y_temp = eval(rename_temp_Y);
    rename_temp_Z=['Data_Teilnehmer_', num2str(i), '_Z'];
    Data_Z_temp = eval(rename_temp_Z);

    rename_temp_X=['Data_Teilnehmer_R', num2str(i), '_X'];
    Data_RotX_temp = eval(rename_temp_X);
    rename_temp_Y=['Data_Teilnehmer_R', num2str(i), '_Y'];
    Data_RotY_temp = eval(rename_temp_Y);
    rename_temp_Z=['Data_Teilnehmer_R', num2str(i), '_Z'];
    Data_RotZ_temp = eval(rename_temp_Z);
    
    index = zeros(1,length(Data_X_temp));
    label = zeros(1,length(Data_X_temp));
    
    for j=2:1:length(Data_X_temp)-1
        
        dis_1 = norm([Data_X_temp(j),Data_Y_temp(j)] - start_position_circ{i}(1,:));
        dis_2 = norm([Data_X_temp(j),Data_Y_temp(j)] - start_position_circ{i}(2,:));
        dis_3 = norm([Data_X_temp(j),Data_Y_temp(j)] - start_position_circ{i}(3,:));

        dis_4 = norm([Data_X_temp(j),Data_Y_temp(j)] - end_position_circ{i}(1,:));
        dis_5 = norm([Data_X_temp(j),Data_Y_temp(j)] - end_position_circ{i}(2,:));
        dis_6 = norm([Data_X_temp(j),Data_Y_temp(j)] - end_position_circ{i}(3,:));
        
        if (Data_Z_temp(j) < Data_Z_temp(j+1) && Data_Z_temp(j) < Data_Z_temp(j-1))
            z_min = 1;
        else 
            z_min = 0;
        end
        
        if (dis_1 < 0.05 && z_min == 1)
            
            index(j) = 1;
            
        elseif (dis_2 < 0.05 && z_min == 1)
            
            index(j) = 2;
        
        elseif (dis_3 < 0.05 && z_min == 1)     
            
            index(j) = 3;
        
        elseif (dis_4 < 0.05 && z_min == 1)     
            
            index(j) = 4;
            
        elseif (dis_5 < 0.05 && z_min == 1)     
            
            index(j) = 5;   
            
        elseif (dis_6 < 0.05 && z_min == 1)     
            
            index(j) = 6;            
        end
        
    end
    
    % Cut trajectories based on starting- and end points (21 categories)
    index_temp = find(index);    %position (on the vector) where trajectory ends
    
    for k = 1:1:length(index_temp)-1
        
        if index(index_temp(k))==1 && index(index_temp(k+1))==4
            
            %Trajectory from data person i, first trajectory.
            
            Traj{i,1}{anzahl(i,1)} = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,1}{anzahl(i,1)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,1) = anzahl(i,1)+1;
            
        elseif index(index_temp(k))==4 && index(index_temp(k+1))==1
            
            Traj{i,2}{anzahl(i,2)}  = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,2}{anzahl(i,2)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,2) = anzahl(i,2)+1; 
            
        elseif index(index_temp(k))==1 && index(index_temp(k+1))==5
            
            Traj{i,3}{anzahl(i,3)}  = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,3}{anzahl(i,3)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,3) = anzahl(i,3)+1; 
                                
        elseif index(index_temp(k))==5 && index(index_temp(k+1))==1
            
            Traj{i,4}{anzahl(i,4)}  = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,4}{anzahl(i,4)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,4) = anzahl(i,4)+1; 
            
        elseif index(index_temp(k))==1 && index(index_temp(k+1))==6
            
            Traj{i,5}{anzahl(i,5)}  = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,5}{anzahl(i,5)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,5) = anzahl(i,5)+1; 
            
        elseif index(index_temp(k))==6 && index(index_temp(k+1))==1
            
            Traj{i,6}{anzahl(i,6)}  = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,6}{anzahl(i,6)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,6) = anzahl(i,6)+1; 
            
        elseif index(index_temp(k))==2 && index(index_temp(k+1))==4
            
            Traj{i,7}{anzahl(i,7)}  = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,7}{anzahl(i,7)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
             
            anzahl(i,7) = anzahl(i,7)+1; 
                       
        elseif index(index_temp(k))==4 && index(index_temp(k+1))==2
            
            Traj{i,8}{anzahl(i,8)}  = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
             TrajRot{i,8}{anzahl(i,8)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
             
            anzahl(i,8) = anzahl(i,8)+1; 
                               
        elseif index(index_temp(k))==2 && index(index_temp(k+1))==5
            
            Traj{i,9}{anzahl(i,9)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,9}{anzahl(i,9)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,9) = anzahl(i,9)+1; 
            
        elseif index(index_temp(k))==5 && index(index_temp(k+1))==2
            
            Traj{i,10}{anzahl(i,10)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,10}{anzahl(i,10)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,10) = anzahl(i,10)+1; 
            
        elseif index(index_temp(k))==2 && index(index_temp(k+1))==6
            
            Traj{i,11}{anzahl(i,11)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,11}{anzahl(i,11)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,11) = anzahl(i,11)+1; 
                        
        elseif index(index_temp(k))==6 && index(index_temp(k+1))==2
            
            Traj{i,12}{anzahl(i,12)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,12}{anzahl(i,12)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,12) = anzahl(i,12)+1; 
                        
        elseif index(index_temp(k))==3 && index(index_temp(k+1))==4
            
            Traj{i,13}{anzahl(i,13)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,13}{anzahl(i,13)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,13) = anzahl(i,13)+1; 
                        
        elseif index(index_temp(k))==4 && index(index_temp(k+1))==3
            
            Traj{i,14}{anzahl(i,14)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,14}{anzahl(i,14)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,14) = anzahl(i,14)+1; 
                                
        elseif index(index_temp(k))==3 && index(index_temp(k+1))==5
            
            Traj{i,15}{anzahl(i,15)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,15}{anzahl(i,15)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,15) = anzahl(i,15)+1; 
            
        elseif index(index_temp(k))==5 && index(index_temp(k+1))==3
            
            Traj{i,16}{anzahl(i,16)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,16}{anzahl(i,16)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,16) = anzahl(i,16)+1; 
            
        elseif index(index_temp(k))==3 && index(index_temp(k+1))==6
            
            Traj{i,17}{anzahl(i,17)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,17}{anzahl(i,17)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,17) = anzahl(i,17)+1; 
                        
        elseif index(index_temp(k))==6 && index(index_temp(k+1))==3
            
            Traj{i,18}{anzahl(i,18)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];            
            TrajRot{i,18}{anzahl(i,18)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,18) = anzahl(i,18)+1; 
            
        elseif index(index_temp(k))==1 && index(index_temp(k+1))==2
            
            Traj{i,19}{anzahl(i,19)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,19}{anzahl(i,19)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,19) = anzahl(i,19)+1; 
            
        elseif index(index_temp(k))==2 && index(index_temp(k+1))==3
            
            Traj{i,20}{anzahl(i,20)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];
            TrajRot{i,20}{anzahl(i,20)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
            
            anzahl(i,20) = anzahl(i,20)+1; 
                        
        elseif index(index_temp(k))==3 && index(index_temp(k+1))==1
            
            Traj{i,21}{anzahl(i,21)}   = [Data_X_temp(index_temp(k):index_temp(k+1)),Data_Y_temp(index_temp(k):index_temp(k+1)),Data_Z_temp(index_temp(k):index_temp(k+1))];    
            TrajRot{i,21}{anzahl(i,21)} = [Data_RotX_temp(index_temp(k):index_temp(k+1)),Data_RotY_temp(index_temp(k):index_temp(k+1)),Data_RotZ_temp(index_temp(k):index_temp(k+1))];
             
            anzahl(i,21) = anzahl(i,21)+1; 
                       
        end
    
    end
    
    
end

%% Testplot
% figure(1)
% index_1 = find(index==1);
% scatter3(Data_X_temp(index_1),Data_Y_temp(index_1),Data_Z_temp(index_1),'b');
% hold on
% index_2 = find(index==2);
% scatter3(Data_X_temp(index_2),Data_Y_temp(index_2),Data_Z_temp(index_2),'r');
% hold on
% index_3 = find(index==3);
% scatter3(Data_X_temp(index_3),Data_Y_temp(index_3),Data_Z_temp(index_3),'y');
% hold on
% index_4 = find(index==4);
% scatter3(Data_X_temp(index_4),Data_Y_temp(index_4),Data_Z_temp(index_4),'m');
% hold on
% index_5 = find(index==5);
% scatter3(Data_X_temp(index_5),Data_Y_temp(index_5),Data_Z_temp(index_5),'g');
% hold on
% index_6 = find(index==6);
% scatter3(Data_X_temp(index_6),Data_Y_temp(index_6),Data_Z_temp(index_6),'k');
% hold on


% m = 1;
% n = 5;

% for i=1:1:length(Traj{m,n}) %m person, doing the n trajectory, attempt i
%     scatter3(TrajRot{m,n}{i}(:,1),TrajRot{m,n}{i}(:,2),TrajRot{m,n}{i}(:,3),'b')
%     hold on
% end

% figure(3)

dt = 0.008333; %time interval

rotationFilter = designfilt('lowpassfir', ...
    'PassbandFrequency',0.04,'StopbandFrequency',0.8, ...
    'PassbandRipple',1,'StopbandAttenuation',20, ...
    'DesignMethod','equiripple');
positionFilter = rotationFilter;


%% position derivatives

d = positionFilter;
for m = 1:5
    
for n = 1:21

for i = 1:1 ...length(Traj{m,n})
Xl = [0; diff(Traj{m,n}{i}(:,1))]/dt; %dx/dt, first velocity starts at 0, also makes sure vectors are the same size
Yl = [0 ;diff(Traj{m,n}{i}(:,2))]/dt;
Zl = [0 ;diff(Traj{m,n}{i}(:,3))]/dt;

Xll = [0; diff(Xl)]/dt;
Yll = [0; diff(Yl)]/dt;
Zll = [0; diff(Zl)]/dt;

% filteredXl = filtfilt(b,a,Xl);
% filteredYl = filtfilt(b,a,Yl);
% filteredZl = filtfilt(b,a,Zl);
% 
% filteredAccX = filtfilt(b,a,Xll);
% filteredAccY = filtfilt(b,a,Yll);
% filteredAccZ = filtfilt(b,a,Zll);
filteredXl = filtfilt(d,Xl);
filteredYl = filtfilt(d,Yl);
filteredZl = filtfilt(d,Zl);

filteredAccX = filtfilt(d,Xll);
filteredAccY = filtfilt(d,Yll);
filteredAccZ = filtfilt(d,Zll);

TrajVelUnfiltered{m,n}{i} = [Xl, Yl, Zl];
TrajAccUnfiltered{m,n}{i} = [Xll, Yll, Zll];
TrajVel{m,n}{i} = [filteredXl, filteredYl, filteredZl];
TrajAcc{m,n}{i} = [filteredAccX, filteredAccY, filteredAccZ];

end
end
end

figure;
hold all
plot(TrajVel{1,1}{1})
plot(TrajVelUnfiltered{1,1}{1})

figure;
hold all
plot(TrajAcc{1,1}{1},'LineWidth',1.5)
plot(TrajAccUnfiltered{1,1}{1})

%% rotation derivatives

d = rotationFilter;

for m = 1:5
    
for n = 1:21

for i = 1:1 ...length(Traj{m,n})
Xl = [0; diff(TrajRot{m,n}{i}(:,1))]/dt; %dx/dt, first velocity starts at 0, also makes sure vectors are the same size
Yl = [0 ;diff(TrajRot{m,n}{i}(:,2))]/dt;
Zl = [0 ;diff(TrajRot{m,n}{i}(:,3))]/dt;

Xll = [0; diff(Xl)]/dt;
Yll = [0; diff(Yl)]/dt;
Zll = [0; diff(Zl)]/dt;

filteredXl = filtfilt(d,Xl);
filteredYl = filtfilt(d,Yl);
filteredZl = filtfilt(d,Zl);

filteredAccX = filtfilt(d,Xll);
filteredAccY = filtfilt(d,Yll);
filteredAccZ = filtfilt(d,Zll);

TrajAngVelUnfiltered{m,n}{i} = [Xl, Yl, Zl];
TrajAngAccUnfiltered{m,n}{i} = [Xll, Yll, Zll];

TrajAngVel{m,n}{i} = [filteredXl, filteredYl, filteredZl];
TrajAngAcc{m,n}{i} = [filteredAccX, filteredAccY, filteredAccZ];

end
end
end

figure;
hold all
plot(TrajAngVel{1,1}{1}(:,1),'LineWidth',1.5)
plot(TrajAngVelUnfiltered{1,1}{1}(:,1))

figure;
hold all
plot(TrajAngAcc{1,1}{1}(:,1),'LineWidth',1.5)
plot(TrajAngAccUnfiltered{1,1}{1}(:,1))

clear a b
m = 1;
a = []; b = a;
for n = 1:21
    temp = [(TrajAcc{m,n}{1} - mean(TrajAcc{m,n}{1}))./max((TrajAcc{m,n}{1} - mean(TrajAcc{m,n}{1})))...
        (TrajVel{m,n}{1} - mean(TrajVel{m,n}{1}))./max((TrajVel{m,n}{1} - mean(TrajVel{m,n}{1}))) ...
        (Traj{m,n}{1} - mean(Traj{m,n}{1}))./max((Traj{m,n}{1} - mean(Traj{m,n}{1})))...
        (TrajAngAcc{m,n}{1} - mean(TrajAngAcc{m,n}{1}))./max((TrajAngAcc{m,n}{1} - mean(TrajAngAcc{m,n}{1}))) ...
        (TrajAngVel{m,n}{1} - mean(TrajAngVel{m,n}{1}))./max((TrajAngVel{m,n}{1} - mean(TrajAngVel{m,n}{1})))...
        (TrajRot{m,n}{1} - mean(TrajRot{m,n}{1}))./max((TrajRot{m,n}{1} - mean(TrajRot{m,n}{1})))];
    a = [a;temp];
    b = [b; n*ones(size(Traj{m,n}{1},1),1)];
end
go  = table(a,b)