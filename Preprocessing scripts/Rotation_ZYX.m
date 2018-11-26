function R_ZYX = Rotation_ZYX( angle_1, angle_2, angle_3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
R_ZYX = [cos(angle_1)*cos(angle_2), -sin(angle_1)*cos(angle_3) + cos(angle_1)*sin(angle_3)*sin(angle_2), sin(angle_3)*sin(angle_1)+ cos(angle_1)*sin(angle_2)*cos(angle_3);
         cos(angle_2)*sin(angle_1), cos(angle_1)*cos(angle_3) + sin(angle_3)*sin(angle_1)*sin(angle_2), -cos(angle_1)*sin(angle_3) + sin(angle_1)*sin(angle_2)*cos(angle_3)
         -sin(angle_2),              cos(angle_2)*sin(angle_3),                                         cos(angle_2)*cos(angle_3);];

end

