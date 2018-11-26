function [euler] = quartenionToEulerZYX(quart)
%argument order [qw qx qy qz] out put [x y z]
%Euler angles ZYX


qw = quart(1);
qx = quart(2);
qy = quart(3);
qz = quart(4);
r11 = 2*(qx*qy + qw*qz);
r12 = qw*qw + qx*qx - qy*qy - qz*qz;
r21 = -2*(qx*qz - qw*qy);
r31 = 2*(qy*qz + qw*qx);
r32 = qw*qw - qx*qx - qy*qy + qz*qz;

x = atan2( r31, r32 )
y = asin ( r21 )
z = atan2( r11, r12)
euler = [x y z];
end

