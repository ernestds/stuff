display(' ')
display("Angulos de Euler no sistema original")
origz = -0.616871928445094;
origy = 3.870267616297426e-04;
origx = -0.008756037604575;
m =[   0.815691930336862   0.578461512163490   0.005380866923830;...
  -0.578486408650728   0.815662683173502   0.006918257337391;...
  -0.000387026751968  -0.008755925064298   0.999961591255665]
display(' ')
display("Vetores unitarios apos rotacao")
    x_vec = m*[1;0;0];
    y_vec = m*[0;1;0];
    z_vec = m*[0;0;1];
    display(' ')
display("Quartenions a partir da matriz de rotacao criada pelos angulos de Euler originais")
qw= sqrt(1 + m(1,1) + m(2,2) + m(3,3)) /2
qx = (m(3,2) - m(2,3))/( 4 *qw)
qy = (m(1,3) - m(3,1))/( 4 *qw)
qz = (m(2,1) - m(1,2))/( 4 *qw)

%the original data follows ZYX; consider though that the coordinate used
%for position is different, as in new = [1 0 0; 0 0 -1; 0 1 0]*old
display(' ')
display("Angulos de Euler a partir dos quartenions")
r11 = 2*(qx*qy + qw*qz);
r12 = qw*qw + qx*qx - qy*qy - qz*qz;
r21 = -2*(qx*qz - qw*qy);
r31 = 2*(qy*qz + qw*qx);
r32 = qw*qw - qx*qx - qy*qy + qz*qz;

x1 = atan2( r31, r32 )
y1 = asin ( r21 )
z1 = atan2( r11, r12)

% w = c1 c2 c3 - s1 s2 s3
% x = s1 s2 c3 +c1 c2 s3
% y = s1 c2 c3 + c1 s2 s3
% z = c1 s2 c3 - s1 c2 s3

%     c1 = cos(heading / 2) heading y
%     c2 = cos(attitude / 2) attitude z
%     c3 = cos(bank / 2) bank x
%     s1 = sin(heading / 2)
%     s2 = sin(attitude / 2)
%     s3 = sin(bank / 2)



display(' ')
display("Transformando para o novo sistema de coordenadas")
R = [1,0,0;0,0,-1;0,1,0];
m = R*m;
display(' ')
display("Vetores unitarios apos rotacao")
    x_vec = m*[1;0;0]
    y_vec = m*[0;1;0]
    z_vec = m*[0;0;1]
        display(' ')
display("Quartenions a partir da matriz de rotacao criada pelos angulos de Euler originais")
qw= sqrt(1 + m(1,1) + m(2,2) + m(3,3)) /2
qx = (m(3,2) - m(2,3))/( 4 *qw)
qy = (m(1,3) - m(3,1))/( 4 *qw)
qz = (m(2,1) - m(1,2))/( 4 *qw)

%the original data follows ZYX; consider though that the coordinate used
%for position is different, as in new = [1 0 0; 0 0 -1; 0 1 0]*old

r11 = 2*(qx*qy + qw*qz);
r12 = qw*qw + qx*qx - qy*qy - qz*qz;
r21 = -2*(qx*qz - qw*qy);
r31 = 2*(qy*qz + qw*qx);
r32 = qw*qw - qx*qx - qy*qy + qz*qz;

x = atan2( r31, r32 )
y = asin ( r21 )
z = atan2( r11, r12)