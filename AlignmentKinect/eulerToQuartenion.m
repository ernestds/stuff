function [quartenions] = eulerToQuartenion(euler)
% Output order [qw qx qy qz]
%input order [X Y Z]
m = Rotation_ZYX(euler(3),euler(2),euler(1));
qw= sqrt(1 + m(1,1) + m(2,2) + m(3,3)) /2;
qx = (m(3,2) - m(2,3))/( 4 *qw);
qy = (m(1,3) - m(3,1))/( 4 *qw);
qz = (m(2,1) - m(1,2))/( 4 *qw);
quartenions = [qw qx qy qz];

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
end

