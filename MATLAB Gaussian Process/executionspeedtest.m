temp = 0;
for i = 1:10000
A = (rand(4,64)); B = (rand(64,4)); 
tic;
A*B;
end
for i = 1:1000
A =(rand(4096,8)); B = (rand(8,4096)); 
tic;
A*B;
temp = temp + toc;
end
temp/10000

