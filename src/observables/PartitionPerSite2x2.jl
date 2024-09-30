


function PartitionPerSite2x2(C,T,tens,gt,i,j)

# d c d c
# b a b a 
# d c d c  
# b a b a 

N = size(gt)[1];
f(x) = mod(x-1,N) + 1;

c1d = getfield(C[1], gt[f(i-1),f(j-1)]);
t1c = getfield(T[1], gt[f(i),f(j-1)]);
t1d = getfield(T[1], gt[f(i+1),f(j-1)]);
c2c = getfield(C[2], gt[f(i+2),f(j-1)]);

t4b = getfield(T[4], gt[f(i-1),f(j)]); 
a1 = getfield(tens, gt[f(i),f(j)]); 
a2 = getfield(tens, gt[f(i+1),f(j)]); 
t2a = getfield(T[2], gt[f(i+2),f(j)]); 

t4d = getfield(T[4], gt[f(i-1),f(j+1)]); 
a3 = getfield(tens, gt[f(i),f(j+1)]); 
a4 = getfield(tens, gt[f(i+1),f(j+1)]); 
t2c = getfield(T[2], gt[f(i+2),f(j+1)]); 

c4b = getfield(C[4], gt[f(i-1),f(j+2)]); 
t3a = getfield(T[3], gt[f(i),f(j+2)]); 
t3b = getfield(T[3], gt[f(i+1),f(j+2)]); 
c3a = getfield(C[3], gt[f(i+2),f(j+2)]); 


C1 = c1d*t1c*t4b*a1;
C2 = t1d*c2c*a2*t2a;
C4 = t4d*a3*c4b*t3a;
C3 = a4*t2c*t3b*c3a;
    
# d c d c
# b a b a 
# d c d c  
# b a b a 

zz = C1*C2*C3*C4;
cc = c1d*c2c*c3a*c4b;

h = c1d*c4b;
h = h*t1c*t3a;
h = h*t1d*t3b;
h = h*c2c*c3a;

v = c4b*c3a;
v = v*t4d*t2c;
v = v*t4b*t2a;
v = v*c1d*c2c;

kk = zz*cc/(v*h);

return kk

end