function m  = massmatrix(rho,area,l,Ip)
%MASSMATRIX mass matrix of each element in local coordinate system
% 
% 
m=rho*area*l*[1/3        0          0          0         0          0    1/6         0          0           0          0          0
               0     13/35          0          0         0   11*l/210      0      9/70          0           0          0  -13*l/420
               0         0      13/35          0 -11*l/210          0      0         0       9/70           0   13*l/420          0
               0         0          0  Ip/3/area         0          0      0         0          0   Ip/6/area          0          0
               0         0  -11*l/210          0   l^2/105          0      0         0  -13*l/420           0   -l^3/140          0
               0  11*l/210          0          0         0    l^2/105      0  13*l/420          0           0          0   -l^3/140
             1/6         0          0          0         0          0    1/3         0          0           0          0          0
               0      9/70          0          0         0   13*l/420      0     13/35          0           0          0  -11*l/210
               0         0       9/70          0 -13*l/420          0      0         0      13/35           0   11*l/210          0
               0         0          0  Ip/6/area         0          0      0         0          0   Ip/3/area          0          0
               0         0   13*l/420          0  -l^3/140          0      0         0   11*l/210           0    l^2/105          0
               0 -13*l/420          0          0         0   -l^3/140      0 -11*l/210          0           0          0    l^2/105];
end

