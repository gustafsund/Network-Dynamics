%% Hand-in 2, Gustaf Sundell, gu0147su-s
clear all
load -ascii traveltime.mat % Time in le, taken from length/60mph
load -ascii traffic.mat % B
load -ascii flow.mat 
load -ascii capacities.mat

B = traffic;
l = traveltime;
f = flow;
c = capacities;
clear traveltime traffic flow capacities


%% 1.a)
% get W from the graph in order to create graph object.

[n,~] = size(B); % # of nodes. 
links = length(l); % # of links. 
W = zeros(n,n); % W to be filled. 

o = 1;
d = 17;

for col = 1:links
    start_node = find(B(:,col)==1);
    end_node = find(B(:,col)==-1);
    W(start_node, end_node) = l(col);    
end
A = sparse(W);

[dist,path,~] = graphshortestpath(A,o,d)


%% 1.b)
% we need to get the capacities in sparse nxn format, analogous to the
% adjacency matrix.
C_square = zeros(links,links);
for col = 1:links
    start_node = find(B(:,col)==1);
    end_node = find(B(:,col)==-1);
    C_square(start_node, end_node) = c(col);    
end
C_sparse = sparse(C_square);

max_flow = graphmaxflow(C_sparse, o, d)


%% 1.c) 
nu = B*f


%% 1.d) 

net_inflow = nu(o);
lambda = zeros(n,1);
lambda(o) = net_inflow;
mu = zeros(n,1);
mu(d) = net_inflow;

f_star = f;
cvx_begin
    variable f_star(links)
    minimize sum(l.*c.*inv_pos(1-f_star./c)-l.*c)
    subject to
        B*f_star == lambda-mu
        0 <= f_star <= c
cvx_end


%% 1.e)

% The integral is done analytically on paper and the evaluation is summed in 
% the goal function
f_wardrop = f;
cvx_begin
    variable f_wardrop(links)
    minimize sum(-l.*c.*log(1-f_wardrop.*inv_pos(c)))
    subject to 
        B*f_wardrop == lambda-mu
        0<=f_wardrop<=c
cvx_end


%% 1.f)
omega = f_star.*(l./c).*(inv_pos(1-f_star./c).^2);

f_tolls = f;
cvx_begin
    variable f_tolls(links)
    minimize sum(-l.*c.*log(1-f_tolls.*inv_pos(c))+f_tolls.*omega)
    subject to 
        B*f_tolls== lambda-mu
        0<=f_tolls<=c
cvx_end


%% 1.g) - getting f_g
f_g = f;
cvx_begin
    variable f_g(links)
    minimize sum(l.*c.*inv_pos(1-f_g./c)-l.*c-l.*f_g)
    subject to 
        B*f_g== lambda-mu
        0<=f_g<=c
cvx_end

% now for the optimal tolls. 
omega_opt =f_g.*c.*l.*inv_pos(pow_p(c-f_g,2)) - l;
f_gtest = f;

cvx_begin
    variable f_gtest(links)
    minimize sum(-l.*c.*log(1-f_gtest.*inv_pos(c))+f_gtest.*omega_opt)
    subject to 
        B*f_gtest== lambda-mu
        0<=f_gtest<=c
cvx_end


