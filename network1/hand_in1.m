%% Hand In No. 1 - Network Dynamics, Gustaf Sundell

%% Loading data 
clear all
close all
clc
% load -ASCII twitter.mat
% load IOdownload.mat
% load -ASCII users.mat
% ^how to load each set of data.
%% 1 - Centrality in input-output network of goods
load IOdownload.mat

n = 47; % number of vertices (sectors). We're dealing with 39 networks, essentially. 

one_vector = ones(n,1);

W_swe = io.swe2000;
W_idn = io.idn2000;

w_out_swe = W_swe*one_vector;
[~, indexes] = sort(w_out_swe, 'descend');
temp = name(indexes);
disp('Out-degree centrality of Sweden')
temp(1:3)
w_out_swe(indexes(1:3))

w_in_swe = W_swe'*one_vector;
[~, indexes] = sort(w_in_swe, 'descend');
temp = name(indexes);
disp('In-degree centrality of Sweden')
temp(1:3)
w_in_swe(indexes(1:3))

w_out_idn = W_idn*one_vector;
[~, indexes] = sort(w_out_idn, 'descend');
temp = name(indexes);
disp('Out-degree centrality of Indonesia')
temp(1:3)
w_out_idn(indexes(1:3))

w_in_idn = W_idn'*one_vector;
[~, indexes] = sort(w_in_idn, 'descend');
temp = name(indexes);
disp('In-degree centrality of Indonesia')
temp(1:3)
w_in_idn(indexes(1:3))
%% 1.b
    
% G_swe = digraph(W_swe, name); %Messy plots when names are included.
% G_idn = digraph(W_idn, name);

G_swe = digraph(W_swe);
G_idn = digraph(W_idn);

plot(G_swe)
title('Sweden, entire network')
figure
plot(G_idn)
title('Indonesia, entire network')

[bin_swe, bin_size_swe] = conncomp(G_swe); % Bin gives a group to each node. bin_size tells how large groups are. 

component_swe = find(bin_swe == find(bin_size_swe == max(bin_size_swe)));
sub_swe = subgraph(G_swe,component_swe);

sub_W_swe = full(adjacency(sub_swe,'weighted'));

[V_swe,D_swe] = eig(sub_W_swe'); 
z_swe = V_swe(:,D_swe == max(eig(sub_W_swe')));

z_swe = z_swe./sum(z_swe);

temp_name = name(component_swe);


disp('Eigen vector centrality of Swedens largest connected component')
[~, indexes] = sort(z_swe, 'descend');
temp = temp_name(indexes);
temp(1:3)
z_swe(indexes(1:3))


[bin_idn, bin_size_idn] = conncomp(G_idn); % for idn, the graph is already connected. 

[V_idn,D_idn] = eig(W_idn'); 
z_idn= V_idn(:,D_idn== max(eig(W_idn'))); %becomes negative? 
%kolla på längst från origo, t.ex. genom att normera med summan. 
z_idn = z_idn./sum(z_idn);
[~, indexes] = sort(z_idn, 'descend');
temp = name(indexes);
disp('Eigen vector centrality of Indonesia')
temp(1:3)
z_idn(indexes(1:3))


%% 1.c Katz 
beta = 0.15;
mu1 = ones(n,1); % OBS lÃ¤gg till andra alternativa mu.

lambda_swe = max(eig(W_swe'));
lambda_idn = max(eig(W_idn'));
katz_swe1 = (eye(n)-((1-beta)/lambda_swe)*W_swe')\mu1*beta;
[~,indexes] = sort(katz_swe1, 'descend'); 
temp = name(indexes);
disp('Katz centrality of Sweden with mu = ones')
temp(1:3)
katz_swe1(indexes(1:3))

katz_idn1 = (eye(n)-((1-beta)/lambda_idn)*W_idn')\mu1*beta;

[~,indexes] = sort(katz_idn1, 'descend'); 
temp = name(indexes);
disp('Katz centrality of Indonesia with mu = ones')
temp(1:3)
katz_idn1(indexes(1:3))

k = 31; % index of Wholesale & retail

mu2 = zeros(n,1);
mu2(k) = 1;

katz_swe2 = (eye(n)-((1-beta)/lambda_swe)*W_swe')\mu2*beta;

[~,indexes] = sort(katz_swe2, 'descend'); 
temp = name(indexes);
disp('Katz centrality of Sweden with mu = 1 only for Wholesale and retail')
temp(1:3)
katz_swe2(indexes(1:3))

katz_idn2 = (eye(n) - ((1-beta)/lambda_idn)*W_idn')\mu2*beta;
[~,indexes] = sort(katz_idn2, 'descend'); 
temp = name(indexes);
disp('Katz centrality of Indonesia with mu = 1 only for Wholesale and retail')
temp(1:3)
katz_idn2(indexes(1:3))
%% 2. - setup
clear all
load -ASCII twitter.mat
load -ASCII users.mat
clc
%% 2.1

W = spconvert(twitter);

s = size(W);
diff = max(s)-min(s); % has more rows than cols, zero pad on cols
W = [W, zeros(max(s),diff)];

n = max(s); % should now be size of the problem, with W as square nxn

beta = 0.15;
mu = sparse(ones(n,1)./n);
% ^gives the same result as when not dividing by n, but slower computation.


w = sum(W,2);


for i = 1:n
   if w(i) == 0
       W(i,i) = 1;
       w(i) = 1;
   end
end
P = sparse(diag(w)\W);

% G_twitter = digraph(W);
% plot(G_twitter);
% [a,b] = conncomp(G_twitter); % Shows that the graph is not connected. 

Ps = sparse(eye(n));
z= sparse(zeros(n,1));
k = 0;
dz = beta*mu;

z = z + dz;
while norm(dz) > 1e-6
   k=k+1;
   Ps = P'*Ps;
   dz = beta*(1-beta)^k*Ps*mu;
   z = z+dz;
end

[~, ind] = sort(z,1,'descend');

top_5_ind = ind(1:5);
disp('5 highest pagerank users')
top_5_users = users(top_5_ind);
for i = 1:5
    top_5_users(i)
end
z(top_5_ind)

%% 2. b

iter = 1000;

stubborn = [3 9];
regular = setdiff(1:n, stubborn);

u = [0 1]';

Q = P(regular, regular);
E = P(regular, stubborn);
x = 0.5.*ones(n, iter);
x(stubborn, 1) = u;

for i = 2:iter
   x(regular,i) = Q*x(regular,i-1) + E*x(stubborn,i-1);
   x(stubborn,i) = x(stubborn,i-1);     
end

n1 = find(round(x(:,end),2) == 0.8); %find index of some nodes tending between 0.5 and 1
n2 = find(round(x(:,end),1) == 0.3);%find index of some nodes tending between 0 and 0.5
plot_nodes = [3; 9; n1; n2;];

close all
plot(x(plot_nodes,:)')
title('Subset of nodes (including stubborn)')


%% 2.c

iter = 1000;

stubborn = [top_5_ind(1) top_5_ind(5)];
regular = setdiff(1:n, stubborn);

u = [0 1]';

Q = P(regular, regular);
E = P(regular, stubborn);
x = 0.5.*ones(n, iter);
x(stubborn, 1) = u;

for i = 2:iter
   x(regular,i) = Q*x(regular,i-1) + E*x(stubborn,i-1);
   x(stubborn,i) = x(stubborn,i-1);     
end

hist(x(:,end),30)
title('Histogram with top 1 and 5 pagerank as stubborn (highest with opinion 0)')


