%% Hand-in 3 Gustaf Sundell, gu0147su-s
clear all
close all
clc

%% Part 1 setup
Lambda = [0 2/5 1/5 0 0;
          0 0 3/4 1/4 0;
          1/2 0 0 1/2 0;
          0 0 1/3 0 2/3;
          0 1/3 0 1/3 0];
[n,~] = size(Lambda);      
omega = Lambda*ones(n,1);
D = diag(omega);      
P = D\Lambda; % this will not be used, rather P_bar is the relevant matrix!
omega_star = max(omega);
P_bar = Lambda./omega_star;
for i = 1:n
   P_bar(i,i) = 1- sum(P_bar(i,1:i)) - sum(P_bar(i,i+1:end));  
end
P_bar;
L = D-Lambda;%Laplacian matrix



%% Visualization setup
iter = 10000; 

figure
set(gcf,'color','white')

startnode = 2;
x = zeros(n,iter); % x is one for the current particle state
x(startnode,1) = 1; % starting state is node 1
currentnode = startnode;
cords = [ 0 0 % sets coordinates of nodes, looks like the network in instrucktions already.
          2 1
          4 1
          4 -1
          2 -1];
      
      
ticker = 0;
T = zeros(1,iter);
S = zeros(1,iter);

% Transition probability matrix for the directed ring
% P = [0 1 0 0 0;
%      0 0 1 0 0;
%      0 0 0 1 0;
%      0 0 0 0 1;
%      1 0 0 0 0];

% Plot the graph and mark the node that the particle is in with red
subplot(211)
gplot(P_bar,cords,'-k');
hold on
for i = 1:n
    if x(i, 1) == 1
        scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'r');
    else
        scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'w');
    end

end
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')

%---- Simulate the particle moving around ----%

zzz = 1e-6; % time that the particle waits in the node before moving
         % on to the next one. Here we have just unit-time and no
         % randomness...   
%% Time-plotting
for k = 2:iter
    U = rand(1);
    S(k) = -log(U)/omega_star; % how long do we wait until we change node
    % here we can figure out the times. Store the jumps (unnecessary?)
    ticker = ticker+ S(k);% Tk = Tk-1 + Sk alt Si - log at what time we changed
    T(k) = ticker; % store all the ticks, poisson clock.
%     t = linspace(0,ticker, 10); % ev for plotting. 
    currentnode = randsample(n,1,true,P_bar(currentnode,:)); % goes to next node. 
    x(currentnode,k)=1;
%     pause(zzz) % sleep for some time
    
    % plot the new location of the node
    subplot(211)
    for i = 1:5
        if x(i, k) == 1
            scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'r');
        else
            scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'w');
        end
    end
    subplot(212)
    tvec = [0 1:k];
    plot(tvec(1:end-1),(x(:, 1:k)'*(1:n)'), '-o')
end

%% simulate without plots
iter = 10000; 

startnode = 2;
x = zeros(n,iter); % x is one for the current particle state
x(startnode,1) = 1; % starting state is node 1
currentnode = startnode;

ticker = 0;
T = zeros(1,iter);
S = zeros(1,iter);



for k = 2:iter
    U = rand(1);
    S(k) = -log(U)/omega_star; % how long do we wait until we change node
    % here we can figure out the times. Store the jumps 
    ticker = ticker+ S(k);% Tk = Tk-1 + Sk alt Si - log at what time we changed
    T(k) = ticker; % store all the ticks, poisson clock.
%     t = linspace(0,ticker, 10); % ev for plotting. 
    currentnode = randsample(n,1,true,P_bar(currentnode,:)); % goes to next node. 
    x(currentnode,k)=1;
end


%% Calculate return time node a
findnode = 2; % for first task, look for return time to node 2 (called 'a')
returns = find(x(findnode,:)==1);
nbr = length(returns);
obs = zeros(nbr-1,1);
for i = 2:nbr
    if returns(i) ~= returns(i-1)
        obs(i-1) = T(returns(i))-T(returns(i-1));
    end
end
tmp = obs(obs~=0); % not accounting for zeros, if it stayed in state. 
exp_return_time = mean(tmp)


% theoretical return time: 
[V,D] = eig(L');
v = V(:,1); % I've checked manually that this is the one corresponding 
% to the zero eigenvalue of the transposed laplacian. 
v = v/sum(v);
theoretical_return_period = 1/v(findnode)

%% calculate hitting time from o to d

origins = find(x(1,:)==1);
destinations = find(x(5,:)==1);

nbr = min(length(origins),length(destinations));
obs = zeros(nbr,1);
for i = 1:nbr 
    s = origins(i);
    t = min(destinations(destinations > s));
    if isempty(t)
        break
    end
    obs(i) = T(t) - T(s);

 end

tmp = obs(obs~=0); % in case of left over zeros
exp_hitting_time = mean(tmp)

theoretical_hitting_time = (eye(4)-P_bar(1:4,1:4))\ones(4,1); % for nodes 1-4
theoretical_hitting_time(1) % this is the one we're after, from o. 

%% Part 2.a

close all 
clear all
clc

%% Part 2 simulation without plots
iter = 1000;
n = 10;
x = zeros(n,iter+1); % 0 is red, 1 is green. 

W = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
currentnode = randi(n);% randomize first node to "wake up";

U = zeros(1,iter+1);
U(1) = potential(W,x,1);

for t = 1:iter
    x(:,t+1) = x(:,t);
    p = exp(-eta(t)*W(currentnode,:)*cost_func(0,x(:,t)));
    p = p/(exp(-eta(t)*W(currentnode,:)*cost_func(0,x(:,t))) + exp(-eta(t)*W(currentnode,:)*cost_func(1,x(:,t))));
    u = rand(1);
    if u<=p
        new_color = 0;
    else 
        new_color = 1;
    end
    
    currentnode = randi(n);%node at time t+1
    x(currentnode,t+1) = new_color;
    U(t+1) = potential(W,x,t+1);
    if U(t+1) == 0
        break
    end
end

U(t+1)
x(:,t+1)
% as this code section is without visualization, this is to verify that 
% potential is zero and that the colors don't "collide"
% visualization, see next section!

%% setup visualization for part 2

iter = 10000; 
n = 10;

x = zeros(n,iter+1); % 0 is red, 1 is green. 

W = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
currentnode = randi(n);% randomize first node to "wake up";

U = zeros(1,iter+1);
U(1) = potential(W,x,1);

figure
set(gcf,'color','white')


x = zeros(n,iter); % x is one for the current particle state

cords = [ 0 0 % sets coordinates of nodes, looks like the network in instrucktions already.
          0 1
          0 2
          0 3
          0 4
          0 5
          0 6
          0 7
          0 8
          0 9];
      
P = W;
P(2:end-1,:) = 0.5.*W(2:end-1,:);
n = 10;
% Plot the graph and mark the node that the particle is in with red
subplot(211)
gplot(P,cords,'-k');
hold on
for i = 1:n
    if x(i, 1) == 0
        scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'r');
    else
        scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'g');
    end

end
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')

%---- Simulate the particle moving around ----%

zzz = 0.01;
%% Time-plotting part 2

iter = 10000;
n = 10;
x = zeros(n,iter+1); % 0 is red, 1 is green. 

W = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
currentnode = randi(n);% initial node to wake up;

U = zeros(1,iter+1);
U(1) = potential(W,x,1);


for t = 2:iter
    x(:,t+1) = x(:,t);
    p = exp(-eta(t)*W(currentnode,:)*cost_func(0,x(:,t)));
    p = p/(exp(-eta(t)*W(currentnode,:)*cost_func(0,x(:,t))) + exp(-eta(t)*W(currentnode,:)*cost_func(1,x(:,t))));
    u = rand(1);
    if u<=p
        new_color = 0;
    else 
        new_color = 1;
    end
    
    currentnode = randi(n);%node at time t+1
    x(currentnode,t+1) = new_color;
    U(t+1) = potential(W,x,t+1);
    
    
    if U(t+1) == 0
        break
    end
end

subplot(211)
for i = 1:n
    if x(i, t+1) == 0
        scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'r');
    else
        scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'g');
    end
end
subplot(212)
tvec = [0 1:t];
plot(tvec(1:end),U(1:t+1), '-o')

%% Part 2 b
close all
clear all
clc
%% setup
load -ascii wifi.mat
load -ascii coord.mat

% plot(digraph(wifi)) % for curiosity

colors = ['r', 'g', 'b', 'y', 'm', 'c','w', 'k']; % storing the colors for plotting

[n,~] = size(wifi);

w = wifi*ones(n,1);
put_selfloops = find(w==0);
W = wifi;
for i = 1:length(put_selfloops)
   W(put_selfloops(i),put_selfloops(i)) = 1; 
end
w(w==0) = 1;
P = diag(w)\W;

iter = 10000;
figure
set(gcf,'color','white')

x = ones(n,iter+1);
cords = coord;

subplot(211)
gplot(P,cords,'-k');
hold on
for i = 1:n
    scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', colors(x(i,1)));
end
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')



% zzz =0.000001; 

         
%% Simulation and plot

iter = 10000;

currentnode = randi(n);% randomize initial node to wake up;

U = zeros(1,iter+1);
U(1) = potential(W,x,1);

for t = 1:iter
    x(:,t+1) = x(:,t);
    currentnode = randi(n);
    p = P_vector(W,x,t,colors, currentnode);
    new_color = randsample(length(colors),1,true,p);
    x(currentnode,t+1) = new_color;
    U(t+1) = potential(W,x,t+1);

end

subplot(211)
for i = 1:n
   scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', colors(x(i,t))); 
end
subplot(212)
tvec = [0 1:t];
plot(tvec,(U), '-o')


function P = P_vector(W,x,t,colors,i)
    % divisor is sum of elements. 
    P = zeros(length(colors),1);
    for c = 1:length(colors)
       P(c) = exp(-eta(t)*W(i,:)*cost_2(c,x(:,t)));  
    end
    P = P./sum(P);
end



function cost = cost_2(s,xj)
    l = length(xj);
    cost = zeros(l,1);
    for i = 1:l
        if xj(i) == s
            cost(i) = 2;
        elseif abs(xj(i)-s) == 1
            cost(i) = 1;
        else 
            cost(i)= 0;
        end
    end   
end

function inv_noise = eta(t)
    inv_noise = 0.01;
end

function cost = cost_func(s,xj)
    cost = xj == s;
end


function pot = potential(W,x,t)
    sum = 0;
    n = length(x(:,1));
    for i =1:n
        for j=1:n
            sum = sum + W(i,j)*cost_func(x(i,t),x(j,t));
        end
    end
    pot = sum/2;
end


