%% Hand-in 4 Gustaf Sundell, gu0147su-s
clear all 
close all
clc


%% Part 1 
close all

% Setting up graph's adjacency matrix

n = 500;
% k = 4;
W = zeros(n);
W = W + diag(ones(n-1,1),1); % add ones on the +1 off-diagonal
W = W + diag(ones(n-1,1),-1); % add ones on the -1 off-diagonal
W = W + diag(ones(n-2,1),2); % add ones on the +2 off-diagonal
W = W + diag(ones(n-2,1),-2); % add ones on the -2 off-diagonal
W = W + diag(ones(1,1),n-1); % add ones on the +n-1 off-diagonal
W = W + diag(ones(1,1),1-n); % add ones on the -n+1 off-diagonal
W = W + diag(ones(2,1),n-2); % add ones on the +n-2 off-diagonal
W = W + diag(ones(2,1),2-n); % add ones on the -n+2 off-diagonal
W = sparse(W); % transform it into a sparse matrix
G = graph(W); % convert W into a graph, note not directed, then would be digraph

plot(G) % for curiosity

% Setup simulation parameters. 
rho = 0.7;
beta = 0.3;
N = 100;
n_weeks = 15;
n_infected = 10;
% Run simulation 
[S_mean, I_mean, I_new_mean, R_mean]= SIR_simulation(W,beta,rho,N,n_weeks,n_infected,false);
% Time vector for plotting
timeline = 1:n_weeks;
% Requested plots
figure 
plot(timeline, S_mean)
hold on 
plot(timeline, I_mean)
hold on 
plot(timeline, R_mean)
title('Mean # Susceptible, Infected and Recovered per week, problem 1.1')
legend('S','I','R')
figure 
plot(timeline, I_new_mean)
title('# of New infected per week, problem 1.1')
disp(strcat('#S(15) = ', string(S_mean(n_weeks))))
disp(strcat('#I(15) = ', string(I_mean(n_weeks))))
disp(strcat('#I_new(15) = ', string(I_new_mean(n_weeks))))
disp(strcat('#R(15) = ', string(R_mean(n_weeks))))

%% Problem 1.2 
% Generate graph by preferential attachment algorithm. 
W = gen_population(1000,8,true);
G = graph(W);
figure
plot(G)
% Verify average degree
disp(strcat('Average degree:', string(mean(W*ones(1000,1)))));

w = W*ones(1000,1);
figure
histogram(w);

%% Part 2. 
close all
clc
n = 500;
k = 6;
beta = 0.3;
rho = 0.7;
n_weeks = 15;
N = 100;
n_infected = 10;

W = gen_population(n,k,true);
[S_mean, I_mean, I_new_mean, R_mean] = SIR_simulation(W,beta,rho,N,n_weeks,n_infected,false);

timeline = 1:n_weeks;

% After combining the two algorithms, plotting stats from simulation
figure
plot(timeline, I_new_mean)
title('New infections per week');
figure 
plot(timeline, S_mean)
hold on 
plot(timeline, I_mean)
hold on 
plot(timeline,R_mean)
legend('S', 'I', 'R')
title('Total # of S, I and R over 15 weeks')

%% Part 3 - introducing vaccinations. 

vacc_scheme = [0, 5, 15, 25, 35, 45, 55, 60,60,60,60,60,60,60,60];

vacc_scheme = vacc_scheme./100; 
% normalize vaccination scheme vector so that they're in percent. 

close all
clc
n = 500;
k = 6;
beta = 0.3;
rho = 0.7;
n_weeks = 15;
N = 100;
n_infected = 10;
plt = false;
timeline = 1:n_weeks;


W = gen_population(n,k,false);

[S_mean, I_mean, I_new_mean, R_mean, V_mean, V_new_mean] = SIRV_simulation(W,beta,rho,vacc_scheme,N,n_weeks,n_infected,plt);


figure
plot(timeline, I_new_mean)
hold on     
plot(timeline, V_new_mean)
legend('new infected', 'new vaccinated')
title('New infections and vaccinations per week');
figure 
plot(timeline, S_mean)
hold on 
plot(timeline, I_mean)
hold on 
plot(timeline,R_mean)
hold on 
plot(timeline, V_mean)
legend('S', 'I', 'R', 'V')
title('Total # of S, I, R and V over 15 weeks')


%% Part 4 - optimizing parameters to fit Sweden's data
close all
clc
I_goal = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0];

vacc_scheme = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60, 60];
vacc_scheme = vacc_scheme./100;
n = 934;
N = 10;
n_weeks = 15;
n_infected = 1; 

k0 = 10;
beta0 = 0.3;
rho0 = 0.7;

dk = 2;
dbeta = 0.2;
drho = 0.2;

params = [k0;beta0;rho0]; % keeping current params in vector, with this ordering. 
iter = 0;
halves = 0;
while true
    old_params = params; % to compare if we stayed
    try_these_params = [repmat([floor(params(1)-dk) params(1) ceil(params(1)+dk)],1,9); 
        repmat([repmat(params(2)-dbeta,1,3) repmat(params(2),1,3) repmat(params(2)+dbeta,1,3)],1,3); 
        [(params(3)-drho)*ones(1,9) params(3)*ones(1,9) (params(3)+drho)*ones(1,9)]];
    save_RMSE = zeros(1,27);
    for i = 1:27
        W = gen_population(n,try_these_params(1,i),false); % Generate new pop, as k may vary. 
        I_ref = SIRV_simulation_streamlined(W,try_these_params(2,i),try_these_params(3,i),vacc_scheme,N,n_weeks,n_infected);
        
        to_square = I_ref-I_goal(2:end);% Excluding initial infected.
        save_RMSE(i) = sqrt((to_square*to_square')./n_weeks);
    end
    best = find(save_RMSE==min(save_RMSE));
    if length(best) >1
       best = randsample(best,1); 
    end
    params = try_these_params(:,best);
    if sum(params==old_params)==3 % i.e. no parameter has changed.
        dk = dk/2;
        dbeta = dbeta/2;
        drho = drho/2;
        halves = halves + 1;
    end
    if halves == 2
        break
    end
    iter = iter +1; 
end
% Some disorderly output. 
params
old_params
I_ref
I_goal
save_RMSE(best)
iter

% timeline = 1:n_weeks;
% plot(timeline, I_goal(2:end))
% hold on 
% plot(timeline, I_ref)
% title('Given and simulated newly infected per week')
% legend('Given', 'Simulated')


%% Plots for Problem 4 
k = params(1);
beta = params(2);
rho = params(3);
N = 100;
W = gen_population(n,k,false);

[S_mean, I_mean, I_new_mean, R_mean, V_mean, V_new_mean] = SIRV_simulation(W,beta,rho,vacc_scheme,N,n_weeks,n_infected,false);
timeline = 1:n_weeks;

figure
plot(timeline, I_new_mean)
hold on     
plot(timeline, V_new_mean)
hold on 
plot(timeline, I_goal(2:end))
legend('new infected', 'new vaccinated', 'I_0')
title('New infections, vaccinations and I_0 per week');
figure 
plot(timeline, S_mean)
hold on 
plot(timeline, I_mean)
hold on 
plot(timeline,R_mean)
hold on 
plot(timeline, V_mean)
legend('S', 'I', 'R', 'V')
title('Total # of S, I, R and V over 15 weeks')





%% Functions


function W = gen_population(n,k, plt)
    
    if nargin < 1
        n = 500;
    end
    if nargin<2
        k=6;
    end
    if nargin <3
        plt = true;
    end
    
    k0 = k+1;
    W = ones(k0,k0) - diag(ones(k0,1));
    W = sparse(W);


    for i = (k0+1):n

       if mod(i,2) == 0
           c = floor(k/2);
       else
           c = ceil(k/2);
       end
       w = sum(W,2);
       P=w./sum(w); % probabilities of attachment. 
       for j = 1:c
          neighbour = randsample(1:k0,1,true,full(P));
          P(neighbour)=0;
          W(k0+1,neighbour) = 1;
          W(neighbour,k0+1) = 1;       
       end
           k0 = k0+1;
    end

    if plt
        G = graph(W);
        plot(G);
    end

end



function [S_mean, I_mean, I_new_mean, R_mean]= SIR_simulation(W,beta,rho,N,n_weeks,n_infected,plt)
    
    if nargin <7
        plt = false;
    end
    if nargin <6
        n_infected = 10;
    end
    if nargin<5 
        n_weeks = 15;
    end
    if nargin <4
        N = 100;
    end
    [n,~] = size(W);
    S = 0;
    I = 1;
    R = 2;
    
    S_save = zeros(N,n_weeks);
    I_save = zeros(N,n_weeks);
    R_save = zeros(N,n_weeks);
    I_new_save = zeros(N,n_weeks);

    for k = 1:N
        X = zeros(n,1);
        initial_infected = randsample(1:n,n_infected);
        X(initial_infected) = I;
        for t = 1:n_weeks 
           m = W*(X==I);
           prob_I = ones(n,1)-(1-beta).^m;
           susceptibles = find(X==S);
           infected = find(X==I);
           u = rand([n,1]);
           new_infecteds = (u(susceptibles) < prob_I(susceptibles));
           X(susceptibles) = new_infecteds;
           X(infected) = (u(infected) < rho) + 1;

           S_save(k,t) = sum(X==S);
           I_save(k,t) = sum(X==I);
           I_new_save(k,t) = sum(new_infecteds);
           R_save(k,t) = sum(X==R);

        end   
    end

    timeline = 1:n_weeks;
    
    if plt
       figure
       subplot(221)
       plot(timeline,mean(S_save));
       hold on 
       subplot(222)
       plot(timeline,mean(I_save));
       hold on 
       subplot(223)
       plot(timeline,mean(I_new_save));
       hold on
       subplot(224)
       plot(timeline,mean(R_save));
    end
    S_mean = mean(S_save);
    I_mean = mean(I_save);
    I_new_mean = mean(I_new_save);
    R_mean = mean(R_save);

end

function [S_mean, I_mean, I_new_mean, R_mean, V_mean, V_new_mean] = SIRV_simulation(W,beta,rho,vacc_scheme,N,n_weeks,n_infected,plt)
    
    if nargin <8
        plt = false;
    end
    if nargin <7
        n_infected = 10;
    end
    if nargin<6 
        n_weeks = 15;
    end
    if nargin <5
        N = 100;
    end
    [n,~] = size(W);
    
    S = 0;
    I = 1;
    R = 2;
    V = 3;
    
    S_save = zeros(N,n_weeks);
    I_save = zeros(N,n_weeks);
    R_save = zeros(N,n_weeks);
    I_new_save = zeros(N,n_weeks);
    V_save = zeros(N, n_weeks);
    V_new_save = zeros(N, n_weeks); % Necessary? deterministic?

    for k = 1:N
        X = zeros(n,1);
        initial_infected = randsample(1:n,n_infected);
        X(initial_infected) = I;
        for t = 1:n_weeks 
           % Här vaccinerar vi!
           if t == 1
               nbr_vacc = vacc_scheme(1)*n;
           else 
               nbr_vacc = n*(vacc_scheme(t)-vacc_scheme(t-1));
           end
           vacc_population = find(X ~=V);
           new_vacc = randsample(vacc_population,ceil(nbr_vacc)); %%%%%%%%% should work
           X(new_vacc) = V;
           
           m = W*(X==I);
           prob_I = ones(n,1)-(1-beta).^m;
           susceptibles = find(X==S);
           infected = find(X==I);
           u = rand([n,1]);
           new_infecteds = (u(susceptibles) < prob_I(susceptibles));
           X(susceptibles) = new_infecteds;
           X(infected) = (u(infected) < rho) + 1;

           S_save(k,t) = sum(X==S);
           I_save(k,t) = sum(X==I);
           I_new_save(k,t) = sum(new_infecteds);
           R_save(k,t) = sum(X==R);
           V_save(k,t) = sum(X==V);
           V_new_save(k,t) = length(new_vacc);

        end   
    end

    timeline = 1:n_weeks;
    
    if plt
       figure
       subplot(221)
       plot(timeline,mean(S_save));
       hold on 
       subplot(222)
       plot(timeline,mean(I_save));
       hold on 
       subplot(223)
       plot(timeline,mean(I_new_save));
       hold on
       subplot(224)
       plot(timeline,mean(R_save));
    end
    S_mean = mean(S_save);
    I_mean = mean(I_save);
    I_new_mean = mean(I_new_save);
    R_mean = mean(R_save);
    V_mean = mean(V_save);
    V_new_mean = mean(V_new_save);

end


function I_new_mean = SIRV_simulation_streamlined(W,beta,rho,vacc_scheme,N,n_weeks,n_infected)
    
    [n,~] = size(W);
    
    S = 0;
    I = 1;
    R = 2;
    V = 3;
    
    I_new_save = zeros(N,n_weeks);
    for k = 1:N
        X = zeros(n,1);
        initial_infected = randsample(1:n,n_infected);
        X(initial_infected) = I;
        for t = 1:n_weeks 
           % Här vaccinerar vi!
           if t == 1
               nbr_vacc = vacc_scheme(1)*n;
           else 
               nbr_vacc = n*(vacc_scheme(t)-vacc_scheme(t-1));
           end
           vacc_population = find(X ~=V);
           new_vacc = randsample(vacc_population,ceil(nbr_vacc)); %%%%%%%%% should work
           X(new_vacc) = V;
           
           m = W*(X==I);
           prob_I = ones(n,1)-(1-beta).^m;
           susceptibles = find(X==S);
           infected = find(X==I);
           u = rand([n,1]);
           new_infecteds = (u(susceptibles) < prob_I(susceptibles));
           X(susceptibles) = new_infecteds;
           X(infected) = (u(infected) < rho) + 1;
           I_new_save(k,t) = sum(new_infecteds);
        end   
    end
    I_new_mean = mean(I_new_save);
end