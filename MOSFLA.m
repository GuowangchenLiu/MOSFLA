clear all; close all;clc;
tic
%% UPDATES of MOSFLA_S2
% 1. memetic optimization
%-----------------------------------
% Thanks for the authors' idea and I rewrote it seriously.
% Reference: Novel Multiobjective Shuffled Frog Leaping Algorithm
% with Application to Reservoir Flood Control Operation
% DOI: 10.1061/(ASCE)WR.1943-5452.0000027
%% initialize parameters
%%%%%%%%%% the entire population %%%%%%%%%%
m = 5; % number of memplexs
n = 10; % number of solution of each memplex
p = m*n; % solutions in the entire population
incyc = 2; % iterations of each memplex = n/m
outcyc = 400; % generations or shuffled iterations
d_max = 0.4; % maximum allowed change step for each frog's position
s = 30; % number of variables
r = 2; % number of minimum objectives
pmax = 1; % x max-limit
pmin = 0; % x min-limit
%% generate population of frogs
x = zeros(p,s+r+2);
x(:,1:s) = pmin+(pmax-pmin)*rand(p,s);
x(:,s+1:s+r) = evaluate_objective(x(:,1:s));
%% initialize archive set
arc_num = 100; % archive size
archive_set = ones(arc_num,s+r+2)*100; % archive size * number of variables
% iteration need
compare_xw = ones(2,s+r+2)*100; % compare between new xw and old xw
%% iterate
for i = 1:outcyc
    %% sorting population
    x(:,s+r+1) = dominate_sort(x,s,r); % domination sort
    x(:,s+r+2) = crowding_distance(x,s,r); % crowding distance
    x = sortrows(x,[s+r+1,-(s+r+2)]); % sorting
    global_rand = randperm(size(x,1)-2,1)+1;
    xg = x(1,:); % global best value
     % partition
     for j=1:m
         x_local = x(j:m:end,:);
         %% local Memetic Evolution within a memeplex
         for k = 1:incyc % iterations of each memplex
             xb = x_local(1,:); % local best value
             xw = x_local(n,:); % local worst value
             % memetic evolution process by ***local best value***
             step = 2*rand.*(xb-xw); % calculate the step
             step(step>d_max) = d_max; % step limited
             new_xw = xw+step; % generating new xw
             new_xw(new_xw(:,1:s)>pmax) = pmax; % max limit
             new_xw(new_xw(:,1:s)<pmin) = pmin; % min limit
             new_xw(s+1:s+r) = evaluate_objective(new_xw);
             % compare old and new xw
             compare_xw(1,:) = xw;
             compare_xw(2,:) = new_xw;
             compare_xw(:,s+r+1) = dominate_sort(compare_xw,s,r);
             % memetic evolution process by ***global best value***
             if compare_xw(2,s+r+1)~=1
                 step = 2*rand.*(xg-xw);
                 step(step>d_max)=d_max; % step limited
                 new_xw = xw+step;
                 new_xw(new_xw(:,1:s)>pmax) = pmax; % max limit
                 new_xw(new_xw(:,1:s)<pmin) = pmin; % min limit
                 new_xw(s+1:s+r) = evaluate_objective(new_xw);
                 % compare old and new xw
                 compare_xw(2,:) = new_xw;
                 compare_xw(:,s+r+1) = dominate_sort(compare_xw,s,r);
                 % if still can not get a better xw
                 if compare_xw(2,s+r+1)~=1
                     if rand(1)>0.5
                         new_xw(1:s)=0.95*xg(1:s)+0.05*rand.*(x(global_rand+1,1:s)-x(global_rand-1,1:s));
                         new_xw(new_xw(:,1:s)>pmax) = pmax; % max limit
                         new_xw(new_xw(:,1:s)<pmin) = pmin; % min limit
                         new_xw(s+1:s+r) = evaluate_objective(new_xw);
                         x_local(n,:) = new_xw;
                     else
                         new_xw(1:s)=mutate(xg(1:s));
                         new_xw(new_xw(:,1:s)>pmax) = pmax; % max limit
                         new_xw(new_xw(:,1:s)<pmin) = pmin; % min limit
                         new_xw(s+1:s+r) = evaluate_objective(new_xw);
                         x_local(n,:) = new_xw;
                     end
                 else
                     x_local(n,:) = new_xw;
                 end
             else
                 x_local(n,:) = new_xw;
             end
         end
         %% shuffled the memeplexes
         x(j:m:end,:)=x_local;
     end
     %% update the archive set
     archive_set = [archive_set;x];
     archive_set(:,s+r+1) = dominate_sort(archive_set,s,r); % domination sort
     archive_set(:,s+r+2) = crowding_distance(archive_set,s,r);
     archive_set = archive_set(archive_set(:,s+r+1)==1,:);
     archive_set = sortrows(archive_set,-(s+r+2));
     if size(archive_set,1) > arc_num
         archive_set=archive_set(1:arc_num,:);
     end
end
% time and plot
toc;                                                % time
output = archive_set(:,s+1:s+r);
plot(output(:,1),output(:,2),'*b');                 % plot
axis([0,1,0,1]);xlabel('F_1');ylabel('F_2');title('ZDT1')
