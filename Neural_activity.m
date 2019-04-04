% function name-> Neural_activity
% Inputs-> N: size of network; T: Time duration of simulation
% trans: transient time steps; randlist: list of random numbers 1 to N 
% alpha: fraction of inhibitory neurons, I,E: inhibitory and excitatory weights
% B_in: Connectivity Matrix; sponr: rate of spontaneous activation
%Output: S, network activity time-series
function S=Neural_activity(N,T,trans,alpha,randlist,I,E,B_in,sponr)
imask1=randlist<=alpha*N; % logical indexing columns for inhibitory neurons
B_in(:,imask1)=-1*B_in(:,imask1);%set outgoing connection from inhibitory neurons to be negative
B_in(:,imask1)=I*B_in(:,imask1);
B_in(:,~imask1)=E*B_in(:,~imask1);
nev=false(N,T);  %initialize matrix for storing activity
%%%%%%%% compute the activity of the network %%%%%%%%%%
%initial condition: activate Ni neurons in first timestep
nev(1:N/2,1)=1;
%evolve dynamics: probabilistic spike propagation
t=1;
while t<T  %stop computing if we reach T steps
    %determine which neurons fire in the next time step
    nev(:, t+1) = B_in*nev(:,t)>rand(N,1);
    %random activation at rate of one spike among all neurons every 100 timesteps
    nev(rand(N,1)<sponr,t+1)=1;
    t=t+1;
end
sumwtrans =sum(nev,1);% complete spike count series
S=sumwtrans(trans+1:end);% removing 1000 time steps of transience data

