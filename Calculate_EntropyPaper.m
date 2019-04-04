task=1;% Set 1 for figure 1 data; 2 for figure 2 and 3 for figure 3
rng('shuffle'); % randomize initial seed
N=10000;  %total number of neurons
T=21000; %duration of simulation
trans = 1000; % transient steps
%spontaneous activation rate
sponr=1/N/100; %one spike per 100 time steps among all neurons
k=0.01;%mean degree N*k
bins = linspace(0,1/N,1);% defining bin size for histogram of neural activity
[~,randlist] = sort(rand(1,N));% initiate random number list, preferably keep it same
B=Connectivity_matrix(N,k,randlist); % Connectivity Matrix
B=abs(B); %make B non-negative to run inside alpha loop
if task==1
    alphalist= 0.09:0.01:0.11;% list of alpha(fraction of inhibitory neurons)
    Ifac= 1.25;% W_E(mean excitatory weight)
    Efac= 1.25;%W_I(mean inhibitory weight)
    for a=1:length(alphalist)
        % Neural acitivity, S
        S(:,a)=Neural_activity(N,T,trans,alphalist(a),randlist,Ifac,Efac,B,sponr)/N;
        h = histc(S(:,a),bins); % make histogram to get S distribution
        % Probability P(S)
        Prob(:,a) = h/(T-trans);
    end
elseif task==2
    % list of alpha(fraction of inhibitory neurons)
    alphalist= [0.1,0.2]; % 0.01:0.01:0.36;
    % W_E(mean excitatory weight)
    Ifac= 1.25:0.25:3.25; % [1.5,2.5];
    %W_I(mean inhibitory weight)
    Efac= 1.25:0.25:3.25; % [1.5,2.5];
    for a=1:length(alphalist)
        for i=1:length(Ifac)
            for e=1:length(Efac)
                % Neural acitivity, S
                S=Neural_activity(N,T,trans,alphalist(a),randlist,Ifac(i),Efac(e),B,sponr)/N;
                h = histc(S,bins); % make histogram to get S distribution
                Prob = h/(T-trans);% Probability P(S)
                H(a,i,e) = Entropy(Prob); % Shannon entropy
            end
        end
    end
elseif task==3
    alphalist= 0.01:0.01:0.7;
    % W_E(mean excitatory weight)
    Ifac= 1.25:0.25:3.25;
    %W_I(mean inhibitory weight)
    Efac= 1.25:0.25:3.25;
    for i=1:length(Ifac)
        for e=1:length(Efac)
            for a=1:length(alphalist)
                % Neural acitivity, S
                S=Neural_activity(N,T,trans,alphalist(a),randlist,Ifac(i),Efac(e),B,sponr)/N;
                h = histc(S,bins); % make histogram to get S distribution
                Prob = h/(T-trans);% Probability P(S)
                H(a) = Entropy(Prob); % Shannon entropy
            end
            [MaxH(e,i),Max_idx] = max(H); % maximum entropy at fixed W_E and W_I
            Alpha_crit(e,i) = alphalist(Max_idx); % crtical alpha at which entropy maximizes
        end
    end
    [Igrid,Egrid] = meshgrid(Efac,Ifac); % make gridpoints
    [nx,ny,nz] = surfnorm(Egrid,Igrid,A_crit);% normal vector at each grid point
    toc
    %points normal*const,c distance away
    c=0.01;
    % above normal
    E1 = Egrid + nx*c; % excitatory weight
    I1 = Igrid + ny*c; % inhibitory weight
    A1 = A_crit + nz*c; % fraction of inhibitory neurons
    % below normal
    E2 = Egrid - nx*c; % excitatory weight
    I2 = Igrid - ny*c; % inhibitory weight
    A2 = A_crit - nz*c;% fraction of inhibitory neurons
    % entropy values at both normal distances
    for e=1:length(Efac)
        for i=1:length(Ifac)
            S=Neural_activity(N,T,trans,A1(e,i),randlist,I1(e,i),E1(e,i),B,sponr)/N;
            h = histc(S,bins); % make histogram to get S distribution
            Prob = h/(T-trans);% Probability P(S)
            H1(e,i) = Entropy(Prob); % Shannon entropy
            S=Neural_activity(N,T,trans,A2(e,i),randlist,I2(e,i),E2(e,i),B,sponr)/N;
            h = histc(S,bins); % make histogram to get S distribution
            Prob = h/(T-trans);% Probability P(S)
            H2(e,i) = Entropy(Prob); % Shannon entropy
        end
    end
    % fragility calculated at all W_E and W_I values
    for e=1:length(Efac)
        for i=1:length(Ifac)
            Fragility(e,i)=((MaxH(e,i)-H1(e,i))+(MaxH(e,i)-H2(e,i)))/2;
        end
    end
else
    disp('Enter 1, 2 or 3 only')
end