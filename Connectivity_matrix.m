function B=Connectivity_matrix(N,k,randlist)
% Setting up Connectivity Matrix
B=rand(N); % initialize as random matrix
imask=randlist<=0.1*N;       %setting fraction to set connectivity matrix at criticality
B(:,imask)=-1*B(:,imask);   %set outgoing connection from inhibitory neurons to be negative
B(rand(N)>k)=0;             %set mean degree
B=B/max(abs(eig(B)));       %enforce largest eigenvalue = 1