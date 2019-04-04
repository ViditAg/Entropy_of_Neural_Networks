%To calculate the network entropy
% function: Entropy 
% input-> Probability distribution of network activity S
% output-> Shannon Entropy
function H=Entropy(Prob_den)
H = - sum(Prob_den(Prob_den>0).*log2(Prob_den(Prob_den>0)));% shannon entropy
