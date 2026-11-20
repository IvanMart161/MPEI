function [expPfa] = monteCarloFalseAlarm(T,noisePwr)
% This function performs a Monte Carlo Importance Sampling to determine
% experimental probabilities of false alarm
%
% INPUTS
% T - This is the threshold
%
% noisePwr - the noise power in linear scales
%
% OUTPUTS
% expPfa: The experimental probability of false alram
%
% ASSUMPTIONS
% Exponential distribution of noise/clutter
%
% Kang K.E. Radar System: Analysis, Design and Simulation. Chapter 9:
% Monte Carlo Method and Function Integration. pp. 247-260

nCFAR = size(T,2);
nRangeGate = length(T); % The number of range gates
monteCarloPwr = noisePwr*20; % Monte Carlo mean 20 times greater
iCFARwin = 1:nRangeGate;
iCFARwin = iCFARwin(T(:,1)~=0); % indicies for non-zero thresholds
nTrial = length(iCFARwin);
% Exponential distribution with the Monte carlo power
pwrMonteCarlo = exprnd(monteCarloPwr,nRangeGate,nCFAR);
% The Monte carlo weight
weight = monteCarloPwr/noisePwr*exp(-(1/noisePwr - 1/monteCarloPwr)*pwrMonteCarlo);
% with iCFAR win, we avoid the zero thresholsd at the ends of the range gate
weight = weight(iCFARwin,:);
T = T(iCFARwin,:);
pwrMonteCarlo = pwrMonteCarlo(iCFARwin,:);
% Logical indicies for when Monte Carlo power greater than threshold
iFalseAlarm = T < pwrMonteCarlo;
% When monte carlo weight greater, we use the corresponding weights to
% determine false alarm.
falseAlarm = zeros(size(T));
%falseAlarm = zeros(size(weight));
falseAlarm(iFalseAlarm) = weight(iFalseAlarm);
% The overall probability of false
expPfa = sum(falseAlarm)/nTrial;
end

