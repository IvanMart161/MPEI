function [P] = simRangeData(N,iN,S,iS)
% [P] = simRangeData(N,iN,S,iS) produces a vector P that corresponds to
% the power over range gates with complex gaussian white noise and
% Swerling I/II targets.
%
% Vector N determines the power of the complex gaussian white noise across
% Range Cells controlled by indices iN.
%
% Vector S determines the power of a Swerling I/II target in Range Cells
% dictated by vector iS

pwrSig = 10.^(S/10);   % Converting signal power from dB to linear
pwrNoise = 10.^(N/10); % Converinting noise power from dB to linear
nRangeGate = iN(end);  % The total number of range gates

%% GENERATING NOISE
voltNoise_I = zeros(nRangeGate,1); % Pre-allocating for speed
voltNoise_Q = zeros(nRangeGate,1); % Pre-allocating for speed
I = length(iN); % Number of range gate segments with different noise power
iLead = 1;
% This loop produces the voltage waveform for the gaussian white noise.
% Every range gate segment between iLead and iLag will have its own power of
% N(i)
for i = 1:I
    iLag = iN(i);
    voltNoise_I(iLead:iLag) = sqrt(pwrNoise(i)/2)*randn(iLag - iLead + 1,1); % voltNoise_I Channel
    voltNoise_Q(iLead:iLag) = sqrt(pwrNoise(i)/2)*randn(iLag - iLead + 1,1); % voltNoise_Q Channel
    iLead = iLag;
end

noise = voltNoise_I + 1i*voltNoise_Q; %Voltage waveform of the noise

%% GENERATING TARGET RETURN
signal = zeros(size(noise)); % Pre-allocating for speed
I = length(iS); % The total number of signal returns
for i = 1:I
    % Power of the signal follows an exponential distribution
    pwrSwerling = exprnd(pwrSig(i),1,1);
    % Voltage waveform fo the signal
    signal(iS(i)) = sqrt(pwrSwerling/2) + 1i*sqrt(pwrSwerling/2);
end

%% GENERATING OVERALL SIGNAL
data = signal + noise; % Overall voltage waveform
P = abs(data).^2; % Power of the overall return square low detector
end