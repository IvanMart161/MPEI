function[T,iCFARwin] = SOCA_CFAR1(varargin)
% [T,iCFARwin]=SOCA_CFAR(Pfa,N,nGuardCell, signal) produces a vector T,
% that is the threshold produced by the Smalelst-Of Cell-Averaging CFAR
% (SOCA-CFAR) method from the vector, signal, that is a one-dimension of
% range power.
%
%
% [T,iCFARwin]=SOCA_CFAR(Pfa,N,signal) produces a vector T, that is the
% threshold produced by the Smallest-Of Cell-Averaging CFAR (SOCA-CFAR
% method from the vector, signal, that is a one-dimension of range power.
% nGuardCell defaults to one
%
% INPUT:
% Pfa - the desired probability of False Alarm
%
% N - the total number of CFAR windows-- with N/2 in the leading window
% and N/2 in the lagging window.
%
% nGuardCell- the number of guard cells to either side of the test cell
%
% Signal - the power of the recieved waveform after having passed through
% a square law detector. It is upon this that the threshold is set.
%
% OUTPUT:
% T is the threshold set according to the functional inputs.
%
% iCFARwin are the indicies corresponding to the range gates where a CFAR
% threshold is produced. When the first few and last few range bins are
% used as test cells, they don't allow enough room for the entire CFAR
% reference window. Subsequently, they are ignored and scipped.

switch nargin
    case 3
        Pfa = varargin{1};
        N = varargin{2};
        signal = varargin{3};
        nGuardCell = 1;
    case 4
        Pfa = varargin{1};
        N = varargin{2};
        nGuardCell = varargin{3};
        signal = varargin{4};
    otherwise
        error('Unexpected Number of Inputs')
end
% Number of Range Gates
L = length(signal);
% Finding the threshold multiplier numerically
a_set = 0.1:0.1:150;
P = zeros(length(a_set),1);
for i = 1:length(a_set)
    a = a_set(i);
    sum_term = zeros(N/2,1);
    for k = 0:N/2 - 1
        sum_term(k + 1,1) = nchoosek(N/2 - 1 + k,k)*(2 + a/(N/2))^-k;
    end
    sum_term = sum(sum_term);
    P(i,1) = 2*(2 + a/(N/2))^(-N/2)*sum_term;
end
a = interp1(P,a_set,Pfa);
% If range gate inputs is a vector, then perform CFAR along entire range
% gate:
if (min(size(signal)) == 1)
    % Organizing linear index into a matrix
    iWindow = 1:N + 1 + 2*nGuardCell;
    nCUT = L - N - 2*nGuardCell;
    shiftCFAR = (0:nCUT - 1).';
    index = repmat(iWindow,nCUT,1) + repmat(shiftCFAR,1,N + 1 + 2*nGuardCell);
    % Linear Index for the reference window
    iLead = index(:,1:N/2);
    iLag = index(:,N/2 + 2*nGuardCell + 2:end);
    % Linear index for the CUT
    i = N/2 + nGuardCell + 1;
    iCUT = index(:,i);
    % Performaing SO-CFAR
    leadWin = signal(iLead);
    sumLead = sum(leadWin,2);
    lagWin = signal(iLag);
    sumLag = sum(lagWin,2);
    B = min(sumLead,sumLag)/(N/2);
    % Outuptted thershold
    T = zeros(L,1); %For CUT too small or too big for complete
    % reference window, we set T to zero
    T(iCUT) = a*B; %Multiplying threshold multiplier
    % If range gate input is a matrix, then treat each row as one experimental
    % trial with one test cell
else 
    % Row index cor the Reference Cells
    leadWin = signal(:,1:N/2);
    sumLead = sum(leadWin,2);
    lagWin = signal(:,N/2 + 2*nGuardCell + 2:end);
    sumLag = sum(lagWin,2);
    B = min(sumLead,sumLag)/(N/2);
    % Calculating the threshold
    T = a*B; %Multiplying, a, the threshold mult
end
% Range Cells where Thrsholds are set:
iCFARwin = 1:L;
iCFARwin(T == 0) = [];
end