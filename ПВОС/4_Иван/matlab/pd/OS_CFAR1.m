function[T,iCFARwin] = OS_CFAR1(varargin)
% [T,iCFARwin]=OS_CFAR(Pfa,N,k,nGuardCell signal) produces a threshold T,
% from signal to limit to probability of false larm Pfa using Ordered
% Statistics CFAR (OS-CFAR). N is the number of reference cells in the
% CFAR window and k is the cell rank used to generate the threshold
% multiplier. nGuardCell is the number of guard cells to either side of
% CUT
%
% [T,iCFARwin]=OS_CFAR(Pfa,N,k,nGuardCell signal) produces a threshold T,
% from signal to limit to probability of false larm Pfa using Ordered
% Statistics CFAR (OS-CFAR). N is the number of reference cells in the
% CFAR window and k is the cell rank used to generate the threshold
% multiplier. nGuardCell, the number of guard cells to either side of CUT
% defaults to 1.
%
% INPUT:
% Pfa - the desired probability of False Alarm
%
% N - the total number of CFAR windows-- with N/2 in the leading window
% and N/2 in the lagging window.
%
% k - the index used in computation of the threshold multiplier used in
% computing the threshold, T.
%
% nGuardCell - the number of guard cells to either side of CUT
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
    % If "k" is omitted as an input, it defaults to round(0.7*N)
    case 4
        Pfa = varargin{1};
        N = varargin{2};
        k = varargin{3};
        signal = varargin{4};
        nGuardCell = 1;
        % If "k" is entered as an input:
    case 5
        Pfa = varargin{1};
        N = varargin{2};
        k = varargin{3};
        nGuardCell = varargin{4};
        signal = varargin{5};
    otherwise
        error('Unexpected Number of Inputs');
end
nLead = N/2; % size of lead window
nLag = N/2; % size of lag window
% Number of Range Gates
L = length(signal);
% The proper threshold multiplier is determined numerically:
a_set = 0.1:0.05:100;
P = zeros(length(a_set),1);
for i = 1:length(a_set)
    a = a_set(i);
    P(i) = k*nchoosek(N,k)*gamma(a + N - k + 1)*gamma(k)/gamma(a + N + 1);
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
    iRefWin = index;
    i = N/2 + (1:2*nGuardCell + 1);
    iRefWin(:,i) = [];
    %Linear index for the CUT
    i = N/2 + nGuardCell + 1;
    iCUT = index(:,i);
    % sorting and taking the kth value
    sortedPwrSignal = sort(signal(iRefWin),2);
    B = sortedPwrSignal(:,k);
    % Outuptted thershold
    T = zeros(L,1); %For CUT too small or too big for complete
    % reference window, we set T to zero
    T(iCUT) = a*B; %Multiplying threshold multiplier
    % If range gate input is a matrix, then treat each row as one experimental
    % trial with one test cell
else
    % Matrix index for the Cell Under Test (CUT)
    iCUT = N/2 + nGuardCell + 1;
    % Row index cor the Reference Cells
    iRefCell = 1:N + 2*nGuardCell + 1;
    notRefCell = N/2 + (1:1 + 2*nGuardCell);
    iRefCell(notRefCell) = [];
    % sorting and taking the kth value
    sortedPwrSignal = sort(signal(:,iRefCell),2);
    B = sortedPwrSignal(:,k);
    T = a*B; %Multiplying, a, the threshold mult
end
% Range Cells where Thrsholds are set:
iCFARwin = 1:L;
iCFARwin(T == 0) = [];
end