function [T,iCFARwin] = CA_CFAR1(varargin)
% [T,iCFARwin] = CA_CFAR(Pfa,N,nGuardCell, signal) produces a vector T,
% that is the threshold produced by the Cell-Averaging CFAR (CACFAR)
% method from the vector, signal, that is a one-dimension of range
% power.
%
%
% [T,iCFARwin]=CA_CFAR(Pfa,N,signal) produces a vector T, that is the
% threshold produced by the Cell-Averaging CFAR (CA-CFAR) method from the
% vector, signal, that is a one-dimension of range power. nGuard-Cell
% defaults to one
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
% CFAR Computations and Parameters
a = N*(Pfa^(-1/N) - 1); % Threshold multiplier
% Number of Range Gates
L = length(signal);
% If range gate inputs is a vector, then perform CFAR along entire range
% gate:
if (min(size(signal)) == 1)
    % Organizing linear index into a matrix
    iWindow = 1:N + 1 + 2*nGuardCell;
    nCUT = L - N - 2*nGuardCell;
    shiftCFAR = (0:nCUT - 1).';
    index = repmat(iWindow, nCUT, 1) + repmat(shiftCFAR, 1, N + 1 + 2*nGuardCell);
    % Linear Index for the reference window
    iRefWin = index;
    i = N/2 + (1:2*nGuardCell + 1);
    iRefWin(:,i) = [];
    % Linear index for the CUT
    i = N/2 + nGuardCell + 1;
    iCUT = index(:,i);
    % Summing the content of each reference window
    B = sum(signal(iRefWin),2)/N;
    % Outuptted thershold
    T = zeros(L,1); %For CUT too small or too big for complete
    % reference window, we set T to zero
    T(iCUT) = a*B; % Multiplying threshold multiplier
    % If range gate input is a matrix, then treat each row as one experimental
    % trial with one test cell
else
    % Row index for the Cell Under Test (CUT)
    iCUT = N/2 + nGuardCell + 1;
    % Row index cor the Reference Cells
    iRefCell = 1:N + 2*nGuardCell + 1;
    notRefCell = N/2 + (1:1 + 2*nGuardCell);
    iRefCell(notRefCell) = [];
    % Calculating the threshold
    B = sum(signal(:,iRefCell),2)/N; % sum of each cell in reference win
    T = a*B; % Multiplying, a, the threshold mult
end
% Range Cells where Thrsholds are set:
iCFARwin = 1:L;
iCFARwin(T == 0) = [];
end