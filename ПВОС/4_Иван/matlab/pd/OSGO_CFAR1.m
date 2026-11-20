function[T,iCFARwin] = OSGO_CFAR1(varargin)
% [T,iCFARwin]=OSGO_CFAR(Pfa,N,k,signal) produces a vector T, that is the
% threshold produced by the Order Statistics Greatest-OF CFAR (OSGO-CFAR)
% method from the vector, signal, based on the k-th ranked index.
%
% [T,iCFARwin]=OSGO_CFAR(Pfa,N,k,signal) produces a vector T, that is the
% threshold produced by the Order Statistics Greatest-OF CFAR (OSGO-CFAR)
% method from the vector, signal. "k" is defaulted to round(0.6*N/2)
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
    % If "k" is omitted as an input, it defaults to round (0.7*N)
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

% Number of Range Gates
L = length(signal);
% Numerically determining the threshold multiplier, a:
a_set = 1:0.05:100;
Pfa_set = ones(size(a_set));
i = (0:N/2 - k)';
j = (0:N/2 - k)';
I = i*ones(1,length(i));
J = ones(length(j),1)*j';
for ia = 1:length(a_set)
    a = a_set(ia);
    A = 2*k^2*nchoosek(N/2,k)^2;
    B = factorial((N/2 - k)*ones(length(j)))./factorial(J)./factorial(N/2 - k - J);
    C = factorial((N/2 - k)*ones(length(j)))./factorial(I)./factorial(N/2 - k - I);
    D = (-1).^(N - 2*k - J - I)./(N/2 - I);
    E = gamma(N - J - I).*gamma(a + 1)./gamma(N - J - I + a + 1);
    Pfa_set(ia) = A*sum(sum(B.*C.*D.*E));
end
a = interp1(Pfa_set,a_set,Pfa);
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
    iLead = iRefWin(:,1:N/2);
    iLag = iRefWin(:,N/2 + 1:N);
    % Linear index for the CUT
    i = N/2 + nGuardCell + 1;
    iCUT = index(:,i);
    % sorting and taking the kth value  
    sortedLead = sort(signal(iLead),2);
    sortedLag = sort(signal(iLag),2);
    B = max(sortedLead(:,k),sortedLag(:,k));
    % Outuptted thershold
    T = zeros(L,1); %For CUT too small or too big for complete
    % reference window, we set T to zero
    T(iCUT) = a*B; % Multiplying threshold multiplier
    % If range gate input is a matrix, then treat each row as one experimental
    % trial with one test cell
else
    % Matrix index for the Cell Under Test (CUT)
    iCUT = N/2 + nGuardCell + 1;
    % Row index cor the Reference Cells
    iLead = 1:N/2;
    iLag = N/2 + 2*nGuardCell + 2:N + 1 + 2*nGuardCell;
    sortedLead = sort(signal(:,iLead),2);
    sortedLag = sort(signal(:,iLag),2);
    B = max(sortedLead(:,k),sortedLag(:,k));
    T = a*B; %Multiplying, a, the threshold mult
end
% Range Cells where Thrsholds are set:
iCFARwin = 1:L;
iCFARwin(T == 0) = [];
end