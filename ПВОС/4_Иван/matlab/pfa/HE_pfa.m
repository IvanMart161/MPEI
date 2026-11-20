clear all; close all; clc

%% Some properties
CA = 1; OS = 2; GOCA = 3; SOCA = 4; OSGO = 5; % CFAR numbering
CFARselection = OSGO;
Pfa_set = 10.^-(1:0.25:7); % Desired probability of false alarm
nPfa = length(Pfa_set);
nRangeGate = 10^6; % Number of trials to run
N = 32;
nGuardCell = 1; % number of Guard Cells
noisePwr = 1;

k_OS = round(3*N/4); % Optimal OS-CFAR k-rank (Richards)
k_OSGO = round(5*N/2/6); % Optimal k_OSGO k_rank (Elia Fustel)

% Legend name
CFARtype{CA}   = 'CA-CFAR';
CFARtype{OS}   = 'OS-CFAR';
CFARtype{GOCA} = 'GOCA-CFAR';
CFARtype{SOCA} = 'SOCA-CFAR';
CFARtype{OSGO} = 'OSGO-CFAR';

%% Theoretical Threshold Multiplier
pwrSignal = exprnd(noisePwr,nRangeGate,1);
a = zeros(size(Pfa_set));

switch CFARselection
    case CA
        a = N*(Pfa_set.^(-1/N)-1); % 7.17
    case OS
        k = k_OS;
        a_set = 0.1:0.05:100;
        P = zeros(length(a_set),1);
        for i = 1:length(a_set)
            a = a_set(i);
            P(i) = k*nchoosek(N,k)*gamma(a + N - k + 1)*gamma(k)/gamma(a + N + 1);  % 7.49
        end
        for i = 1:length(Pfa_set)
            Pfa = Pfa_set(i);
            a(i) = interp1(P,a_set,Pfa);
        end
    case SOCA
        a_set = 0.1:0.1:150;
        P = zeros(length(a_set),1);
        for i = 1:length(a_set)
            a = a_set(i);
            sum_term = zeros(N/2,1);
            for k = 0:N/2 - 1
                sum_term(k + 1,1) = nchoosek(N/2 - 1 + k,k)*(2 + a/(N/2))^-k;
            end
            sum_term = sum(sum_term);
            P(i,1) = 2*(2 + a/(N/2))^(-N/2)*sum_term;   % 7.35
        end
        for i = 1:length(Pfa_set)
            Pfa = Pfa_set(i);
            a(i) = interp1(P,a_set,Pfa);
        end
    case GOCA
        k = k_OSGO;
        a_set = 0.1:0.05:50;
        P = zeros(length(a_set),1);
        for i = 1:length(a_set)
            a = a_set(i);
            sum_term = zeros(N/2,1);
            for k = 0:N/2 - 1
                sum_term(k + 1,1) = nchoosek(N/2 - 1 + k,k)*(2 + a/(N/2))^-k;
            end
            sum_term = sum(sum_term);
            P(i,1) = 2*((1 + a/(N/2))^(-N/2) - (2 + a/(N/2))^(-N/2)*sum_term);  % 7.37
        end
        for i = 1:length(Pfa_set)
            Pfa = Pfa_set(i);
            a(i) = interp1(P,a_set,Pfa);
        end
    case OSGO
        k = k_OSGO;
        % Numerically determining the threshold multiplier, a:
        a_set = 1:0.05:100;
        P = zeros(size(a_set));
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
            P(ia) = A*sum(sum(B.*C.*D.*E));
        end
        for i = 1:length(Pfa_set)
            Pfa = Pfa_set(i);
            a(i) = interp1(P,a_set,Pfa);
        end
end % for CFARselection

%% CFAR Thresholding
threshold = zeros(nRangeGate,nPfa);
h = waitbar(0,'Calculating...');

for i = 1:nPfa
    waitbar(i/nPfa);
    Pfa = Pfa_set(i);
    switch (CFARselection)
        case CA
            [threshold(:,i) iCFAR] = CA_CFAR1(Pfa,N,nGuardCell,pwrSignal);
        case SOCA
            [threshold(:,i) iCFAR] = SOCA_CFAR1(Pfa,N,nGuardCell,pwrSignal);
        case GOCA
            [threshold(:,i) iCFAR] = GOCA_CFAR1(Pfa,N,nGuardCell,pwrSignal);
        case OS
            [threshold(:,i) iCFAR] = OS_CFAR1(Pfa,N,k,nGuardCell,pwrSignal);
        case OSGO
            [threshold(:,i) iCFAR] = OSGO_CFAR1(Pfa,N,k_OSGO,nGuardCell,pwrSignal);
    end % switch CFARselection
end % for i = 1:length(Pfa_set)

close(h);
expPfa = monteCarloFalseAlarm(threshold(iCFAR,:),noisePwr);

%% Plotting 
% CFAR plot color
clr = zeros(5,3);
clr(CA,:) = [0 0.447 0.741];        % blue
clr(OS,:) = [0.85 0.325 0.098];     % red
clr(GOCA,:) = [0.929 0.694 0.125];  % yellow
clr(SOCA,:) = [0.494 0.184 0.556];  % violet
clr(OSGO,:) = [0.466 0.674 0.188];  % green

figure
    set(gca,'FontName','Times New Roman','FontSize',24);
    semilogy(a,Pfa_set,'-','color',clr(CFARselection,:),'linewidth',2);
    hold on
    semilogy(a,expPfa,'*','color',clr(CFARselection,:),'linewidth',2);
    grid on
    legend('Theoretical P_F_A','Experimental P_F_A')
    title(['Theoretical vs Experimental P_F_A - \color{blue}' CFARtype{CFARselection}]);
    xlabel('Threshold Multiplier');
    ylabel('Probability of False Alarm')
    set(gca,'FontName','Times New Roman','FontSize',24);