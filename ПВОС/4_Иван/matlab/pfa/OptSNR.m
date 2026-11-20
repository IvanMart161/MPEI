function [SNRdB] = OptSNR(Pd,Pfa)
    SNR = log(Pfa./Pd)./log(Pd);
    SNRdB = pow2db(SNR);
end
