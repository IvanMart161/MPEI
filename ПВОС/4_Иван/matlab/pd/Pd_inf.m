function [Pd] = Pd_inf(Pfa,SNRdB)
    SNR = 10.^(SNRdB./10);
    Pd = Pfa.^(1./(1 + SNR));
end

