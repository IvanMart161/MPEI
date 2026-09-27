function [nco_signal,nco_gain] = get_nco(nobg,fg,t,nco_type,amv,pmv)

    % Get signal nco
    nco_real = cos(2*pi*fg*t).';
    % Get Phase mismatch IQ
    nco_imag = sin(2*pi*fg*t + deg2rad(pmv)).';
    % Get am[lide mismatch IQ
    nco_real = max(real(nco_real))*db2mag(amv)*nco_real;
    % Calulate nco
    nco = nco_real - 1i*nco_imag;
    
    % Calcilate gain nco
    nco_gain = 2^(nobg - 1) - 1;

    % Fromat conversion
    switch nco_type
        case 'fixed'
            nco = nco_gain.*nco;
            nco_signal = sfi(nco,nobg,0);
        case 'single'
            nco_signal = single(nco);
        case 'double'
            nco_signal = double(nco);
    end
end