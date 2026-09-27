function [f,s] = get_spectrum(in,fs,weighting)
    f = linspace(-fs/2,fs/2,length(in));
    if (weighting)
        s = fftshift(fft(double(in).*hann(length(in))));
    else
        s = fftshift(fft(double(in)));
    end
end