function dem_signal = get_dem(adc_signal,nco_signal,noba,nobg,dem_type)
    switch dem_type
        case 'fixed'
            dem_signal = adc_signal.*nco_signal;
            dem_signal = sfi(dem_signal,noba + nobg,0);
        case 'single'
            dem_signal = single(adc_signal).*single(nco_signal);
            dem_signal = sfi(dem_signal,noba + nobg,0);
        case 'double'
            dem_signal = double(adc_signal).*double(nco_signal);
            dem_signal = sfi(dem_signal,noba + nobg,0);  
    end
end