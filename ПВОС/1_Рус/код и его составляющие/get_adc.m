function out = get_adc(in,bits)
    out = quantize(int32(in),1,bits,0,'Floor','Saturate');
end