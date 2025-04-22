function postprocess(iqs, oran, fs)
    prbs= oran |> RadioDownLink
    r=fs ÷  sample_frequency(prbs)

    mx= -frequency_offset(oran) << subcarrier_spacing_configuration(oran)
    ω= oscillator(fs,mx,0:length(iqs)-1)
    result= ω .* iqs
    
    cp=length_of_cyclic_prefix(prbs)
    fr=r*from(prbs)+r*cp
    ϕ=first(oscillator(fs,mx,fr:fr))
    result ./= ϕ
    
    boi=band_of_interest(prbs)
    gb=guardband(prbs)
    n_gb=ceil(Int64,gb/fs*length(result))
    n_boi=ceil(Int64,boi/fs*length(result))
    
    result=fft(result)
    result[n_boi+n_gb+1:end-n_gb] .= 0
    result=ifft(result)
        
    cp2= length_of_cyclic_prefix(prbs) >> 1
    fr= r*from(prbs)+r*cp2
    th= r*from(prbs)+r*cp2+r*length_of_symbol(prbs)-1
    result=result[fr+1:th+1]

    delay=zeros(length(result))
    delay[r*cp2+1]=1
    ω=fft(delay)
    result=fft(result) ./ ω

    expected=oran.iqs
    result=result[1:length(expected)]
        
    return result, expected
end


