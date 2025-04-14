function symbol_with_cp(y,fs,o::OranType3C)
    r=fs ÷ sample_frequency(o |> oran2prbs)

    fr=from(o |> oran2prbs |> prbs2bins |> bins2symbol |> with_cyclic_prefix)
    th=thru(o |> oran2prbs |> prbs2bins |> bins2symbol |> with_cyclic_prefix)

    n=r*(th-fr+1)
    fr *= r
    th=fr+n-1

    iqs=y[fr+1:th+1]

    return (oran=o, fs=fs, iqs=iqs)
end

function odd_shift(y)
    oran=y.oran
    offset=frequency_offset(oran)
    fs=sample_frequency(y)
    r=fs ÷ sample_frequency(oran |> oran2prbs)
    
    cp=cyclic_prefix(oran |> oran2prbs) * r

    iqs=inphase_n_quadratures(y)
    n=length(iqs)-cp
    
    if isodd(offset)
        println("odd shift active")
        ω=exp.(-π*im*(0-cp:length(iqs)-cp-1)/n)
        iqs .*= ω
    else
        println("odd shift in-active")
    end
    
    return (oran=oran, fs=fs, iqs=iqs)
end

function without_cyclicprefix(y)
    oran=y.oran
    r=fs ÷ sample_frequency(oran |> oran2prbs)
    iqs=inphase_n_quadratures(y)

    cp=cyclic_prefix(oran |> oran2prbs)
    cp2=cp>>1

    cp2 *= r
    iqs=iqs[cp2+1:end-cp2]
    iqs=vcat(iqs[cp2+1:end],iqs[1:cp2])

    #cp *= r
    #iqs=iqs[cp+1:end]

    return (oran=oran, fs=fs, iqs=iqs)
end

function time2bins(y)
    oran=y.oran
    fs=sample_frequency(y)
    iqs=inphase_n_quadratures(y)
    iqs=fft(iqs)
    
    return (oran=oran, fs=fs, iqs=iqs)
end

function shift_frequency_offset(y)
    oran=y.oran
    offset2=-(frequency_offset(oran) >> 1)

    iqs=inphase_n_quadratures(y)

    iqs=circshift(iqs, offset2)
    return (oran=oran, fs=fs, iqs=iqs)
end

# o=a1
# fs=sample_frequency(data)
# s=symbol_with_cp(y,fs,o) |> odd_shift |> without_cyclicprefix   |> time2bins |> shift_frequency_offset
# expected = o |>  inphase_n_quadratures
# result=s |>  inphase_n_quadratures

# @test all(isapprox.(result[1:length(expected)],expected, atol=0.01))


