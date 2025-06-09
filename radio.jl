mutable struct RadioDownLink
    μ
    nbins
    from
    cp
    symb
    thru
    boi
    gb
    fo
    iqs
    flt
    mix
    fs

    function RadioDownLink(args...;kwargs...)
        μ, nbins, from, cp, symb, thru, boi, gb, fo, iqs, flt, mix, fs = radioDownLink(args...;kwargs...)
        return new(μ, nbins, from, cp, symb, thru, boi, gb, fo, iqs, flt, mix, fs)
    end
end

function radioDownLink(rdl::RadioDownLink)
    return rdl.μ, rdl.nbins, rdl.from, rdl.cp, rdl.symb, rdl.thru, rdl.boi, rdl.gb, rdl.fo, rdl.iqs, rdl.flt, rdl.mix, rdl.fs
end

function radioDownLink()
    return nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
end

Base.show(io::IO,rdl::RadioDownLink)=print(io,"mix=$(mixer_frequency(rdl)), boi=$(band_of_interest(rdl)) @ fs=$(sample_frequency(rdl)), $(from(rdl)):$(thru(rdl))")

function with_this_sample_frequency(k,o)
    μ= subcarrier_spacing_configuration(o)
    symbol_length_at_30M72Hz=symbol_length(μ)
    symbol_length_at_fs=number_of_bins(o)
    ratio=symbol_length_at_fs/symbol_length_at_30M72Hz
    nshifts=round(Int64,log2(ratio))
    result=k << nshifts
    return result
end

function radioDownLink(o::OranType3C; lowpassfilter=nothing, bw_table=BANDWIDTH_TABLE, sg_table=SIGNAL_GENERATION_TABLE)
    μ= subcarrier_spacing_configuration(o)
    nbins= number_of_bins(o)

    scs= subcarrier_spacing(o)
    nprbs= number_of_prbs(o)

    bw= bw_table[(scs=scs, nprbs=nprbs)]
    
    gb= in_units_of_7k5Hz( guardband_in_steps_of_1kHz(scs,bw) )
    boi= band_of_interest_in_units_of_7k5Hz(o)

    fo= frequency_offset_in_units_of_7k5Hz(o)

    subFrameIdMod2 = o.subFrameId % 2
    slotId= o.slotId
    symbolId= o.startSymbolId
    extended= isextended(o)

    from, length_of_cyclic_prefix, length_of_symbol= sg_table[(μ = μ, extended = extended, subFrameId = subFrameIdMod2, slotId = slotId, symbolId = symbolId)]

    from= with_this_sample_frequency(from, o)
    cp= with_this_sample_frequency(length_of_cyclic_prefix ,o)
    symb = with_this_sample_frequency(length_of_symbol ,o)

    iqs= copy(inphase_n_quadratures(o))

    lpf= lowpassfilter
    mix= nothing
    fs= nothing
    thru= nothing

    return μ, nbins, from, cp, symb, thru, boi, gb, fo, iqs, lpf, mix, fs
end

hertz(n)=7500n

subcarrier_spacing_configuration(scs::Int)=findfirst(e -> e==scs, [15,30,60])-1
subcarrier_spacing_configuration(rdl::RadioDownLink)=rdl.μ
subcarrier_spacing_configuration!(rdl::RadioDownLink,μ)=(rdl.μ=μ)

number_of_bins(rdl::RadioDownLink) = rdl.nbins
number_of_bins!(rdl::RadioDownLink,nbins)=(rdl.nbins=nbins)

inphase_n_quadratures(rdl::RadioDownLink)=rdl.iqs
inphase_n_quadratures!(rdl::RadioDownLink, iqs)=(rdl.iqs=iqs)

frequency_offset(rdl::RadioDownLink) = rdl.fo
frequency_offset!(rdl::RadioDownLink, fo) = (rdl.fo=fo)

mixer_frequency(rdl::RadioDownLink) = rdl.mix
mixer_frequency!(rdl::RadioDownLink, mix) = (rdl.mix=mix)

sample_frequency(rdl::RadioDownLink)= isnothing(rdl.fs) ? rdl.nbins << (rdl.μ+1) : rdl.fs
sample_frequency!(rdl::RadioDownLink,fs)= (rdl.fs=fs)

length_of_cyclic_prefix(rdl::RadioDownLink)=rdl.cp
length_of_cyclic_prefix!(rdl::RadioDownLink, cp)=(rdl.cp=cp)

from(rdl::RadioDownLink)=rdl.from
from!(rdl::RadioDownLink, from)=(rdl.from=from)

thru(rdl::RadioDownLink)=rdl.thru
thru!(rdl::RadioDownLink, thru)=(rdl.thru=thru)

band_of_interest(rdl::RadioDownLink)=rdl.boi
band_of_interest!(rdl::RadioDownLink, boi)=(rdl.boi=boi)

guardband(rdl::RadioDownLink)=rdl.gb
guardband!(rdl::RadioDownLink, gb)=(rdl.gb=gb)

length_of_symbol(rdl::RadioDownLink)=rdl.symb
length_of_symbol!(rdl::RadioDownLink, symb)=(rdl.symb=symb)

lowpassfilter(rdl::RadioDownLink)=rdl.flt
lowpassfilter!(rdl::RadioDownLink, flt)=(rdl.flt=flt)

halfbandfilter = lowpassfilter
halfbandfilter! = lowpassfilter!


static_chain(rdl::RadioDownLink)= rdl |>    prbs2bins |> phase_correction |> create_lowpassfilter |> 
                                            amplitude_correction  |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier |> 
                                            out_of_band_suppression



function bins_with_iq(iqs) 
    n= length(iqs)
    n2= n>>1
    m= nextpow(2,n)
    result=vcat(collect(m-n2+1:m),collect(1:n2))
    return result
end

@test bins_with_iq(1:10)==[12,13,14,15,16,1,2,3,4,5]

function prbs2bins(rdl::RadioDownLink)
    μ= subcarrier_spacing_configuration(rdl)
    nbins= number_of_bins(rdl)
    fo= frequency_offset(rdl)
    
    iqs= inphase_n_quadratures(rdl)
    bins= zeros(eltype(iqs), nbins)
    bins[bins_with_iq(iqs)] .= iqs
    
    mix = fo + length(iqs) << μ
    
    mixer_frequency!(rdl, mix)
    inphase_n_quadratures!(rdl, bins)

    return rdl
end

function  oscillator(fs, mix, range)
    ms = 0.001
    number_of_samples_in_2ms = hertz(fs)*2ms
    number_of_7500hz_periods_per_2ms = 15

    ω= number_of_7500hz_periods_per_2ms * mix
    result= exp.(1im*2π*(range)*ω/number_of_samples_in_2ms)
    return result
end

function phase_correction(rdl::RadioDownLink)
    fs = sample_frequency(rdl)
    mix = mixer_frequency(rdl)
    fr = from(rdl)
    cp = length_of_cyclic_prefix(rdl)
    iqs= inphase_n_quadratures(rdl)
    
    phase = first(oscillator(fs, mix, fr+cp:fr+cp))
    
    iqs ./= phase
    inphase_n_quadratures!(rdl, iqs)

    return rdl
end

function create_lowpassfilter(rdl::RadioDownLink)
    lpf= lowpassfilter(rdl)
    fs= sample_frequency(rdl)
    boi= band_of_interest(rdl)
    gb= guardband(rdl)

    lpf= isnothing(lpf) ? remezfind(boi/fs, boi/fs+gb/fs; Rs=db2amp(-26), Rp=db2amp(1.0)-1) |> FIRFilter : lpf

    lowpassfilter!(rdl, lpf)
    boi= boi + 2gb
    boi=band_of_interest!(rdl, boi)
    guardband!(rdl,nothing)

    return rdl
end

function amplitude_correction(rdl::RadioDownLink)
    lpf= lowpassfilter(rdl)
    b= coef(lpf)

    nbins= number_of_bins(rdl)
    z= exp.(1im *2π * (0.5:nbins-0.5) ./ nbins)
    ampl= abs.(tf(b;z=z))
    ampl[ampl .< 1e-12] .= 1
    
    iqs= inphase_n_quadratures(rdl)
    iqs ./= ampl
    inphase_n_quadratures!(rdl,iqs)

    return rdl
end

function bins2symbol(rdl::RadioDownLink)
    iqs=inphase_n_quadratures(rdl)
    time_domain_samples = ifft(iqs)
    inphase_n_quadratures!(rdl, time_domain_samples)

    return rdl
end

function with_cyclic_prefix(rdl::RadioDownLink)
    cp = length_of_cyclic_prefix(rdl)
    sym=length_of_symbol(rdl)
    fr=from(rdl)
    
    iqs=inphase_n_quadratures(rdl)
    cp_n_symbol = vcat(iqs[end-cp+1:end], iqs)
    inphase_n_quadratures!(rdl, cp_n_symbol)

    th = fr+cp+sym-1
    thru!(rdl,th)

    length_of_symbol!(rdl, nothing)

    return rdl
end

function shift_half_subcarrier(rdl::RadioDownLink)
    μ= subcarrier_spacing_configuration(rdl)
    fs= sample_frequency(rdl)
    fr= from(rdl)
    th= thru(rdl)

    nshifts= 1 << μ
    phase=oscillator(fs, nshifts, fr:th)

    iqs=inphase_n_quadratures(rdl)
    iqs .*= phase
    inphase_n_quadratures!(rdl,iqs)

    mix=mixer_frequency(rdl)
    mix -= nshifts
    mixer_frequency!(rdl,mix)

    return rdl
end

coef(f::FIRFilter)=f.h

filter_delay(rdl::RadioDownLink)= rdl |> lowpassfilter |> coef |> length |> r -> r >> 1

function out_of_band_suppression(rdl::RadioDownLink)
    lpf= lowpassfilter(rdl)
    delay= filter_delay(rdl)

    fr=from(rdl)
    th=thru(rdl)
    fr -= delay
    th -= delay
    from!(rdl,fr)
    thru!(rdl,th)

    iqs=inphase_n_quadratures(rdl)

    filtered_iqs= filt(lpf, iqs)
    lowpassfilter!(rdl, lpf)
    inphase_n_quadratures!(rdl, filtered_iqs)
    
    return rdl
end

function filter_w_meta!(rdl::RadioDownLink)
    result=RadioDownLink(rdl)

    th=thru(rdl)
    d=filter_delay(rdl)
    fr=th+d+1

    from!(result,fr)
    thru!(result,nothing)

    lowpassfilter!(rdl,nothing)
    return result
end

function upsample(rdl::RadioDownLink)
    fs=sample_frequency(rdl)
    fr= from(rdl)
    th= thru(rdl)
    iqs= inphase_n_quadratures(rdl)

    n=2(th-fr+1)
    fs *= 2
    fr *= 2

    th = fr + n -1

    iqs=reshape(iqs,1,length(iqs))
    zs=zeros(eltype(iqs),size(iqs))
    iqs=vcat(iqs,zs)[:]

    inphase_n_quadratures!(rdl,iqs)
    thru!(rdl, th)
    from!(rdl,fr)
    sample_frequency!(rdl,fs)

    return rdl
end

function mixer(rdl::RadioDownLink)
    fs= sample_frequency(rdl)
    mix= mixer_frequency(rdl)
    
    fr= from(rdl)
    th= thru(rdl)
    
    iqs=inphase_n_quadratures(rdl)
    
    ϕ=oscillator(fs,mix,fr:th)
    iqs .*= ϕ
    mix=0

    inphase_n_quadratures!(rdl,iqs)
    mixer_frequency!(rdl,mix)

    return rdl
end

function create_halfbandfilter(rdl::RadioDownLink)
    hbf= halfbandfilter(rdl)
    fs = sample_frequency(rdl)
    boi= band_of_interest(rdl)

    Ws=Wp=boi/fs
    hbf= isnothing(hbf) ? find_halfband(Wp,Ws) |> FIRFilter : hbf

    halfbandfilter!(rdl, hbf)

    return rdl
end

function flush(flt)
    n=2*filter_delay(flt)
    T=eltype(flt.flt.history)
    
    iqs=zeros(T,n)
    inphase_n_quadratures!(flt,iqs)
    
    fr=from(flt)
    th=fr+n-1
    thru!(flt,th)
    
    result=out_of_band_suppression(flt)
        
    return result
end

function mix_n_boi(a::RadioDownLink,b::RadioDownLink)
    a,b= mixer_frequency(a)<mixer_frequency(b) ? (a,b) : (b,a)

    mixa=mixer_frequency(a)
    mixb=mixer_frequency(b)
    boia=band_of_interest(a)
    boib=band_of_interest(b)

    mix= (mixa-boia>>1)>>1+(mixb+boib>>1)>>1
    boi= 2max(mix-mixa+boia>>1,mixb+boib>>1-mix)
    return mix, boi        
end

function align_end(a,b)
    th= min(thru(a),thru(b))

    if thru(b)>thru(a)
        a, b= b, a
    end

    x=nothing
    
    if thru(a)>thru(b)
        x=RadioDownLink(a)
        iqs=inphase_n_quadratures(x)
        fr=from(x)
        n=th-fr+1

        from!(x,th+1)
        inphase_n_quadratures!(x,iqs[n+1:end])
        thru!(a,th)
        inphase_n_quadratures!(a,iqs[1:n])
    end

    return a, b, x
end

function align_start(a,b)
    fr= min(from(a),from(b))
    th= thru(a)

    if th != thru(b)
        error("thru messed up")
    end
    
    for d in [a,b]
        iqs=inphase_n_quadratures(d)
        n=max(from(d)-fr,0)
        zs=zeros(eltype(iqs),n)
        iqs=vcat(zs,iqs)
        inphase_n_quadratures!(d,iqs)
        from!(d,fr)
    end
    
    return a, b
end

function mix_n_merge(a,b)
    a, b, x= align_end(a, b)
    a, b= align_start(a, b)
    
    fs=sample_frequency(a)
    fr= from(a)
    th= thru(a)
    mix,boi=mix_n_boi(a,b)
    
    mixed_n_merged= similar(inphase_n_quadratures(a)) .=0
 
    for d in [a,b]
        iqs= inphase_n_quadratures(d)
        old= mixer_frequency(d)
        phase= oscillator(fs, old-mix, fr:th)

        mixed_n_merged .+= iqs .* phase
    end

    result=RadioDownLink()
    sample_frequency!(result,fs)
    from!(result,fr)
    thru!(result,th)
    band_of_interest!(result,boi)
    mixer_frequency!(result,mix)
    inphase_n_quadratures!(result,mixed_n_merged)

    return result, x
end

suppress_mirror(rdl::RadioDownLink)=out_of_band_suppression(rdl::RadioDownLink)

function output_buffer!(y,data)
    iqs= inphase_n_quadratures(data)
    
    fr=from(data)
    th=thru(data)
    indexes=collect(fr:th)
    inrange=findall(e -> e in 0:length(y)-1, indexes)
    y[indexes[inrange].+1]+=iqs[inrange]
    return y
end

freq_n_from(o::OranType3C)= o |> RadioDownLink |> prbs2bins |> r -> (frequency_offset(r), from(r))

function process_orans(orans)
    datas=[]
    flts=[]
    
    for o in orans
        fo,fr=freq_n_from(o)
    
        p=findfirst(e -> frequency_offset(e) == fo && from(e)==fr, flts)
    
        if isnothing(p)
            data= RadioDownLink(o) |> static_chain
        else
            flt=lowpassfilter(flts[p])
            deleteat!(flts,p)       
            data=RadioDownLink(o; lowpassfilter=flt) |> static_chain
        end
    
        flt=filter_w_meta!(data)
        
        push!(datas,data)
        push!(flts,flt)
    end
    
    while length(flts)>0
        flt=flts[1]
        deleteat!(flts,1) 
        data=flush(flt)
        push!(datas,data)
    end

    return datas,flts
end

function can_merge(a,b)
    isnothing(a) && return false
    isnothing(b) && return false
    (from(a)>thru(b) || from(b)>thru(a)) && return false

    boi=last(mix_n_boi(a,b))
    fs=sample_frequency(a)
    result= boi/fs < 0.97
    return result    
end

function find_data(datas,fs, t)
    isempty(datas) && return nothing
    ok= sample_frequency.(datas) .== fs
    ok .&= (from.(datas).<=t)
    v= (ok.==0)*+Inf .+ mixer_frequency.(datas)
    p= last(findmin(v))
    isinf(v[p]) && return nothing
    return p
end

function find_data!(datas,fs, t)
    p=find_data(datas,fs,t)
    if isnothing(p)
        return nothing
    else
        result=datas[p]
        deleteat!(datas,p)
        return result
   end
end

function find_next!(datas,fs, t, a)
    isempty(datas) && return nothing
    ok= sample_frequency.(datas) .== fs
    ok .&= (from.(datas).<=t)
    fr=from(a)
    th=thru(a)
    ok .&= (thru.(datas).<fr .|| from.(datas).>th) .==0
    
    v= (ok.==0)*+Inf .+ mixer_frequency.(datas)
    p= last(findmin(v))
    isinf(v[p]) && return nothing
    result=datas[p]
    deleteat!(datas,p)

    return result
end

function find_flt!(flts, fs, t)
    obsoletes = sample_frequency.(flts) .== fs .&& from.(flts) .< t
    p= findfirst(obsoletes)
    if isnothing(p)
        return nothing
    else
        result=flts[p]
        deleteat!(flts,p)
        return result
    end
end

sample_frequencies(fs_in,fs_out)=2 .^ (Int(log2(fs_in)):Int(log2(fs_out)))
issomething(u) = ! isnothing(u)

function process_data!(y,datas, flts, fs_in, fs_out)
    fs=fs_in
    fs_in_time=0
    t=fs_in_time
    while !isempty(datas)
        println("t:$t")
        for fs in sample_frequencies(fs_in,fs_out)
            while (a=find_data!(datas,fs,t)) |> issomething
            
                b=find_next!(datas,fs,t,a)
            
                if can_merge(a,b)
                    print("$a + $b -> ")
                    a,b=mix_n_merge(a,b)
                    push!(datas,a)
                    push!(datas,b)
                    println("$a + $b")
                elseif fs < fs_out
                    issomething(b) && push!(datas,b)
                    print("$a -> ")
                    data=upsample(a)
                    p=findfirst(e -> mixer_frequency(e) == mixer_frequency(data) && from(e)==from(data), flts)
                    
                    if isnothing(p)
                        data=create_halfbandfilter(data) |> suppress_mirror
                    else
                        flt=lowpassfilter(flts[p])
                        deleteat!(flts,p)       
                        halfbandfilter!(data,flt)
                        data=suppress_mirror(data)
                    end
                    flt=filter_w_meta!(data)
                    push!(flts,flt)
                    push!(datas,data)
                    println("$data")
                elseif fs == fs_out
                    issomething(b) && push!(datas,b)
                    print("$a -> ")
                    a=mixer(a) 
                    output_buffer!(y, a)
                    println("out($a)")
                end
            end
            
            while (flt=find_flt!(flts, fs, t)) |> issomething
                data=flush(flt)
                push!(datas,data)
            end
            
            fs *= 2
            t *= 2
        end
        
        fs_in_time+=512
        t=fs_in_time
        fs=fs_in
    end
end
