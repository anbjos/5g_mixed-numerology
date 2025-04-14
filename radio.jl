guardband(scs, bw; table=MIN_GUARDBAND_KHZ)=round(Int64, table[(scs=scs,bw=bw)]/7.5)
subcarrier_spacing(μ::Int)=15<<μ

function oran2prbs(o::OranType3C; lowpassfilter=nothing, bw_table=BANDWIDTH_TABLE, sg_table=SIGNAL_GENERATION_TABLE)
    μ=subcarrier_spacing_configuration(o)
    nbins=number_of_bins(o)

    scs=subcarrier_spacing(o)
    nprbs=number_of_prbs(o)

    bw= bw_table[(scs=scs, nprbs=nprbs)]
    
    gb=guardband(scs,bw)
    boi=band_of_interest(o)

    fo=frequency_offset_7k5Hz(o)

    subFrameIdMod2 = o.subFrameId % 2
    slotId=o.slotId
    symbolId=o.startSymbolId
    extended=isextended(o)

    sgt=sg_table[(μ = μ, extended = extended, subFrameId = subFrameIdMod2, slotId = slotId, symbolId = symbolId)]

    fs_adjust=convert(Int64,log2(nbins/symbol_length(μ)))

    from=sgt.from << fs_adjust
    cp=sgt.cp  << fs_adjust
    symb = sgt.symbol  << fs_adjust

    iqs=copy(inphase_n_quadratures(o))

    return (μ=μ, nbins=nbins, from=from, cp=cp, symbol=symb, boi=boi, gb=gb, fo=fo, iqs=iqs, lpf=lowpassfilter)
end

hertz(n)=7500n

import Base: getindex

function getindex(nt::NamedTuple, key::Symbol, default)
    return haskey(nt, key) ? nt[key] : default
end

subcarrier_spacing_configuration(scs::Int)=findfirst(e -> e==scs, [15,30,60])-1
subcarrier_spacing_configuration(nt::NamedTuple)=nt.μ

sample_frequency(data::NamedTuple)= haskey(data, :fs) ? data.fs : data.nbins << (data.μ+1)
number_of_bins(data::NamedTuple) = data.nbins

guardband(data::NamedTuple)=getindex(data,:gb,0) 
band_of_interest(data::NamedTuple)=data.boi
mixer_frequency(data::NamedTuple)=data.mix
frequency_offset(data::NamedTuple) = data.fo

from(data::NamedTuple)=data.from
thru(data::NamedTuple)=data.thru

length_of_cyclic_prefix(data::NamedTuple)=data.cp
length_of_symbol(data::NamedTuple)=data.symbol   

lowpassfilter(data::NamedTuple)=data.lpf

inphase_n_quadratures(data::NamedTuple)=data.iqs

function iqmap(iqs)
    n= length(iqs)
    n2= n>>1
    m= nextpow(2,n)
    result=vcat(collect(m-n2+1:m),collect(1:n2))
    return result
end

@test iqmap(1:10)==[12,13,14,15,16,1,2,3,4,5]

function prbs2bins(data)
    μ= subcarrier_spacing_configuration(data)
    nbins= number_of_bins(data)
    lpf= lowpassfilter(data)

    boi= band_of_interest(data)
    gb= guardband(data)
    fo= frequency_offset(data)
    
    fr= from(data)
    cp= length_of_cyclic_prefix(data)
    sym= length_of_symbol(data)
    
    iqs= inphase_n_quadratures(data)

    bins= zeros(eltype(iqs), nbins)
    bins[iqmap(iqs)] .= iqs

    mix = fo + length(iqs) << μ

    return (μ=μ, nbins=nbins, from=fr, cp=cp, symbol=sym, boi=boi, gb=gb, mix=mix, lpf=lpf, iqs=bins)
end

function  oscillator(fs, mix, range)
    from= first(range)
    thru= last(range)

    ms= 0.001
    number_of_samples_in_2ms= hertz(fs)*2ms
    number_of_7500hz_periods_per_2ms=15

    ω= number_of_7500hz_periods_per_2ms * mix
    result= exp.(1im*2π*(from:thru)*ω/number_of_samples_in_2ms)
    return result
end

function phase_correction(data)
    μ=subcarrier_spacing_configuration(data)
    nbins=number_of_bins(data)
    fs=sample_frequency(data)
    boi=band_of_interest(data)
    gb=guardband(data)
    lpf=lowpassfilter(data)
    mix=mixer_frequency(data)
    
    fr=from(data)
    cp=length_of_cyclic_prefix(data)
    sym=length_of_symbol(data)
    
    phase = first(oscillator(fs, mix, fr+cp:fr+cp))
    
    iqs= inphase_n_quadratures(data) ./ phase
    
    return (μ=μ, nbins=nbins, from=fr, cp=cp, symbol=sym, boi=boi, gb=gb, mix=mix, iqs=iqs, lpf=lpf)
end


function create_lowpassfilter(data)
    μ=subcarrier_spacing_configuration(data)
    nbins=number_of_bins(data)
    fs=sample_frequency(data)
    boi=band_of_interest(data)
    gb=guardband(data)
    lpf=lowpassfilter(data)
    mix=mixer_frequency(data)

    lpf = isnothing(lpf) ? remezfind(boi/fs, boi/fs+gb/fs; Rs=db2amp(-40)) |> FIRFilter : lpf.flt

    fr=from(data)
    cp=length_of_cyclic_prefix(data)
    sym=length_of_symbol(data)

    iqs= inphase_n_quadratures(data)

    return (μ=μ, nbins=nbins, from=fr, cp=cp, symbol=sym, boi=boi, gb=gb, mix=mix, iqs=iqs, lpf=lpf)
end

function amplitude_correction(data)
    μ=subcarrier_spacing_configuration(data)
    nbins=number_of_bins(data)
    boi=band_of_interest(data)
    gb=guardband(data)
    lpf=lowpassfilter(data)
    mix=mixer_frequency(data)
    
    fr= from(data)
    cp= length_of_cyclic_prefix(data)
    sym= length_of_symbol(data)
    
    iqs= inphase_n_quadratures(data)

    z=exp.(1im *2π * (0.5:nbins-0.5) ./ nbins)
    b=coef(lpf)
    
    ampl=abs.(tf(b;z=z))
    ampl[ampl .< 1e-12] .= 1
    
    iqs ./= ampl
        
    return (μ=μ, nbins=nbins, from=fr, cp=cp, symbol=sym, boi=boi, gb=gb, mix=mix, iqs=iqs, lpf=lpf)
end

function bins2symbol(data)
    μ=subcarrier_spacing_configuration(data)
    nbins=number_of_bins(data)
    boi=band_of_interest(data)
    gb=guardband(data)
    lpf=lowpassfilter(data)
    mix=mixer_frequency(data)
    
    fr = from(data)
    cp = length_of_cyclic_prefix(data)
    sym = length_of_symbol(data)
    
    iqs=inphase_n_quadratures(data)

    time_domain_samples = ifft(iqs)

    return (μ=μ, nbins=nbins, from= fr, cp= cp, symbol= sym, boi=boi, gb=gb, mix=mix, lpf=lpf, iqs=time_domain_samples)
end

function with_cyclic_prefix(data)
    μ= subcarrier_spacing_configuration(data)
    nbins= number_of_bins(data)
    boi= band_of_interest(data)
    gb= guardband(data)
    mix= mixer_frequency(data)
    lpf= lowpassfilter(data)
    
    fr = from(data)
    cp = length_of_cyclic_prefix(data)
    sym = length_of_symbol(data)
    thru = fr+cp+sym-1
    
    iqs=inphase_n_quadratures(data)

    samples = vcat(iqs[end-cp+1:end], iqs)

    return (μ=μ, nbins=nbins, from=fr, thru=thru, boi=boi, gb=gb, mix=mix, lpf=lpf, iqs=samples)
end

function shift_half_subcarrier(data)
    μ= subcarrier_spacing_configuration(data)
    fs=sample_frequency(data)
    boi=band_of_interest(data)
    gb=guardband(data)
    mix=mixer_frequency(data)
    lpf= lowpassfilter(data)
    
    fr=from(data)
    th=thru(data)
    
    iqs=inphase_n_quadratures(data)

    nshifts= 1 << μ
    phase=oscillator(fs, nshifts, fr:th)
    iqs .*= phase
    #println("shs")
    #println(phase[1:3])
    
    #println("shs $mix, $nshifts")
    mix -= nshifts

    lpf=(fs=fs, from=fr, thru=th, boi=boi, gb=gb, mix=mix, flt=lpf)
    
    return (fs=fs, from=fr, thru=th, boi=boi, gb=gb, mix=mix, iqs=iqs), lpf
end

coef(f::FIRFilter)=f.h

# function out_of_band_suppression(in)
#     fs=sample_frequency(in)
#     boi=band_of_interest(in)
#     gb=guardband(in)
    
#     flt=remezfind(boi/fs, boi/fs+gb/fs) |> FIRFilter 

#     return out_of_band_suppression(in, (flt = flt,))
# end

function out_of_band_suppression(in, lpf)
    b= lpf.flt
    boi=band_of_interest(in)
    gb=guardband(in)
    fs=sample_frequency(in)
    iqs=inphase_n_quadratures(in)
    mix=mixer_frequency(data)
    
    d=length(coef(b))>>1

    fr=from(in)-d
    th=thru(in)-d

    u=iqs
    y=filt(b,u)
    boi= boi + 2gb
    
    out=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=y)
    
    o=length(coef(b))-1

    fr=th+1
    th=fr+o-1
    boi=band_of_interest(in)
    
    lpf=(fs=fs, from=fr, thru=th, boi=boi, gb=gb, mix=mix, flt=b)
    return out, lpf
end

function mixer(data)
    fs=sample_frequency(data)
    boi=band_of_interest(data)
    mix=mixer_frequency(data)
    
    fr= from(data)
    th= thru(data)
    
    iqs=inphase_n_quadratures(data)
    
    phase=oscillator(fs,mix,fr:th)
    
    iqs .*= phase
    mix=0

    return (fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs)
end

function upsample(data)
    fs=sample_frequency(data)
    boi=band_of_interest(data)
    mix=mixer_frequency(data)
    fr= from(data)
    th= thru(data)
    iqs= inphase_n_quadratures(data)

    n=2(th-fr+1)
    fs *= 2
    fr *= 2

    th = fr + n -1

    iqs=reshape(iqs,1,length(iqs))
    zs=zeros(eltype(iqs),size(iqs))
    iqs=vcat(iqs,zs)[:]
    

    return (fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs)
end

function suppress_mirror(in)
    fs= sample_frequency(in)
    boi= band_of_interest(in)
    
    Ws=Wp=boi/fs
    flt=find_halfband(Wp,Ws) |> FIRFilter 

    return suppress_mirror(in, (flt = flt,))
end

function suppress_mirror(in, hbf)
    b= hbf.flt
    boi=band_of_interest(in)

    fs=sample_frequency(in)
    iqs=inphase_n_quadratures(in)
    mix=mixer_frequency(data)
    
    d=length(coef(b))>>1

    fr=from(in)-d
    th=thru(in)-d

    u=iqs
    y=filt(b,u)
    
    out=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=y)

    o=length(coef(b))-1

    fr=th+1
    th=fr+o-1
    boi=band_of_interest(in)
    
    hbf=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, flt=b)
    return out, hbf
end

function flush(flt)
    b= flt.flt
    boi=band_of_interest(flt)
    gb=guardband(flt)
    fs=sample_frequency(flt)
    mix=mixer_frequency(flt)
    fr=from(flt)
    th=thru(flt)
    
    T=eltype(b.history)
    u=zeros(T,th-fr+1)
    y=filt(b,u)
    
    
    out=(fs=fs, from=fr, thru=th, boi=boi + 2gb, mix=mix, iqs=y)
    return out,nothing
end

function mix_n_boi(adata,bdata)
    fs= sample_frequency(adata)
    if fs!=sample_frequency(bdata)
        error("problex with fs.")
    end
    mixa=mixer_frequency(adata)
    mixb=mixer_frequency(bdata)
    boia=band_of_interest(adata)
    boib=band_of_interest(bdata)
    
    boi_x2=2abs(mixb-mixa)+boia+boib
    boi_x2+=isodd(boi_x2)
    boi=boi_x2 >> 1
    
    if mixa>mixb
        mix_x2=2mixa+boia+2mixb-boib
    else
        mix_x2=2mixb+boib+2mixa-boia
    end
    mix=mix_x2>>1
    return mix,boi        
end

function align_end(adata,bdata)
    th= min(thru(adata),thru(bdata))
    fs= sample_frequency(adata)

    x=nothing
    
    if thru(adata)>thru(bdata)
        iqs=inphase_n_quadratures(adata)
        fr=from(adata)
        mix=mixer_frequency(adata)
        boi=band_of_interest(adata)
        
        n=th-fr+1
        x=(fs=fs, from=th+1, thru=thru(adata), boi=boi, mix=mix, iqs=iqs[n+1:end])
        adata=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs[1:n])
        
    elseif thru(bdata)>thru(adata)
        iqs=inphase_n_quadratures(bdata)
        fr=from(bdata)
        mix=mixer_frequency(bdata)
        boi=band_of_interest(bdata)
        
        n=th-fr+1
        x=(fs=fs, from=th+1, thru=thru(bdata), boi=boi, mix=mix, iqs=iqs[n+1:end])
        bdata=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs[1:n])
        
    end

    return adata, bdata, x
end

function align_start(adata,bdata)
    fr= min(from(adata),from(bdata))
    fs= sample_frequency(adata)
    th= thru(adata)

    if th != thru(bdata)
        error("thru messed up")
    end
    
    boi=band_of_interest(adata)
    mix=mixer_frequency(adata)
    
    iqs=inphase_n_quadratures(adata)
    n=max(from(adata)-fr,0)
    zs=zeros(eltype(iqs),n)
    iqs=vcat(zs,iqs)
    
    adata=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs)

    boi=band_of_interest(bdata)
    mix=mixer_frequency(bdata)
    
    iqs=inphase_n_quadratures(bdata)
    n=max(from(bdata)-fr,0)
    zs=zeros(eltype(iqs),n)
    iqs=vcat(zs,iqs)
    
    bdata=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs)
    return adata, bdata
end

function mix_n_merge(adata,bdata)
    adata, bdata, x= align_end(adata, bdata)
    adata, bdata= align_start(adata, bdata)
    
    fs=sample_frequency(adata)
    fr= from(adata)
    th= thru(adata)
    mix,boi=mix_n_boi(adata,bdata)
    
    iqs=inphase_n_quadratures(adata)
    mixa=mixer_frequency(adata)
    phase=oscillator(fs, mixa-mix, fr:th)
    
    mixed_n_merged = iqs .* phase
    
    iqs=inphase_n_quadratures(bdata)
    mixb=mixer_frequency(bdata)
    phase=oscillator(fs, mixb-mix, fr:th)
    
    mixed_n_merged .+= iqs .* phase

    result=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=mixed_n_merged)

    return result, x
end

function output_buffer!(y,data)
    iqs= inphase_n_quadratures(data)
    
    fr=from(data)
    th=thru(data)
    indexes=collect(fr:th)
    inrange=findall(e -> e in 0:length(y)-1, indexes)
    y[indexes[inrange].+1]+=iqs[inrange]
    return y
end
