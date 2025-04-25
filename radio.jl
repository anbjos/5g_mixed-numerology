guardband(scs, bw; table=MIN_GUARDBAND_KHZ)=round(Int64, table[(scs=scs,bw=bw)]/7.5)
subcarrier_spacing(μ::Int)=15<<μ
hertz(n)=7500n
subcarrier_spacing_configuration(scs::Int)=findfirst(e -> e==scs, [15,30,60])-1

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

function radioDownLink(o::OranType3C; lowpassfilter=nothing, bw_table=BANDWIDTH_TABLE, sg_table=SIGNAL_GENERATION_TABLE)
    μ= subcarrier_spacing_configuration(o)
    nbins= number_of_bins(o)

    scs= subcarrier_spacing(o)
    nprbs= number_of_prbs(o)

    bw= bw_table[(scs=scs, nprbs=nprbs)]
    
    gb= guardband(scs,bw)
    boi= band_of_interest(o)

    fo= frequency_offset(o) << subcarrier_spacing_configuration(o)

    subFrameIdMod2 = o.subFrameId % 2
    slotId= o.slotId
    symbolId= o.startSymbolId
    extended= isextended(o)

    sgt= sg_table[(μ = μ, extended = extended, subFrameId = subFrameIdMod2, slotId = slotId, symbolId = symbolId)]

    fs_adjust= convert(Int64,log2(nbins/symbol_length(μ)))

    from= sgt.from << fs_adjust
    cp= sgt.cp  << fs_adjust
    symb = sgt.symbol  << fs_adjust

    iqs= copy(inphase_n_quadratures(o))

    lpf= lowpassfilter
    mix= nothing
    fs= nothing
    thru= nothing

    return μ, nbins, from, cp, symb, thru, boi, gb, fo, iqs, lpf, mix, fs
end


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



function iqmap(iqs)
    n= length(iqs)
    n2= n>>1
    m= nextpow(2,n)
    result=vcat(collect(m-n2+1:m),collect(1:n2))
    return result
end

@test iqmap(1:10)==[12,13,14,15,16,1,2,3,4,5]

function prbs2bins(rdl::RadioDownLink)
    μ= subcarrier_spacing_configuration(rdl)
    nbins= number_of_bins(rdl)
    fo = frequency_offset(rdl)
    
    iqs= inphase_n_quadratures(rdl)
    bins= zeros(eltype(iqs), nbins)
    bins[iqmap(iqs)] .= iqs
    
    mix = fo + length(iqs) << μ
    
    mixer_frequency!(rdl, mix)
    inphase_n_quadratures!(rdl, bins)

    return rdl
end

function  oscillator(fs, mix, range)
    ms= 0.001
    number_of_samples_in_2ms= hertz(fs)*2ms
    number_of_7500hz_periods_per_2ms= 15

    ω= number_of_7500hz_periods_per_2ms * mix
    result= exp.(1im*2π*(range)*ω/number_of_samples_in_2ms)
    return result
end

function phase_correction(rdl::RadioDownLink)
    fs= sample_frequency(rdl)
    mix= mixer_frequency(rdl)
    fr= from(rdl)
    cp= length_of_cyclic_prefix(rdl)
    
    phase= first(oscillator(fs, mix, fr+cp:fr+cp))
    
    iqs= inphase_n_quadratures(rdl) ./ phase
    inphase_n_quadratures!(rdl, iqs)

    return rdl
end

function create_lowpassfilter(rdl::RadioDownLink)
    lpf= lowpassfilter(rdl)
    fs = sample_frequency(rdl)
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

    fr=from(rdl)-delay
    th=thru(rdl)-delay
    from!(rdl,fr)
    thru!(rdl,th)

    iqs=inphase_n_quadratures(rdl)

    y= filt(lpf, iqs)
    lowpassfilter!(rdl, lpf)
    inphase_n_quadratures!(rdl, y)
    
    return rdl
end

# function filter_w_meta(rdl::RadioDownLink)
#     boi=band_of_interest(rdl)
#     gb=guardband(rdl)
#     fs=sample_frequency(rdl)
#     mix=mixer_frequency(rdl)
    
#     lpf= lowpassfilter(rdl)
#     delay=length(coef(lpf))>>1
#     fr=from(rdl)-delay
#     th=thru(rdl)-delay

#     fr=thru(rdl)+1
#     n=length(coef(lpf))-1
#     th=fr+n-1
    
#     return (fs=fs, from=fr, thru=th, boi=boi, gb=gb, mix=mix, flt=lpf)
# end

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
    
    phase=oscillator(fs,mix,fr:th)

    println(size(iqs),"x",size(phase)," $fr:$th")
    
    iqs .*= phase
    mix=0

    inphase_n_quadratures!(rdl,iqs)
    mixer_frequency!(rdl,mix)

    return rdl
end

function create_lowpassfilter(rdl::RadioDownLink)
    lpf= lowpassfilter(rdl)
    fs = sample_frequency(rdl)
    boi= band_of_interest(rdl)
    gb= guardband(rdl)

    lpf= isnothing(lpf) ? remezfind(boi/fs, boi/fs+gb/fs; Rs=db2amp(-26), Rp=db2amp(1.0)-1) |> FIRFilter : lpf

    lowpassfilter!(rdl, lpf)
    boi= boi + 2gb
    boi=band_of_interest!(rdl, boi)
    guardband!(rdl,nothing)

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



# function create_halfbandfilter(rdl::RadioDownLink)
#     fs= sample_frequency(rdl)
#     boi= band_of_interest(rdl)
    
#     Ws=Wp=boi/fs
#     hbf=find_halfband(Wp,Ws) |> FIRFilter 

#     return hbf
# end

# function flush(flt)
#     b= flt.flt
#     boi=flt.boi
#     fs=flt.fs
#     mix=flt.mix
#     fr=flt.from
#     th=flt.thru
    
#     T=eltype(b.history)
#     u=zeros(T,th-fr+1)
#     y=filt(b,u)
    
#     out=RadioDownLink()
#     sample_frequency!(out,fs)
#     from!(out,fr)
#     thru!(out,th)
#     band_of_interest!(out,boi)
#     mixer_frequency!(out,mix)
#     inphase_n_quadratures!(out,y)
    
#     return out, nothing
# end

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

function mix_n_boi(adata::RadioDownLink,bdata::RadioDownLink)
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

    x=RadioDownLink()
    
    if thru(adata)>thru(bdata)
        iqs=inphase_n_quadratures(adata)
        fr=from(adata)
        mix=mixer_frequency(adata)
        boi=band_of_interest(adata)
        
        n=th-fr+1
        sample_frequency!(x,fs)
        from!(x,th+1)
        band_of_interest!(x,boi)
        mixer_frequency!(x,mix)
        inphase_n_quadratures!(x,iqs[n+1:end])
        thru!(x,thru(adata))

        # x=(fs=fs, from=th+1, thru=thru(adata), boi=boi, mix=mix, iqs=iqs[n+1:end])
        #adata=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs[1:n])
        thru!(adata,th)
        inphase_n_quadratures!(adata,iqs[1:n])

        
    elseif thru(bdata)>thru(adata)
        iqs=inphase_n_quadratures(bdata)
        fr=from(bdata)
        mix=mixer_frequency(bdata)
        boi=band_of_interest(bdata)
        
        n=th-fr+1
        sample_frequency!(x,fs)
        from!(x,th+1)
        band_of_interest!(x,boi)
        mixer_frequency!(x,mix)
        inphase_n_quadratures!(x,iqs[n+1:end])
        thru!(x,thru(bdata))

        # x=(fs=fs, from=th+1, thru=thru(bdata), boi=boi, mix=mix, iqs=iqs[n+1:end])
#        bdata=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs[1:n])
        thru!(bdata,th)
        inphase_n_quadratures!(bdata,iqs[1:n])
    else
        x=nothing
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
    
    # boi=band_of_interest(adata)
    # mix=mixer_frequency(adata)
    
    iqs=inphase_n_quadratures(adata)
    n=max(from(adata)-fr,0)
    zs=zeros(eltype(iqs),n)
    iqs=vcat(zs,iqs)
    inphase_n_quadratures!(adata,iqs)
    from!(adata,fr)
    
    # adata=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs)

    # boi=band_of_interest(bdata)
    # mix=mixer_frequency(bdata)
    
    iqs=inphase_n_quadratures(bdata)
    n=max(from(bdata)-fr,0)
    zs=zeros(eltype(iqs),n)
    iqs=vcat(zs,iqs)
    inphase_n_quadratures!(bdata,iqs)
    from!(bdata,fr)
    
    # bdata=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=iqs)
    return adata, bdata
end

function mix_n_merge(adata,bdata)
    adata, bdata, x= align_end(adata, bdata)
    adata, bdata= align_start(adata, bdata)
    
    fs=sample_frequency(adata)
    fr= from(adata)
    th= thru(adata)
    mix,boi=mix_n_boi(adata,bdata)
    
    iqs= inphase_n_quadratures(adata)
    mixa= mixer_frequency(adata)
    phase= oscillator(fs, mixa-mix, fr:th)
    
    mixed_n_merged = iqs .* phase
    
    iqs= inphase_n_quadratures(bdata)
    mixb= mixer_frequency(bdata)
    phase= oscillator(fs, mixb-mix, fr:th)
    
    mixed_n_merged .+= iqs .* phase

    result=RadioDownLink()
    sample_frequency!(result,fs)
    from!(result,fr)
    thru!(result,th)
    band_of_interest!(result,boi)
    mixer_frequency!(result,mix)
    inphase_n_quadratures!(result,mixed_n_merged)

    # result=(fs=fs, from=fr, thru=th, boi=boi, mix=mix, iqs=mixed_n_merged)

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
