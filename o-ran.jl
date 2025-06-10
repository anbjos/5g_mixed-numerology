mutable struct OranType3C
    frameId
    subFrameId
    slotId
    startSymbolId
    frameStructure
    cpLength
    startPrbc
    numPrbs
    freqOffset
    iqs
end

Base.show(io::IO,o::OranType3C)=print(io,"fr=$(o.frameId), sf=$(o.subFrameId),sl=$(o.slotId),sy=$(o.startSymbolId),fs=$(o.frameStructure),st=$(o.startPrbc),n=$(o.numPrbs),fo=$(o.freqOffset)")

function iq_values(modulation, scs, bw, nrb=NRB)
    n=12nrb[( scs=scs, bw=bw )]
    r=length(modulation)
    result=[modulation[rand(1:r)] for k in 0:n-1]
    #result=[modulation[mod(k,r)+1] for k in 0:n-1]
    return result
end

function resolution_for_cp(n,μ)
    if μ==2
        result=max(256,n)
    else
        result=max(128,n)
    end
    return result
end

function number_of_bins(o::OranType3C)
    n_sc=number_of_subcarriers(o)
    μ=subcarrier_spacing_configuration(o)
    np=nextpow(2,n_sc)
    result=resolution_for_cp(np, μ)
    return result
end

frequency_offset_in_half_subcarriers(o::OranType3C)=o.freqOffset
frequency_offset(o::OranType3C)=frequency_offset_in_half_subcarriers(o)
frequency_offset_in_units_of_7k5Hz(o::OranType3C)=frequency_offset_in_half_subcarriers(o)  << subcarrier_spacing_configuration(o) 

number_of_prbs(o::OranType3C)=o.numPrbs
number_of_subcarriers(o::OranType3C)=12number_of_prbs(o)

subcarrier_spacing(o::OranType3C)=15 << subcarrier_spacing_configuration(o)
subcarrier_spacing_configuration(o::OranType3C)= o.frameStructure.μ

cp_length(o)=o.cpLength

sample_frequency(o::OranType3C)= (1000subcarrier_spacing(o) * number_of_bins(o)) ÷ 7500
isextended(o::OranType3C)= (cp_length(o) in 512 .>> (0:4))

inphase_n_quadratures(o::OranType3C)=o.iqs

band_of_interest_in_half_subcarriers(o::OranType3C)=2number_of_subcarriers(o)
band_of_interest_in_units_of_7k5Hz(o::OranType3C)=band_of_interest_in_half_subcarriers(o) << subcarrier_spacing_configuration(o) 

bandwidth(o::OranType3C; table=BANDWIDTH_TABLE)=table[(scs=subcarrier_spacing(o), nprbs=number_of_prbs(o))]
guardband(o::OranType3C)=guardband(subcarrier_spacing(o), bandwidth(o))

function align_with_nummerology(μ, n)
    p=n
    while p % (2^μ) !=0
        p+=1
    end
    return p
end

function total_bandwidth(allocations, NRB=NRB,MIN_GUARDBAND_KHZ=MIN_GUARDBAND_KHZ)
    n7k5=0
    old_gb=0
    for k in allocations
        scs=k.scs
        bw=k.bw
        μ=subcarrier_spacing_configuration(scs)
        n7k5+=12*NRB[(scs=scs, bw=bw)]<<(μ+1)
        new_gb=round(Int64,2MIN_GUARDBAND_KHZ[(scs=scs, bw=bw)] ÷ 15)
        n7k5+=max(new_gb,old_gb)
        old_gb=new_gb
    end
    return n7k5
end


function allocate_spectrum(bandwidths, fs_out, NRB=NRB, MIN_GUARDBAND_KHZ=MIN_GUARDBAND_KHZ)
    result=[]

    tbw=total_bandwidth(bandwidths)
    if tbw > fs_out
        error("illegal frequency allocation, bw: $tbw, fs: $fs_out")
    end

    freq_offset_7k5=-(tbw >> 1)

    old_gb=0
    for (k, nummerology) in enumerate(bandwidths)
        scs=nummerology.scs
        bw=nummerology.bw
        μ=subcarrier_spacing_configuration(scs)

        new_gb=round(Int64,2MIN_GUARDBAND_KHZ[(scs=scs, bw=bw)] ÷ 15)
        freq_offset_7k5 += k>1 ? max(old_gb, new_gb) : 0
        old_gb = new_gb

        freq_offset_7k5=align_with_nummerology(μ, freq_offset_7k5)
        fo = freq_offset_7k5 >> μ
        push!(result,(scs=scs, bw=bw, fo=fo))
        freq_offset_7k5 += 12*NRB[(scs=scs, bw=bw)]<<(μ+1)
    end

    return result
end
