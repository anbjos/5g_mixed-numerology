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

frequency_offset(o::OranType3C)=o.freqOffset

number_of_prbs(o::OranType3C)=o.numPrbs

number_of_subcarriers(o::OranType3C)=12number_of_prbs(o)

subcarrier_spacing(o::OranType3C)=15 << subcarrier_spacing_configuration(o)

subcarrier_spacing_configuration(o::OranType3C)= o.frameStructure.μ

cp_length(o)=o.cpLength

sample_frequency(o::OranType3C)= (1000subcarrier_spacing(o) * number_of_bins(o)) ÷ 7500

isextended(o::OranType3C)= (cp_length(o) in 512 .>> (0:4))

frequency_offset_7k5Hz(o::OranType3C)=frequency_offset(o) << subcarrier_spacing_configuration(o)

inphase_n_quadratures(o::OranType3C)=o.iqs

band_of_interest(o::OranType3C)=2number_of_subcarriers(o)<< subcarrier_spacing_configuration(o)

bandwidth(o::OranType3C; table=BANDWIDTH_TABLE)=table[(scs=subcarrier_spacing(o),nprbs=number_of_prbs(o))]

guardband(o::OranType3C)=guardband(subcarrier_spacing(o),bandwidth(o))

