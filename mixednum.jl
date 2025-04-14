using Pkg
Pkg.activate(".")

using DSP
using PyPlot; pygui(true)
using PolynomialRoots
using FFTW
using Test
using Printf

include("./filter_design.jl")
include("./etsi.jl")
include("./o-ran.jl")
include("./radio.jl")
include("./user_equipment.jl")


function postprocess(iqs, oran, fs)
    prbs= oran |> oran2prbs
    r=fs ÷  sample_frequency(prbs)

    mx= -frequency_offset_7k5Hz(oran)
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
    # result=result[1:length(expected)]
        
    return result, expected
end



#iqs=o1.iqs
#iqs.=0
#iqs[end>>1+1]=1


scs=15
extended=0
bw=10
nrb=NRB[(scs=scs,bw=bw)]
fo=47
μ=subcarrier_spacing_configuration(scs)

symb=0
a1=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))
symb+=1
a2=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))

scs=30
extended=0
bw=10
nrb=NRB[(scs=scs,bw=bw)]
μ=subcarrier_spacing_configuration(scs)
fo=-24nrb-4

symb=0
b1=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))
symb+=1
b2=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))

y=zeros(ComplexF64,2(2048+512))

data, a_lpf=a1 |> oran2prbs |> prbs2bins |> phase_correction |> create_lowpassfilter |> amplitude_correction  |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier
data, a_lpf= out_of_band_suppression(data, a_lpf)
adata, a_hbf = suppress_mirror( data |> upsample )

# plot(pow2db.(abs2.(tf(coef(a_lpf.flt)))))

data, b_lpf=b1 |> oran2prbs |> prbs2bins |> phase_correction |> create_lowpassfilter |> amplitude_correction  |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier
data, b_lpf= out_of_band_suppression(data, b_lpf)
bdata, b_hbf = suppress_mirror( data |> upsample )

data, x= mix_n_merge(adata ,bdata)
data=mixer(data)
output_buffer!(y, data)

data, a_lpf=flush(a_lpf)
data=upsample(data)
adata, a_hbf = suppress_mirror(data, a_hbf)

data, b_lpf=flush(b_lpf)
data=upsample(data)
bdata, b_hbf = suppress_mirror(data, b_hbf)

data, x= mix_n_merge(bdata,x)
data=mixer(data)
output_buffer!(y, data)

data, x= mix_n_merge(adata,x)

data=mixer(data)
output_buffer!(y, data)


adata, a_hbf=flush(a_hbf)
bdata, b_hbf=flush(b_hbf)

data=mixer(adata)
output_buffer!(y, data)

data=mixer(bdata)
output_buffer!(y, data)

data=mixer(x)
output_buffer!(y, data) 


#figure(3)
#plot((0:length(y)-1)/length(y)*hertz(fs),10*log10.(abs2.(fft(hanning(length(y)) .* y))))

fs_out=sample_frequency(data)


result, expected= postprocess(y, a1, fs_out)
@test all(isapprox.(result[1:length(expected)],expected, atol=0.05))


y

plot(abs.(result[1:length(expected)].-expected),".")
plot(angle.(s[1:number_of_subcarriers(o)] ./ o.iqs))

#result=result[1:length(expected)]
figure(10)

figure()
plot(real(result),imag(result),".")

figure(2)
plot(abs.(result))
plot(abs.(expected),"--")

plot(abs.(result[1:length(expected)].-expected))




figure(2)
plot(abs.(result))
plot(abs.(expected),"--")

result=result[1:length(expected)]
figure(10)
plot(real(result),imag(result),".")



#################################################################

y=zeros(ComplexF64,2(2048+512))


scs=15
extended=0
bw=10
nrb=NRB[(scs=scs,bw=bw)]
fo=47
μ=subcarrier_spacing_configuration(scs)

symb=0
a1=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))
symb+=1
a2=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))

scs=30
extended=0
bw=10
nrb=NRB[(scs=scs,bw=bw)]
μ=subcarrier_spacing_configuration(scs)
fo=-24nrb-5

symb=0
b1=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))
symb+=1
b2=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))

data = a1 |> oran2prbs |> prbs2bins |> phase_correction |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier
data, a_lpf = out_of_band_suppression(data)
adata, a_hbf = suppress_mirror( data |> upsample )

data = b1 |> oran2prbs |> prbs2bins |> phase_correction |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier
data, b_lpf = out_of_band_suppression(data)
bdata, b_hbf = suppress_mirror( data |> upsample )

data, x= mix_n_merge(adata ,bdata)
data=mixer(data)
output_buffer!(y, data)

data, a_lpf=flush(a_lpf)
data=upsample(data)
adata, a_hbf = suppress_mirror(data, a_hbf)

data, b_lpf=flush(b_lpf)
data=upsample(data)
bdata, b_hbf = suppress_mirror(data, b_hbf)

data, x= mix_n_merge(bdata,x)
data=mixer(data)
output_buffer!(y, data)

data, x= mix_n_merge(adata,x)

data=mixer(data)
output_buffer!(y, data)


adata, a_hbf=flush(a_hbf)
bdata, b_hbf=flush(b_hbf)

data=mixer(adata)
output_buffer!(y, data)

data=mixer(bdata)
output_buffer!(y, data)

data=mixer(x)
output_buffer!(y, data)

o=a1

fs=sample_frequency(data)
s=symbol_with_cp(y,fs,o) |> odd_shift |> without_cyclicprefix  |> time2bins  |> shift_frequency_offset

expected = o |>  inphase_n_quadratures
result=s |>  inphase_n_quadratures
@test all(isapprox.(result[1:length(expected)],expected, atol=0.05))

figure(2)
plot(abs.(result))
plot(abs.(expected),"--")
