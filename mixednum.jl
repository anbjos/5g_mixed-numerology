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

data=a1 |>  RadioDownLink |> prbs2bins |> phase_correction |> create_lowpassfilter |> 
            amplitude_correction  |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier |> 
            out_of_band_suppression

a_lpf=filter_w_meta(data)
data=upsample(data)
a_hbf=create_halfbandfilter(data)
lowpassfilter!(data,a_hbf)
adata=suppress_mirror(data)
a_hbf=filter_w_meta(adata)

figure()
plot(pow2db.(abs2.(tf(coef(a_lpf.flt)))))

data=b1 |>  RadioDownLink |> prbs2bins |> phase_correction |> create_lowpassfilter |> 
            amplitude_correction  |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier |> 
            out_of_band_suppression

b_lpf=filter_w_meta(data)

data=upsample(data)
b_hbf=create_halfbandfilter(data)
lowpassfilter!(data,b_hbf)
bdata=suppress_mirror(data)
b_hbf=filter_w_meta(bdata)

data, x= mix_n_merge(adata ,bdata)

data=mixer(data)
output_buffer!(y, data)

data, a_lpf=flush(a_lpf)
data=upsample(data)
lowpassfilter!(data, a_hbf.flt)
adata= suppress_mirror(data)
a_hbf=filter_w_meta(adata)

data, b_lpf=flush(b_lpf)
data=upsample(data)
lowpassfilter!(data,b_hbf.flt)
bdata = suppress_mirror(data)
b_hbf=filter_w_meta(adata)

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

fs_out=sample_frequency(data)
 
result, expected= postprocess(y, a1, fs_out)
@test all(isapprox.(result,expected, atol=0.02))




#figure(3)
#plot((0:length(y)-1)/length(y)*hertz(fs),10*log10.(abs2.(fft(hanning(length(y)) .* y))))

figure()
plot(abs.(result.-expected),".")

figure()
plot(angle.(result ./ expected))

figure()
plot(real(result),imag(result),".")

figure()
plot(abs.(result))
plot(abs.(expected),"--")
