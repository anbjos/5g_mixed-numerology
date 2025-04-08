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

####  complete chain:




#iqs=o1.iqs
#iqs.=0
#iqs[end>>1+1]=1
#iqs=o2.iqs
#iqs.=0
#iqs[end>>1+1]=1
#iqs[1]=1


function shift_to_dc(y)
    oran=y.oran
    offset=number_of_subcarriers(oran)>>1
    println("offset $offset")

    iqs=inphase_n_quadratures(y)

    iqs=circshift(iqs, offset)
    return (oran=oran, fs=fs, iqs=iqs)
end



y=zeros(ComplexF64,2(2048+512))

scs=15
extended=0
bw=10
nrb=NRB[(scs=scs,bw=bw)]
fo=40
μ=subcarrier_spacing_configuration(scs)
symb=0

a1=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))



#iqs=a1.iqs
#iqs.=0
#iqs[end>>1+1]=1
#iqs[1]=1


data, a_lpf=a1 |> oran2prbs |> prbs2bins |> phase_correction |> create_lowpassfilter |> amplitude_correction  |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier
data, a_lpf= out_of_band_suppression(data, a_lpf)
data=mixer(data)
output_buffer!(y, data)

data, a_lpf=flush(a_lpf)
data=mixer(data)
output_buffer!(y, data)

o=a1
fs=sample_frequency(data)
s=symbol_with_cp(y,fs,o) |> odd_shift |> without_cyclicprefix   |> time2bins |> shift_frequency_offset
expected = o |>  inphase_n_quadratures
result=s |>  inphase_n_quadratures

@test all(isapprox.(result[1:length(expected)],expected, atol=0.01))

plot(abs.(result[1:length(expected)].-expected))




figure(2)
plot(abs.(result))
plot(abs.(expected),"--")

result=result[1:length(expected)]
figure(10)
plot(real(result),imag(result),".")


#mixer_frequency(o::OranType3C)=frequency_offset_7k5Hz(o) + number_of_subcarriers(o)  << subcarrier_spacing_configuration(o)
#center_offset(o::OranType3C)= 1 << subcarrier_spacing_configuration(o)



o=a1

#move to dc

fs=sample_frequency(data)
prbs= o |> oran2prbs
r=fs ÷  sample_frequency(prbs)

mx=-frequency_offset_7k5Hz(o)
ω=oscillator(fs,mx,0:length(y)-1)
yy=ω .* y

#figure(3)
#plot((0:length(y)-1)/length(y)*hertz(fs),10*log10.(abs2.(fft(hanning(length(y)) .* y))))
#plot((0:length(y)-1)/length(y)*hertz(fs),10*log10.(abs2.(fft(hanning(length(y)) .* yy))))

cp=cyclic_prefix(prbs)
fr=r*from(prbs)+r*cp
ϕ=first(oscillator(fs,mx,fr:fr))
yy=yy ./ ϕ


# mx=-number_of_subcarriers(o) << subcarrier_spacing_configuration(o)
# ω=oscillator(fs,mx,1:length(y))
# yy=ω .* yy

# plot((0:length(y)-1)/length(y)*hertz(fs),10*log10.(abs2.(fft(hanning(length(y)) .* yy))))

# boi=band_of_interest(o)
# gb=guardband(o)

# b=remezfind(boi/fs, boi/fs+gb/fs; Rs=db2amp(-20)) 
# plot(10log10.(abs2.(tf(b))))
# d=length(b)>>1
# yy=filt(b,vcat(yy,zeros(d)))[1+d:end]


# figure(3)
# plot((0:length(y)-1)/length(y)*hertz(fs),10*log10.(abs2.(fft(hanning(length(y)) .* yy))))



# ω=oscillator(fs,-mx,1:length(y))
# yy=ω .* yy

cp2=cyclic_prefix(prbs)>>1
fr=r*from(prbs)+r*cp2
th=r*from(prbs)+r*cp2+r*symbol(prbs)-1

yy=yy[fr+1:th+1]
#yy=vcat(yy[1:end-cp],yy[end-cp+1:end])

cp2
n=length(yy)
delay=zeros(n)
delay[cp2+1]=1
ω=fft(delay)

s=fft(yy) ./ ω
plot(abs.(s))
plot(abs.(o.iqs))

plot(real.(s))
plot(real.(o.iqs))

result=s
expected=o.iqs

plot(abs.(result[1:length(expected)].-expected),".")

@test all(isapprox.(result[1:length(expected)],expected, atol=0.01))

plot(angle.(s[1:number_of_subcarriers(o)] ./ o.iqs))

##############################################################



plot(real.(s))
plot(real.(o.iqs))


s=symbol_with_cp(yy,fs,o) |> odd_shift |> without_cyclicprefix   |> time2bins

nbins=length(s.iqs)
iqs= inphase_n_quadratures(s)

z=exp.(1im *2π * (0.5:nbins-0.5) ./ nbins)

ampl=abs.(tf(b;z=z))
ampl[ampl .< 1e-12] .= 1

iqs ./= ampl

plot(ampl)

number_of_subcarriers(o)



expected = o |>  inphase_n_quadratures
result=s |>  inphase_n_quadratures
@test all(isapprox.(result[1:length(expected)],expected, atol=0.02))

figure()
plot(abs.(result))

# s=symbol_with_cp(yy,fs,o) |> odd_shift |> without_cyclicprefix   |> time2bins

figure(22)
plot(abs.(result))
plot(abs.(expected),"--")





frequency_offset(s.oran)


-(frequency_offset(s.oran) >> 1)





#data = a1 |> oran2prbs |> prbs2bins |> phase_correction |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier
#data, a_lpf = out_of_band_suppression(data)
data, a_hbf = suppress_mirror( data |> upsample )
data=mixer(data)
output_buffer!(y, data)

data, a_lpf=flush(a_lpf)
data, a_hbf = suppress_mirror( data |> upsample, a_hbf)
data=mixer(data)
output_buffer!(y, data)

data, a_hbf=flush(a_hbf)
data=mixer(data)
output_buffer!(y, data)

scs=30
extended=0
bw=10
nrb=NRB[(scs=scs,bw=bw)]
μ=subcarrier_spacing_configuration(scs)
fo=-24nrb-12

symb=0
b1=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))

data = b1 |> oran2prbs |> prbs2bins |> phase_correction |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier
data, b_lpf = out_of_band_suppression(data)
data, b_hbf = suppress_mirror( data |> upsample )

data=mixer(data)
output_buffer!(y, data)

data, b_lpf=flush(b_lpf)
data, b_hbf = suppress_mirror( data |> upsample, b_hbf)
data=mixer(data)
output_buffer!(y, data)

data, b_hbf=flush(b_hbf)
data=mixer(data)
output_buffer!(y, data)

fs=sample_frequency(data)
figure(1)
plot((0:length(y)-1)/length(y)*hertz(fs),10*log10.(abs2.(fft(hanning(length(y)) .* y))))
grid(true)

o=a1
s=symbol_with_cp(y,fs,o) |> odd_shift |> without_cyclicprefix   |> time2bins |> shift_frequency_offset
expected = o |>  inphase_n_quadratures
result=s |>  inphase_n_quadratures
@test all(isapprox.(result[1:length(expected)],expected, atol=0.02))

figure(2)
plot(abs.(result))
plot(abs.(expected),"--")

result=result[1:length(expected)]
plot(real(result),imag(result),".")


plot(pow2db.(abs2.(tf(coef(a_hbf.flt)))))

figure()
plot(pow2db.(abs2.(tf(coef(a_lpf.flt)))))

#plot(imag(y[1:1:end]))
mixer_frequency(o::OranType3C)=frequency_offset_7k5Hz(o) + number_of_subcarriers(o)  << subcarrier_spacing_configuration(o)

1-db2amp(0.01)

1-db2amp(0.1)

o=a1
fs=sample_frequency(data)

mx=mixer_frequency(o)
ω=oscillator(fs,-mx,1:length(y))

yy=ω .* y

figure(3)
plot((0:length(y)-1)/length(y)*hertz(fs),10*log10.(abs2.(fft(hanning(length(y)) .* yy))))



boi=band_of_interest(o)
gb=guardband(o)

flt=remezfind(boi/fs, boi/fs+gb/fs; Rs=db2amp(-20)) |> FIRFilter

d=length(coef(flt))>>1

yy=filt(flt,yy)
yy=vcat(yy[d+1:end],yy[1:d])

ω=oscillator(fs,mx,1:length(y))
yy=yy .* ω

s=symbol_with_cp(yy,fs,o) |> odd_shift |> without_cyclicprefix   |> time2bins |> shift_frequency_offset
expected = o |>  inphase_n_quadratures
result=s |>  inphase_n_quadratures
@test all(isapprox.(result[1:length(expected)],expected, atol=0.02))

figure()
plot(isapprox.(result[1:length(expected)],expected, atol=0.02))



figure(4)
plot(abs.(result))
plot(abs.(expected),"--")



#mix=mixer_frequency(data)

#number_of_subcarriers(o)
#frequency_offset(o)
#mixer_frequency(o)
#mixer_frequency(o |> oran2prbs |> prbs2bins)
#subcarrier_spacing(b1)
#typeof(b1)
#oran2prbs(b1)
#subcarrier_spacing(b1)
#number_of_subcarriers(b1)

s=symbol_with_cp(y,fs,o) |> odd_shift |> without_cyclicprefix   |> time2bins |> shift_frequency_offset
expected = o |>  inphase_n_quadratures
result=s |>  inphase_n_quadratures
@test all(isapprox.(result[1:length(expected)],expected, atol=0.02))

figure(3)
plot(abs.(result))
plot(abs.(expected),"--")


#plot(isapprox.(result[1:length(expected)],expected, atol=0.02))


plot(imag(s.iqs))

figure()
plot(abs.(s.iqs))

figure(4)
result=result[1:length(expected)]
plot(real(result),imag(result),".")

figure(15)
plot(real(result),".")

#25 tops, +/- 0.04

figure(6)
plot(real(result),".")
grid(true)


ref=copy(y)



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
