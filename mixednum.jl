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
include("./get_symbol.jl")

fs_in=512
fs_out=4096

bandwidths=[[(scs=15, bw=10),(scs=30,bw=10),(scs=60,bw=10),],[(scs=15,bw=5),(scs=30,bw=5),(scs=60,bw=20),]]
orans=[]
slot_range(scs)=[0:0,0:2:2,0:2][subcarrier_spacing_configuration(scs)+1]

for sub_frame in 0:length(bandwidths)-1
    allocations=allocate_spectrum(bandwidths[sub_frame+1], fs_out)
    for allocation in allocations
        scs=allocation.scs
        bw=allocation.bw
        fo=allocation.fo
        μ=subcarrier_spacing_configuration(scs)
        extended=0
        nrb=NRB[(scs=scs,bw=bw)]
        for slot in slot_range(scs)
            for symb in 0:13
                e=OranType3C(0,sub_frame,slot,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM64,scs,bw))
                push!(orans,e)
            end
        end
    end
end

y=zeros(ComplexF64,hertz(fs_out) ÷ 500) # 2 ms
datas,flts=process_orans(orans)
process_data!(y,datas, flts, fs_in, fs_out)


result,expected=nothing,nothing
for o in orans
    result, expected= get_symbol(y, o, fs_out)
    if !all(isapprox.(result,expected, atol=0.05))
        figure()
        plot(real(result),imag(result),".")
        println("$(o.frameId) $(o.subFrameId) $(o.slotId) $(o.startSymbolId) $(o.frameStructure) $(o.cpLength) st: $(o.startPrbc) n:$(o.numPrbs) fo:$(o.freqOffset)")
        error("x")
    end
    @test all(isapprox.(result,expected, atol=0.05))
end


s=welch_pgram(y[1:end>>1], 4096; fs=hertz(fs_out), window=hanning)
figure()
plot(fftshift(freq(s)).+hertz(fs_out>>1),pow2db.(fftshift(power(s))))
title("Power spectrum, 15, 30 and 60 kHz nummerologies\n bw=3 × 10 MHz")
xlabel("Hz")
ylabel("dB")
grid(true)
#savefig("./Figures/spectrum.png")

o=orans[1]

flt=static_chain(RadioDownLink(o)).flt
b=coef(flt)
figure()
plot(0:1/512:1,amp2db.(abs.(tf(b))))
title("out of band suppression filter\nscs=15 kHz, bw=10MHz")
ylabel("dB")
xlabel("freq. (rel. to Nyquist/2)")
grid(true)
#savefig("./Figures/oob_tf.png")

zplane(b)
gcf()
title("z-plane for out of band suppression filter")
#savefig("./Figures/oob_zp.png")

d=static_chain(RadioDownLink(o))
d2=upsample(d)
halfbandfilter!(d2, nothing)
flt=halfbandfilter(create_halfbandfilter(d2))
b=coef(flt)
figure()
plot(0:1/512:1,amp2db.(abs.(tf(b))))
title("halfband filter\nscs=15 kHz, bw=10MHz")
ylabel("dB")
xlabel("freq. (rel. to Nyquist/2)")
grid(true)
#savefig("./Figures/hb_tf.png")

zplane(b)
gcf()
title("z-plane for halfband filter")
#savefig("./Figures/hb_zp.png")

result, expected= get_symbol(y, o, fs_out)

figure()
plot(real(result),imag(result),".")
title("constalation diagram\nscs=15 kHz, bw=10MHz")
xlabel("real")
ylabel("imag")
#savefig("./Figures/constalation.png")

figure()
plot(0:length(result)-1,abs.(result.-expected))
title("absolute error\nscs=15 kHz, bw=10MHz")
xlabel("error")
ylabel("subcarrier")
grid(true)
#savefig("./Figures/error.png")
