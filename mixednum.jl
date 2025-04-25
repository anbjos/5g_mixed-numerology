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
a1=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM256,scs,bw))
symb+=1
a2=OranType3C(0,0,0,symb,(μ=μ,),extended,0,nrb,fo,iq_values(QAM256,scs,bw))

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

freq_n_from(o::OranType3C)= o |> RadioDownLink |> prbs2bins |> r -> (frequency_offset(r), from(r))



from_div_fs(u)= u |> RadioDownLink |> r -> from(r)/sample_frequency(r)

Base.show(io::IO,rdl::RadioDownLink)=print(io,"mix=$(mixer_frequency(rdl)) @ fs=$(sample_frequency(rdl)), $(from(rdl)):$(thru(rdl))")

static_chain(rdl::RadioDownLink)= rdl |>    prbs2bins |> phase_correction |> create_lowpassfilter |> 
                                            amplitude_correction  |> bins2symbol |> with_cyclic_prefix |> shift_half_subcarrier |> 
                                            out_of_band_suppression

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
                                            
                                            
############################################################################################

y=zeros(ComplexF64,2(2048+512))

orans=[a1,a2,b1,b2]
orans=[a1]
datas,flts=process_orans(orans)

for data in datas
    print("$data -> ")
    data=upsample(data)
    fr=from(data)
    mx=mixer_frequency(data)
    p=findfirst(e -> mixer_frequency(e) == mx && from(e)==fr, flts)
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
    data=mixer(data)
    println("out($data)")
    output_buffer!(y, data)
end

for p in 1:length(flts)
    flt=flts[1]
    deleteat!(flts,1) 
    data=flush(flt)
    print("$data -> ")
    data=mixer(data)
    println("out($data)")
    output_buffer!(y, data)
end

fs_out=4096
 
result, expected= postprocess(y, a1, fs_out)
@test all(isapprox.(result,expected, atol=0.03))


###########################################################################################
function can_merge(a,b)
    isnothing(a) && return false
    isnothing(b) && return false
    ok=(sample_frequency(a)==4096)
    ok &=!(from(a)>thru(b) || from(b)>thru(a))
    ok &= mixer_frequency(a) in [-610,670]
    ok &= mixer_frequency(b) in [-610,670]
    return ok    
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



y=zeros(ComplexF64,4(2048+512))

orans=[a1,a2,b1,b2]
datas,flts=process_orans(orans)

fs_in=512
fs_in=2048
fs_out=4096

fs=fs_in
t0=0
t=t0
 
sample_frequencies(fs_in,fs_out)=2 .^ (Int(log2(fs_in)):Int(log2(fs_out)))

issomething(u) = ! isnothing(u)

while !isempty(datas)
    println("t:$t")
    for fs in sample_frequencies(fs_in,fs_out)
        while (a=find_data!(datas,fs,t)) |> issomething
        
            b=find_next!(datas,fs,t,a)
        
            if can_merge(a,b)
                # deleteat!(datas,p)
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
        # p=find_data(datas,fs,t)
    end
    
    t0+=1024
    t=t0
    fs=fs_in
    # p=find_data(datas,fs,t)
end


result, expected= postprocess(y, a2, fs_out)
@test all(isapprox.(result,expected, atol=0.03))

figure()
plot(abs.(result))
plot(abs.(expected),"--")






a


fs
sample_frequency.(datas)
from.(datas)























while !isnothing(p0) && !isnothing(p1)
    a=datas[p0]
    b=datas[p1]
    
    if can_merge(a,b)
        deleteat!(datas,max(p1,p0))
        deleteat!(datas,min(p1,p0))
        a,b=mix_n_merge(a,b)
        push!(datas,a)
        push!(datas,b)
        datas=datas[sortperm(mixer_frequency.(datas))]
        p0=find_data(datas,fs,t,1)
        p1=find_data(datas,fs,t,p0+1)
        println("*")
    else
        p0=p1
        p1=find_data(datas,fs,t,p0+1)
        println("-")
    end
end

#interpolation
p=findfirst(d -> (sample_frequency(d)==fs && from(d) <=t), datas)
while !isnothing(p)
    data=datas[p]
    deleteat!(datas,p)
    data=upsample(data)
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
    p=findfirst(d -> (sample_frequency(d)==fs && from(d) <=t), datas)
end

for p in length(flts):-1:1
    flt=flts[p]
    if from(flt)<t
        deleteat!(flts,p) 
        data=flush(flt)
        push!(datas,data)
        print("-")
    else
        print("o")
    end
end




[(sample_frequency(d), from(d), thru(d)) for d in datas]

from.(datas)

plot(abs.(y))
plot(abs.(ref),"--")


y=zeros(ComplexF64,2(2048+512))



ref=copy(y)







d1=nothing
for data in datas[1:1]
    data=upsample(data)
    fr=from(data)
    mx=mixer_frequency(data)
    p=findfirst(e -> mixer_frequency(e) == mx && from(e)==fr, flts)
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
    d1=data
    # data=mixer(data)
    # output_buffer!(y, data)
end

d2=nothing
for data in datas[2:2]
    data=upsample(data)
    fr=from(data)
    mx=mixer_frequency(data)
    p=findfirst(e -> mixer_frequency(e) == mx && from(e)==fr, flts)
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
    d2=data
    # data=mixer(data)
    # output_buffer!(y, data)
end

d2



d,x=mix_n_merge(d1,d2)
d=mixer(d)
output_buffer!(y,d)
x=mixer(x)
output_buffer!(y,x)



for data in datas[3:end]
    data=upsample(data)
    fr=from(data)
    mx=mixer_frequency(data)
    p=findfirst(e -> mixer_frequency(e) == mx && from(e)==fr, flts)
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
    data=mixer(data)
    output_buffer!(y, data)
end



flts

data=flush(flts[1])
data=mixer(data)
output_buffer!(y, data)

data=flush(flts[2])
data=mixer(data)
output_buffer!(y, data)




fs_out=sample_frequency(data)
 
result, expected= postprocess(y, a1, fs_out)
@test all(isapprox.(result,expected, atol=0.03))

figure()
plot(real(result),imag(result),".")




data=datas[1]
data=upsample(data)
data=create_halfbandfilter(data) |> suppress_mirror
flt=filter_w_meta!(data)
push!(flts,flt)
data=mixer(data)
output_buffer!(y, data)


data=datas[2]






halfbandfilter!(data,halfbandfilter(flt))
date=suppress_mirror(data)
flt=filter_w_meta!(data)
data=mixer(data)
output_buffer!(y, data)








for data in datas
    data=mixer(data)
    output_buffer!(y, data)
end




t=0
p=from.(datas).<t

vdatas=datas[p]

frequency_offset.(vdatas)

vdatas=vdatas[sortperm(frequency_offset.(vdatas))]


v1=vdatas[1]
v2=vdatas[2]

v1=upsample(v1)
v2=upsample(v2)

v1=create_halfbandfilter(v1) |> suppress_mirror
flt=filter_w_meta!(v1)
push!(flts,flt)

v2=create_halfbandfilter(v2) |> suppress_mirror
flt=filter_w_meta!(v2)
push!(flts,flt)



data, x= mix_n_merge(v1 ,v2)
data=mixer(data)
output_buffer!(y, data)




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
