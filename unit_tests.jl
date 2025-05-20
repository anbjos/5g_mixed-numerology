# This code is a simplified processing of the orans without scheduling
# Usefull for generating unit tests

y=zeros(ComplexF64,2(2048+512))

orans=[a1,a2,b1,b2]
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
 
result, expected= get_symbol(y, a1, fs_out)
@test all(isapprox.(result,expected, atol=0.03))
