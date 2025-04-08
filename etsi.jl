
all_bit_vectors(n)=Iterators.product(Tuple([[false,true] for k in 1:n])...)

# 3GPP 38.211, 16QAM etc.
qam16(b)=((1-2b[1+0])*(2-(1-2b[1+2]))+im*(1-2b[1+1])*(2-(1-2b[1+3])))/sqrt(10)
qam64(b)=((1-2b[1+0])*(4-(1-2b[1+2])*(2-(1-2b[1+4])))+im*(1-2b[1+1])*(4-(1-2b[1+3])*(2-(1-2b[1+5]))))/sqrt(42)
qam256(b)=((1-2b[1+0])*(8-(1-2b[1+2])*(4-(1-2b[1+4])*(2-(1-2b[1+6]))))+1.0im*(1-2b[1+1])*(8-(1-2b[1+3])*(4-(1-2b[1+5])*(2-(1-2b[1+7])))))/sqrt(170)

const QAM16=[qam16(v) for v in all_bit_vectors(4)][:]
const QAM64=[qam64(v) for v in all_bit_vectors(6)][:]
const QAM256=[qam256(v) for v in all_bit_vectors(8)][:]

#ETSI 38104 section 5.3.2
const NRB=Dict( (scs=15, bw=3)=>15,  (scs=15, bw=5)=>25,   (scs=15, bw=10)=>52,  (scs=15, bw=15)=>79, 
                (scs=15, bw=20)=>106,(scs=15, bw=25)=>133, (scs=15, bw=30)=>160, (scs=15, bw=35)=>188,
                (scs=15, bw=40)=>216,(scs=15, bw=45)=>242, (scs=15, bw=50)=>270,
                (scs=30, bw=5)=>11,  (scs=30, bw=10)=>25,  (scs=30, bw=15)=>38,
                (scs=30, bw=20)=>51, (scs=30, bw=25)=>65,  (scs=30, bw=30)=>78,  (scs=30, bw=35)=>92,
                (scs=30, bw=40)=>106,(scs=30, bw=45)=>119, (scs=30, bw=50)=>133,
                (scs=30, bw=60)=>162,(scs=30, bw=70)=>189, (scs=30, bw=80)=>217, (scs=30, bw=90)=>245, (scs=30, bw=100)=>273,
                (scs=60, bw=10)=>11, (scs=60, bw=15)=>18,  (scs=60, bw=20)=>24,
                (scs=60, bw=25)=>31, (scs=60, bw=30)=>38,  (scs=60, bw=35)=>44,
                (scs=60, bw=40)=>51, (scs=60, bw=45)=>58,  (scs=60, bw=50)=>65,
                (scs=60, bw=60)=>79, (scs=60, bw=70)=>93,  (scs=60, bw=80)=>107, (scs=60, bw=90)=>121, (scs=60, bw=100)=>135)


 #ETSI 38104 section 5.3.3              
const MIN_GUARDBAND_KHZ=    Dict(
                                (scs=15, bw=3)=>142.5, (scs=15, bw=5)=>242.5, (scs=15, bw=10)=>312.5, (scs=15, bw=15)=>382.5, 
                                (scs=15, bw=20)=>452.5, (scs=15, bw=25)=>522.5, (scs=15, bw=30)=>592.5,(scs=15, bw=35)=>572.5,
                                (scs=15, bw=40)=>552.5,(scs=15, bw=45)=>712.5,(scs=15, bw=50)=>692.5,
                        
                                (scs=30, bw=5)=>505,(scs=30, bw=10)=>665,(scs=30, bw=15)=>645,(scs=30, bw=20)=>805,
                                (scs=30, bw=25)=>785,(scs=30, bw=30)=>945,(scs=30, bw=35)=>925,(scs=30, bw=40)=>905,
                                (scs=30, bw=45)=>1065,(scs=30, bw=50)=>1045,(scs=30, bw=60)=>825,(scs=30, bw=70)=>965,
                                (scs=30, bw=80)=>925,(scs=30, bw=90)=>885,(scs=30, bw=100)=>845,
                        
                                (scs=60, bw=10)=>1010,(scs=60, bw=15)=>990,(scs=60, bw=20)=>1330,(scs=60,bw=25)=>1310,
                                (scs=60, bw=30)=>1290,(scs=60, bw=35)=>1630,(scs=60, bw=40)=>1610,(scs=60, bw=45)=>1590,
                                (scs=60, bw=50)=>1570,(scs=60, bw=60)=>1530,(scs=60, bw=70)=>1490,(scs=60, bw=80)=>1450,
                                (scs=60, bw=90)=>1410,(scs=60, bw=100)=>1370
                            )


function bandwidth_table(nrb=NRB)
    result=Dict{NamedTuple{(:scs, :nprbs,), Tuple{Int64, Int64}}, Int64}()
    for k in keys(nrb)
        n=nrb[k]
        result[(scs=k.scs, nprbs=n)]=k.bw
    end
    return result
end

const BANDWIDTH_TABLE=bandwidth_table()

# Functions to provide range in terms of κ, see ETSI TS 138 211, '5.3.1 OFDM baseband signal generation for all channels except PRACH'
# Covers two subFrames

cyclic_prefix_length(μ,l;extended=false)= extended ? 512 >> μ : 144 >> μ + (l==0 || l==7<<μ) * 16
symbol_length(μ)=2048 >> μ

function  signal_generation_table()
    slot_range=Dict((μ=0) => 0:0, (μ=1) => 0:2:2, (μ=2) => 0:3)
    symbol_range=Dict(false => 0:13, true => 0:11)

    result=Dict{NamedTuple{(:μ, :extended, :subFrameId, :slotId, :symbolId), Tuple{Int64, Bool, Int64, Int64, Int64}}, NamedTuple{(:from, :cp, :symbol), Tuple{Int64, Int64, Int64}}}()

    for (μ, extended) in [(0,false),(1,false),(2,false),(2,true)]
        count=0
        for subFrameId in 0:1
            for slotId in slot_range[μ]
                for symbolId in symbol_range[extended]
                    cp=cyclic_prefix_length(μ,symbolId;extended=extended)
                    s=symbol_length(μ)
                    result[(μ=μ, extended=extended, subFrameId=subFrameId, slotId=slotId, symbolId=symbolId)]=(from=count, cp=cp, symbol=s)
                    count += cp+s
                end 
            end 
        end
    end
    return result
end

const SIGNAL_GENERATION_TABLE=signal_generation_table()

