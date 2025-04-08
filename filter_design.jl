horner(poly, x) = foldr((left,right) -> left+right*x,poly; init=0) # poly[i] is the coefficient of x^(i-1)
poly(p,x)=horner(p[:],x)
poly(p,x::AbstractVector)=[poly(p,u) for u in x]

tf(b; z=exp.(im*π*(0:1/512:1)) )=poly(b,z) # [z^0,..z^n]


function findclusters(u::Vector,tol=1e-12)
    clusters=Any[]
    for e in u
        distances=abs.(e .- (clusters .|> c -> c.c))
    
        distance,p = isempty(distances) ? (NaN, 0) : findmin(distances)
        if distance<tol
            n=clusters[p].n
            c=clusters[p].c
            n+=1
            c=((n-1)*c+e)/n
            clusters[p]=(c=c,n=n)
        else
            cluster=(c=e, n=1)
            push!(clusters,cluster)
        end
    end
    return clusters
end

function with_multiplicity( r, mark, tol=1e-12 )
    clusters=findclusters(r, tol)
    for (center,n) in clusters
        if mark=="O"
            plot(real(center), imag(center), "bo")
            plot(real(center), imag(center), "w.")            
        else
            plot(real(center), imag(center), mark)
        end
        if n > 1
            annotate(string(n, "×"),
                     xy=(real(center), imag(center)),
                     xytext=(5,5), textcoords="offset points", color="b")
        end
    end
end

function resultT(b,r)
    c=isa(r[1],Complex)
    p=typeof(real(b[1]))
    
    lu=Dict((c=true, p=Float16) => ComplexF16, (c=true, p=Float32) => ComplexF32, (c=true, p=Float64) => ComplexF64, 
            (c=false, p=Float16) => ComplexF16, (c=false, p=Float32) => ComplexF32, (c=false, p=Float64) => ComplexF64,)
    return lu[(c=c, p=p)]    
end

function zplane(b,a=vcat(1,zeros(length(b)-1)))
    plt=figure(figsize=(6, 6))
    plt.gca().set_aspect("equal")
    θ = 0: 2π/2048: 2π
    circle=(x=cos.(θ),y=sin.(θ))
    plot(circle.x, circle.y;color="gray", linewidth=1, linestyle="--")
    z = roots(BigFloat.(b[end:-1:1]))
    p = roots(BigFloat.(a[end:-1:1]))
    
    T=resultT(b,z)

    z=convert.(T,z)
    p=convert.(T,p)

    with_multiplicity(z,"O")
    with_multiplicity(p,"x")

    xlabel("real")
    ylabel("imag")

    axvline(0, color="gray", linewidth=1, linestyle="--")
    axhline(0, color="gray", linewidth=1, linestyle="--")

    title("z-plane")

    return plt
end

function tf_on_arc(b,from,thru;grid=1024)
    z=exp.(im*π*(from:(thru-from)/grid:thru))
    h=poly(b,z)
    return h
end

function ripple(args...)
    f=abs.(tf_on_arc(args...))
    return max(abs(maximum(f)-1),abs(1-minimum(f)))  
end

function suppression(args...)
    f=abs.(tf_on_arc(args...))
    return maximum(f)
end

function _remezfind(Wp, Ws; Rp=db2amp(0.1)-1, Rs=db2amp(-20))
    b=vcat(zeros(100),1,zeros(100))
    return b
end

function remezfind(Wp, Ws; Rp=db2amp(0.1)-1, Rs=db2amp(-20))
    n=DSP.Filters.remezord(Wp,Ws,Rp,Rs)
    wp=1.0
    ws=1.0
    
    b=nothing
    e=(false,false)

    n= iseven(n) ? n+1 : n
    
    while n<512 && e != (true,true)
        bands=[(0, Wp)=>(1.0, wp), (Ws, 1.0)=>(0, ws)]
        b=DSP.Filters.remez(n,bands;Hz=2)

        
        e=(ripple(b, 0, Wp)<Rp, suppression(b,Ws,1)<Rs)

        if e==(true,false) && wp>1
            e=(false,false)
        elseif e==(true,false) && ws>1
            e=(false,false)
        end

        if e==(false,false) || ws>1000000 || wp>1000000
            n+=2
            wp=ws=1.0
        elseif e==(true,false)
            ws*=sqrt(2)
        elseif e==(false,true)
            wp*=sqrt(2)
        end
    end
    return b
end

function halfband(n, boi)
    o=n-1 # mult of two, but not 4.
    o+=[2 1 0 3][o%4+1]
    o2=o>>1

    bands=[(0, boi)=>0.5, (0.5, 0.5)=>0]
    p=DSP.Filters.remez(o2+1, bands; Hz=1)

    b=zeros(o+1)
    b[o2+1]=0.5
    b[1:2:end].=p
    return b
end

function find_halfband(Wp, Ws; Rp=db2amp(0.1)-1, Rs=db2amp(-60))    
    ok=false
    k=0
    b=nothing
    while !ok
        k+=1
        n=2+k*4+1
        b=halfband(n, Ws)
        ok= ripple(b, 0, Wp)<Rp
        ok &= suppression(b,1-Ws,1)<Rs
    end
    return b
end


