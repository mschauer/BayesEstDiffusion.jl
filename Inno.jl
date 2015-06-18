module Inno
export iter, Atan, OU, generate, runexamples
using Winston
#using Debug
#using Distributions
include("SDE.jl")
using .SDE
include("misc.jl")
import .SDE.sigma, .SDE.b, .SDE.a #, SDE.UvPro, SDE.CTPro
import .SDE.dotVs, .SDE.Vs



####### dimension

const d = 1

######## data?
generated = false


####### Atan
  
type Atan <: UvPro
    thC
    thnonC
    d::Int   
    Atan(B, beta, sigma) = new([B, beta], [sigma], 1)
end
nphi(_::Atan) = 2
nth(_::Atan) = 3


b(t,x, P::Atan) =  P.thC[1]*phi1(t,x, P) + P.thC[2]*phi2(t, x, P)
sigma(t,x, P::Atan) = P.thnonC[1]

Bof(P::Atan) = P.thC[1]*(cos(-P.thC[2]/P.thC[1]))^2
betaof(P::Atan) = 1/(2P.thC[1])*sin(2*P.thC[2]/P.thC[1])

#Bof(P::Atan) = P.thC[1]
#betaof(P::Atan) = P.thC[2]

muof(P::Atan) = 0.

phi1(t, x, _::Atan) = atan(x)
phi2(t, x, _::Atan) = 1.
phi(t, x, _::Atan) = [atan(x), 1.]





####### Prior
import Base.clamp;
clamp(x::Vector,l::Vector, u::Vector) = map(i -> clamp(x[i], l[i], u[i]), 1:length(x))



function generate(Model_; thtrue = [-2., 0.0, 0.75], u = 0.1, Nobs = 101, N = 400000 + 1, Delta=0.3, SV=true)
    global Model = Model_
    global Ptrue = Model(thtrue...)
    global model = string(Model)

    global M = Nobs - 1 # number of bridges
    # N - full observations i*M*n+1 for natural number i
    global Ttotal = Delta*M      #total time span

    tag = "$model$(M)T$Ttotal"
    println(tag)

    ############

    ############
    print("Generate X")
    ############

    srand(3)
    tt = collect(linspace(0., Ttotal, N))
    W = sample(tt, Wiener())
    xx = euler(u, W, Ptrue).yy
    println(".")
    
    global ttf = tt[1:(N-1)//(M*100):end]
    global xxf = xx[1:(N-1)//(M*100):end]
    global ttd = tt[1:(N-1)//M:end]
    global xxd = xx[1:(N-1)//M:end]

    global quxx = (sum(diff(xx[:]).^2)/Ttotal)
    global quxxd = (sum(diff(xxd[:]).^2)/Ttotal)
    println("si2 $( (thtrue[end])^2) qufull $quxx quobs $quxxd")
    global generated = true
    println("beta ", (xxd[end] - xxd[1])/Ttotal)
    
    
    if SV
        writecsv("xtrue$(tag).csv",[linspace(0, Ttotal, length(xx)) xx])
        writecsv("xobs$(tag).csv",[linspace(0, Ttotal, length(xxd)) xxd])
    end
end


#randq(thold, dl) = clamp( thold .* [ones(Nphi), exp(dl*(2rand(Nth - Nphi)-1))], [-Inf, -Inf, 0.0], [-0.001, Inf, Inf])
function randq(thold, delta)
        th = copy(thold)
        th[end] *= exp(delta*(2rand()-1))
#        th[1] = min(th[1], -0.001)
#        th[3] += delta*(2rand()-1)
        return th
end

function iter(K, n; TC = false, SV = true, DRAW = false, VERBOSE = 1, proptype = :lin, th1 = [-0.1, -0.1, 2.0], delta=0.1)
# n - timepoints each bridge (excluding excluding endpoint)
# K - number of iterations
    if !generated error("Use Inno.generate(Inno.Atan) or Inno.generate(Inno.OU) to generate observation process.") end

    tag = "$model$K$proptype$(M)x$(n)T$Ttotal"*(TC ? "TC":"")
    
    TC2 = TC
    
    Nth = nth(Ptrue)
    Nphi = nphi(Ptrue)

## PRIOR
    #uniform on 1/sigma 
    pth(th) = 1. #pdf(InverseGamma(0.001, 0.001), (th[end])^2)
    qth(th, thold) = 1.

    xi = 0.2*ones(Nphi) # si^2 = 5

    ############
    VERBOSE > 0 && println("Generate bridges (using proposals $proptype): $tag.")
    ############
    srand(3)
    global thetas = zeros(Nth, K)
    global bb = zeros(Bool, M, K)
    global bbp = zeros(M, K)
    global aacc = zeros(Bool, K)
    global aaccp = zeros(K)

    llold = zeros(M)
    ll = zeros(M)
    ll2 = 0.


    global llmax = zeros(M)
    global ll2max = 0

    TT = cell(M)
    TTU = cell(M)
    W = Array(CTPath{1},M)
    W° = Array(CTPath{1},M)
    Y = Array(CTPath{1},M)
    Y° = Array(CTPath{1},M)
    Z = Array(CTPath{1},M)
    Z° = Array(CTPath{1},M)
    U = Array(CTPath{1},M)
    U° = Array(CTPath{1},M)
    Bth = zeros(Nth, M)
    Bth° = zeros(Nth, M)

     
    # allocate arrays for bridges
    dt = Ttotal/M/n
    for m = 1:M 
        Tmin = ttd[m]
        Tmax = ttd[m+1]
        TTU[m] = collect(linspace(0.0, Tmax-Tmin, n))
        if TC
                TT[m] = map(s -> tofs(s, Tmin, Tmax-Tmin), TTU[m])
        else
                TT[m] = collect(linspace(Tmin, Tmax, n))
        end
        W[m] = UvPath(TT[m])
        W°[m] = UvPath(TT[m])    
        Y[m] = UvPath(TT[m])  
        Y°[m] = UvPath(TT[m])
        U[m] = UvPath(TTU[m])
        U°[m] = UvPath(TTU[m])    
        Z[m] = UvPath(TTU[m])  
        Z°[m] = UvPath(TTU[m])

    end

    th = copy(th1)

    srand(3)
    bbsum = 0
    aaccsum = 0
    m1th = 0*th
    m2th = eps2 .+ 0*th

    # conjugate posterior distribution N(mu, Si)
    mu = zeros(Nphi)
    Si = zeros(Nphi, Nphi)

    for k = 1:K
     
        m1th[:] = m1th .+ th
        m2th[:] = m2th .+ th.^2
        VERBOSE > 0 && k > 1 && print("k $k $(mc2str(k, m1th, m2th)) $(repr(round(th,3)))")    
        
        
        #### generate bridges
        P = Model(th...)
    
        for m = 1:M
        
            u = xxd[m]
            v = xxd[m+1]
            
            Tmin = ttd[m]
            Tmax = ttd[m+1]
            T = Tmax-Tmin
     
            if proptype == :brown || proptype == :brownbrown
                Pt = UvAffPro(muof(P), sigma(Tmax,v, P))
                #Pt  = UvLinPro(-1E-7, 0.0, sigma(Tmax,v, P))
            elseif proptype == :lin || proptype == :linlin
                Pt = UvLinPro(Bof(P),betaof(P), sigma(Tmax,v,P))
            end
            
            
            Bth°[:,m] = th 
            if TC
                resample!(Z°[m], Wiener())
#                println(Z°[1])
                u0 = uofx(0., u, T, v, Pt)
                eulerU!(U°[m], u0, Z°[m], Tmin, T, v, Pt, P)
                XofU!(Y°[m], U°[m], Tmin, T, v, Pt)
                ll[m] = llikeliU(U°[m], Tmin, T, v, Pt, P)
                llold[m] = llikeliU(U[m], Tmin, T, v, Pt, P)
            else
                resample!(W°[m], Wiener())
                guidedeuler!(Y°[m], u, W°[m], Tmax, v, Pt, P)
                ll[m] = llikeliXcirc(Y°[m], Pt, P)
                llold[m] = llikeliXcirc(Y[m], Pt, P)
            end

            llmax = map(max, ll, llmax)
            if(k == 1) 
                    Bth[m] = Bth°[m]
                    Z[m].yy[:] = Z°[m].yy
                    U[m].yy[:] = U°[m].yy
                    Y[m].yy[:] = Y°[m].yy
                    Y[m].tt[:] = Y°[m].tt
                    W[m].yy[:] = W°[m].yy
                    llold[m] = ll[m]
            end            
            
            acc = min(1.0, exp(ll[m]-llold[m]))
            VERBOSE > 1 && println("\t acc ", round(acc, 3), " ", round(llold[m],4), " ", round(ll[m],4))
            bbp[m, k] = acc
            bbsum += acc
            bb[m, k] = false
            if rand() <= acc
                bb[m, k] = true
                Bth[m], Bth°[m] = Bth°[m], Bth[m]
	            Y[m], Y°[m] = Y°[m], Y[m]
                W[m], W°[m] = W°[m], W[m]
	            Z[m], Z°[m] = Z°[m], Z[m]
                U[m], U°[m] = U°[m], U[m]
            end 
        end #for m
        
        if DRAW
             I = 5
             p = FramedPlot()
             setattr(p, "xrange", (Y[1].tt[1], Y[I].tt[end]))
             setattr(p, "yrange", (-5.,5.))
             add(p, Points(ttd[1:I+1],xxd[1:I+1], "color","red"))
             add(p, Curve(ttf[1:100*I+1],xxf[1:100I+1], "color","red","linewidth", 0.5))
        end
        
        if DRAW 
           for i in 1:I
                    if bb[m,k]
                        add(p, Curve(Y[i].tt[:],Y[i].yy[:], "color","black", "linewidth", 0.5))
                        add(p, Curve(Y[i].tt[:],Y°[i].yy[:], "color","grey", "linewidth", 0.5))
                    else
                        add(p, Curve(Y[i].tt[:],Y[i].yy[:], "color","grey", "linewidth", 0.5))
                        add(p, Curve(Y[i].tt[:],Y°[i].yy[:], "color","orange", "linewidth", 0.5))
                    end
             end
             Winston.display(p)
        end
        

        P = Model(th...)
        th° = randq(th, delta)
		P° = Model(th°...)
		
		VERBOSE > 0 && k > 1 && print(" ", round(100*bbsum/k/M))    

        #### update sigma and theta
        laccth = 0.

        for m in 1:M
            u = xxd[m]
            v = xxd[m+1]
            
            Tmin = ttd[m]
            Tmax = ttd[m+1]
            T = Tmax-Tmin
     
            if proptype == :brown || proptype == :brownbrown
                Pt = UvAffPro(muof(P), sigma(Tmax,v, P))
                Pt° = UvAffPro(muof(P°),  sigma(Tmax,v, P°))
                #Pt = UvLinPro(-1E-6, 0.0, sigma(Tmax,v, P))
                #Pt° = UvLinPro(-1E-6, 0.0,  sigma(Tmax,v, P°))
            elseif proptype == :lin || proptype == :linlin
                Pt = UvLinPro(Bof(P),betaof(P), sigma(Tmax,v,P))
                Pt° = UvLinPro(Bof(P°),betaof(P°), sigma(Tmax,v, P°))
            end
           

            # compute acceptance probability of change in theta with equivalent innovations
              
            ### bridge Y° (U) using original W (Z)
            Bth°[:,m] = th° 
            if TC2
                u0 = uofx(0., u, T, v, Pt)
                eulerU!(U[m], u0, Z[m], Tmin, T, v, Pt, P)
                XofU!(Y[m], U[m], Tmin, T, v, Pt)

                u0 = uofx(0., u, T, v, Pt°)
                eulerU!(U°[m], u0, Z[m], Tmin, T, v, Pt°, P°)
                XofU!(Y°[m], U°[m], Tmin, T, v, Pt°)
          
                m3 = +llikeliU(U°[m], Tmin, T, v, Pt°, P°)
                m4 = -llikeliU(U[m], Tmin, T, v, Pt, P) # assert that U is a Pt-P-proposal
#                m3 = +llikeliXcirc(Y°[m], Pt°, P°) # for testing: this should yield the same indirectly
#                m4 = -llikeliXcirc(Y[m], Pt, P)
            else
                guidedeuler!(Y[m], u, W[m], Tmax, v, Pt, P)
                guidedeuler!(Y°[m], u, W[m], Tmax, v, Pt°, P°)
                 
                m3 = +llikeliXcirc(Y°[m], Pt°, P°)
                m4 = -llikeliXcirc(Y[m], Pt, P)
            end
            m1 = +lp(Tmin, u, Tmax, v, Pt°) 
            m2 = -lp(Tmin, u, Tmax, v, Pt) 
            laccth += m1 + m2 + m3 + m4   
        end #for m
        ll2max = max(laccth, ll2max)
        accth = min(1., exp(laccth) * pth(th°)/pth(th)*qth(th, th°)/qth(th°, th))
        aacc[k] = false
        aaccp[k] = accth
        VERBOSE > 0 && print(" th° $(repr(round(th°,3))) Pacc $(repr(round(accth,3))) avg ", round(100aaccsum/k,2))  

        if rand() <= accth
            aacc[k] = true
            aaccsum += 1

            th[:] = th°
            for m = 1:M
                # W and Z do not change
                Bth[m] = Bth°[m]
                Y[m], Y°[m] = Y°[m], Y[m]
                U[m], U°[m] = U°[m], U[m]
            end 

        end
        if false
             for i in 1:I
                    if aacc[k]
                        add(p, Curve(Y[i].tt[:],Y[i].yy[:], "color","black", "linewidth", 0.5))
                        add(p, Curve(Y[i].tt[:],Y°[i].yy[:], "color","grey", "linewidth", 0.5))
                    else
                        add(p, Curve(Y[i].tt[:],Y[i].yy[:], "color","grey", "linewidth", 0.5))
                        add(p, Curve(Y[i].tt[:],Y°[i].yy[:], "color","orange", "linewidth", 0.5))
                    end
             end
             Winston.display(p)
        end


        VERBOSE > 0 && print(aacc[k])
         
        #### update conjugate parameters
        P = Model(th...)
        mu[:] = 0.
        Si[:] = 0.

        for m in 1:M
            t = Y[m].tt
            yy = Y[m].yy
            if Nphi == 2
                for i in 1:length(t)-1
                    ainv = 1./sigma(t[i], yy[i], P)^2
                    gdy = (yy[i+1]-yy[i])*ainv
                    gdt = (t[i+1]-t[i])*ainv
                    p1 = phi1(t[i], yy[i], P)
                    p2 = phi2(t[i], yy[i], P)
                    mu[1] += p1*gdy 
                    mu[2] += p2*gdy
                    Si[1,1] +=  p1*p1 * gdt 
                    Si[1,2] +=  p1*p2 * gdt 
                    Si[2,2] +=  p2*p2 * gdt  
                end
                Si[2,1] = Si[1,2]
            else
                #stop()
                for i in 1:length(t)-1
                    ainv = 1/sigma(t[i], yy[i], P)^2
                    phii = phi(t[i], yy[i], P)
                    dy = yy[i+1]-yy[i]
                    mu[:] += phii*dy*ainv
                    Si[:] = Si +  phii*phii' * (t[i+1]-t[i])*ainv
                end
            end
        end #for m

        # sampling conditional posterior
        WW = Si + diagm(xi)
        WL = chol(WW)
        th[1:Nphi] = WL\(WL'\mu + randn(Nphi))
        thetas[:, k] = th
        # th changed, and so the innovations process

        P = Model(th...)
        for m in 1:M        
            u = xxd[m]
            v = xxd[m+1]
            
            Tmin = ttd[m]
            Tmax = ttd[m+1]
            T = Tmax-Tmin

            if proptype == :brown || proptype == :brownbrown
                Pt = UvAffPro(muof(P), sigma(Tmax,v, P))
                #Pt  = UvLinPro(-1E-7, 0.0, sigma(Tmax,v, P))
            elseif proptype == :lin || proptype == :linlin
                Pt = UvLinPro(Bof(P), betaof(P), sigma(Tmax,v,P))
            end


            Bth°[:,m] = th
            if TC
                Z°[m].yy[1] = 0.

                for i in 1:length(U[m].tt)-1
                    s = U[m].tt[i]
                    u1 = U[m].yy[i]
                    t, y = txofsu(s, u1, Tmin, T, v, Pt)
                    ds = U[m].tt[i+1] - s
                    T1 = sqrt(T*(T-s)/2)                    
                    T2 = T1*2/T 
                    T3 = sqrt((T/2)/(T-s))  
                    Z°[m].yy[i+1] = Z°[m].yy[i] -inv(sigma(t, y, P))*(T1*(U[m].yy[i+1] - u1) - 
                                T2*(dotVs(s, T, v, Pt) - b(t,y, P))*ds +
                                T3*(u1 - 2*a(t, y, P)*SDE.Ju(s, T, Pt, u1))*ds)
                end
                
               
                
            else
                W°[m].yy[1] = 0.
                for i in 1:length(Y[m].tt)-1
                    t = Y[m].tt[i]
                    dt = Y[m].tt[i+1] - t
                    y = Y[m].yy[i]
                    W°[m].yy[i+1] = W°[m].yy[i] + inv(sigma(t, y, P))*((Y[m].yy[i+1] - y) - b(t,y, P)*dt - a(t, y, P)*r(t, y, Tmax, v, Pt)*dt)
                end
            end        
        end
        if false #DRAW && TC
            for i in 1:I
                add(p, Curve(Y[i].tt[:],Z[i].yy[:], "color","green", "linewidth", 0.5))
                add(p, Curve(Y[i].tt[:],Z°[i].yy[:], "color","blue", "linewidth", 0.5))
            end
            Winston.display(p)
        end
        for m in 1:M        
            Bth[m] = Bth°[m]
            # Y and U do not change
            W[m], W°[m] = W°[m], W[m]
            Z[m], Z°[m] = Z°[m], Z[m]
        end
        println()


    end # for K

    if SV
        if !isdir(tag)
             mkdir(tag)
        end
        writecsv("$tag/thetas$tag.csv", thetas') 
        accs = [aacc aaccp bb' bbp']
        writecsv("$tag/accs$tag.csv", accs)
        open("$tag/acc$tag.txt", "w") do f
            println(f,"qufull $quxx quobs $quxxd")
            println(f,"estimate: $(repr(mean(thetas,2)[:]))")
            println(f,"est: $(mc2str(K, m1th, m2th))") 
            println(f,"bridge acc ",  repr(round(100*mean(1.0bb),2)), "%, avg prob  ", repr(round(100*mean(bbp),2)) )
            println(f, "llmax ", maximum(llmax))
            #println(f,"bridge acc ",  repr(round(100*mean(1.0bb,2),2)), "%, avg prob  ", repr(round(mean(bbp, 2),2)) )
            println(f,"sigma acc ",  round(100*sum(aacc)/K,2), "% avg prob ", round(100*sum(aaccp)/K, 2), "% " )
        end
    end

    println("estimate: $(repr(mean(thetas,2)[:]))")
    println("bridge acc ",  round(100*mean(1.0bb,2),2), "%, avg prob  ", repr(round(100*mean(bbp, 2),2)) )
    println("sigma acc ",  round(100*sum(aacc)/K,2), "% avg prob ", round(100*sum(aaccp)/K, 2), "% " )
    return thetas, accs
end


function runexamples()
    generate(Atan)
    for n in [10, 100,1000], t in [:lin], tc in [true, false]
        println("$n $t $tc")
        iter(10000,n; proptype=t, DRAW=true, TC=tc, th1=[-.1, -.1, 2.])
    end    
end




   
end #module
                       
