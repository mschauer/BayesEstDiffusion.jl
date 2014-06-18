module ProkarC #prokaryotic auto-regulation using conjugate proposals
export iter, PAR, generate, runexamples, load
using Winston
include("SDE.jl")
using .SDE
include("misc.jl")
import .SDE.sigma, .SDE.b, .SDE.a

######## data?
status = :nodata

####### VOU
  
immutable PAR <: MvPro
    th::Vector{Float64}
    d
    d2
    k::Float64   # number of copies of gene
    PAR(th) = new(th, 4, 8, 10.)
end
nthI(_::PAR) = 8
nthII(_::PAR) = 8
nth(_::PAR) = 8

const S =
   [ 0. 0  1  0  0  0 -1  0
     0  0  0  1 -2  2  0 -1
    -1  1  0  0  1 -1  0  0
    -1  1  0  0  0  0  0  0] 
const St = S'

h(x, P::PAR) = max([max(x[3],0.)*x[4], P.k - x[4], x[4], x[1], 0.5*abs(x[2])*(x[2]-1),  x[3],  x[1],  x[2]], 0.)
thh(x, P::PAR) = P.th .* h(x,  P)
b(t,x, P::PAR) = S*thh(x, P)
sigma(t,x, P::PAR) = sqrt(max(thh(x, P),0.))' .* S
function a(t,x, P::PAR) 
    s = sigma(t,x,P)
    s*s' 
end

phi(t, x, P::PAR) = S .* h(x, P)'



#const U2 = 8. 
#const I2 = -3. 
#const U3 = 4.
#const U4 = 6.
#const I34 =0.

# 29.7586 5.02829 6.044613
# -52.47444 14.83634

#const U2 = 15.
#const I2 = -52.
#const U4 = 5.
#const U3 = 6.
#const I34 = 30.

const U2 = 18.
const I2 = -82.
const U4 = 5.
const U3 = 6.
const I34 = -32.



Bof(P) = 
[   (-P.th[7])    0                   0                   P.th[3]
    P.th[4]     (-P.th[8]-U2*P.th[5]) 2P.th[6]            0
    0           0.5*U2*P.th[5]      (-P.th[6]-U4*P.th[1]) (-P.th[2] - U3*P.th[1])
    0           0                   -U4*P.th[1]           (-P.th[2] - U3*P.th[1])]

betaof(P) = [0, -I2*P.th[5],0.5I2*P.th[5] + P.k*P.th[2] - I34*P.th[1], P.k*P.th[2] - I34*P.th[1] ]

muof(P) = zeros(4)

####### Prior
import Base.clamp;
clamp(x::Vector,l::Vector, u::Vector) = map(i -> clamp(x[i], l[i], u[i]), 1:length(x))


function generate(Model_; thtrueI = [],
                          thtrueII = [0.1, 0.7, 0.35, 0.2, 0.1, 0.9, 0.3, 0.1],
u = [7.5, 10., 5., 5.], Nobs = 51, N = 400000 + 1, Delta=1., SV=false)
    global Model = Model_
    th = [thtrueI, thtrueII]
    global Ptrue = Model(th)
    global model = string(Model)

    global M = Nobs - 1 # number of bridges
    global Ttotal = Delta*M      #total time span

    tag = "$model$(M)T$Ttotal"
    println(tag)
    d = Ptrue.d
    d2 = Ptrue.d2
    ############

    ############
    print("Generate X")
    ############

    srand(3)
    tt = linspace(0., Ttotal, N)
    W = sample(tt, Wiener(d2))
    xx = euler(u, W, Ptrue).yy
    println(".")
    
    global ttf = tt[1:(N-1)//M/100:end]
    global xxf = xx[:, 1:(N-1)//M/100:end]
    global ttd = tt[1:(N-1)//M:end]
    global xxd = xx[:, 1:(N-1)//M:end]

    global status = :generated
    
    if SV
        writecsv("xtrue$(tag).csv",[linspace(0, Ttotal, size(xx,2)) xx'])
        writecsv("xobs$(tag).csv",[linspace(0, Ttotal, size(xxd,2)) xxd'])
    end
end

function load(file="autoreg50fo", Model_=PAR; thtrueI = [],
            thtrueII = [0.1, 0.7, 0.35, 0.2, 0.1, 0.9, 0.3, 0.1],
            Delta=1.)
    th = [thtrueI, thtrueII]

    const global Model = Model_
    const global Ptrue = Model(th)
    const global model = string(Model)

    const global xxd = readcsv(string(file,".csv"))'
    const global Nobs = size(xxd,2)
    const global M = Nobs - 1 # number of bridges
    const global Ttotal = Delta*M      #total time span
    const global ttd = linspace(0., Ttotal, Nobs)

    tag = "$model$file"
    println(tag)

    global status = :empir
    
end


function randq(th, delta, P)
        th .* exp(delta*(randn(length(th))))
end

function iter(K, n;             # numer of iterations, number of points for imputed bridges, including endpoints
proptype = :lin,                # type of proposals, :lin or :brown
TC = false,                     # apply timechange
PRECOMP = false,                # use tabulated Phim
SV = true, DRAW = false, VERBOSE = 1, REJNEG = true,
th1 = 0.05ones(8),
delta=0.05)
    NEEDX = true
    if status == :nodata error("Use generate(Model) or load(...) to generate or load observation process.") end

    tag = "C$model$K$proptype$(M)x$(n)T$Ttotal"*(TC ? "TC":"")

    de_ = 0.
    Nth = nth(Ptrue)
    NthII = nthII(Ptrue)
    NthI = nthI(Ptrue)
    d = Ptrue.d
    d2 = Ptrue.d2

## PRIOR
    pth(th) = prod(1./th) 
    qth(th, thold) = 1.


    xi = 0.1*ones(NthI)

    ############
    VERBOSE > 0 && println("Generate bridges (using proposals $proptype): $tag.")
    ############
    srand(3)
    thetas = zeros(Nth, K)
    bb = zeros(Bool, M, K)
    bbp = zeros(M, K)
    aacc = zeros(Bool, K)
    aaccp = zeros(K)
    Td = zeros(M+1)

    llold = zeros(M)
    ll = zeros(M)
    ll2 = 0.

    laccth = 0.
    llmax = zeros(M)
    ll2max = 0
    rejected = 0
    rejectedth = 0

    TT = cell(M)
    TTU = cell(M)
    W = Array(CTPath{2},M)
    W° = Array(CTPath{2},M)
    Y = Array(CTPath{2},M)
    Y° = Array(CTPath{2},M)
    Z = Array(CTPath{2},M)
    Z° = Array(CTPath{2},M)
    U = Array(CTPath{2},M)
    U° = Array(CTPath{2},M)
    Phim = cell(n-1)
    Phim° = cell(n-1)

     
    # allocate arrays for bridges
    dt = Ttotal/M/n
    Tall = Ttotal/M
    TTUall = linspace(0.0, Tall, n)
    for m = 1:M 
        Tmin = ttd[m]
        Tmax = ttd[m+1]
        TT[m] = linspace(Tmin, Tmax, n)
        TTU[m] = linspace(0.0, Tmax-Tmin, n)
        W[m] = MvPath(TT[m],d2)
        W°[m] = MvPath(TT[m],d2)    
        Y[m] = MvPath(TT[m],d)  
        Y°[m] = MvPath(TT[m],d)
        U[m] = MvPath(TTU[m],d)
        U°[m] = MvPath(TTU[m],d)    
        Z[m] = MvPath(TTU[m],d2)  
        Z°[m] = MvPath(TTU[m],d2)

    end

    P1 = Model(th1)
    th = copy(th1)
    th° = copy(th)
    println("th1: ", repr(th))
    qq = 1.0;

    srand(3)
    bbsumall = 0
    bbsum = zeros(M)
    aaccsum = 0
    mth = 0*th

    # conjugate posterior distribution N(mu, Si)
    mu = zeros(NthI)
    Si = zeros(NthI, NthI)

    for k = 1:K
        mth[:] = mth .+ th
        VERBOSE > 0 && k > 1 && print("k $k ", repr(round(mth./k,3)), " ",repr(round(th,3)))    
        
        
        #### generate bridges
        P = Model(th)
        if PRECOMP
            for i = 1:n-1
                Phim[i] = Phims(TTUall[i], Tall, Bof(P))
            end 
        end
        for m = 1:M
        
            u = xxd[:, m]
            v = xxd[:, m+1]
            
            Tmin = ttd[m]
            Tmax = ttd[m+1]
            T = Tmax-Tmin
     
            if proptype == :brown || proptype == :brownbrown
                Pt = MvAffPro(muof(P), sigma(Tmax,v, P))
                #Pt  = MvLinPro(-1E-7, 0.0, sigma(Tmax,v, P))
            elseif proptype == :lin || proptype == :linlin
                Pt = MvLinPro(Bof(P),betaof(P), sigma(Tmax,v,P))
            end

            posi = false
            while !posi
                        
                if TC && PRECOMP
                    resample!(Z°[m], Wiener(d2))
                    u0 = uofx(0., u, T, v, Pt, Phim[1])
                    eulerU!(U°[m], u0, Z°[m], Tmin, T, v, Pt, P, Phim)
    NEEDX &&        XofU!(Y°[m], U°[m], Tmin, T, v, Pt)
                    ll[m] = llikeliU(U°[m], Tmin, T, v, Pt, P, Phim)
                    llold[m] = llikeliU(U[m], Tmin, T, v, Pt, P, Phim)
                elseif TC
                    resample!(Z°[m], Wiener(d2))
                    u0 = uofx(0., u, T, v, Pt)
                    eulerU!(U°[m], u0, Z°[m], Tmin, T, v, Pt, P)
    NEEDX &&        XofU!(Y°[m], U°[m], Tmin, T, v, Pt)
                    ll[m] = llikeliU(U°[m], Tmin, T, v, Pt, P)
					llold[m] = llikeliU(U[m], Tmin, T, v, Pt, P)
                else
                    resample!(W°[m], Wiener(d2))
                    guidedeuler!(Y°[m], u, W°[m], Tmax, v, Pt, P)
                    ll[m] = llikeliXcirc(Y°[m], Pt, P)
                    llold[m] = llikeliXcirc(Y[m], Pt, P)
                end
                posi = true
                if REJNEG && any(Y°[m].yy .<= 0)
                    posi = false 
                    print("\$")
                    rejected += 1
                end

            end

            llmax = map(max, ll, llmax)
            if(k == 1) 
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
            bbsum[m] += acc
            bbsumall += acc
            bb[m, k] = false
            if rand() <= acc
                bb[m, k] = true
	            Y[m], Y°[m] = Y°[m], Y[m]
                W[m], W°[m] = W°[m], W[m]
	            Z[m], Z°[m] = Z°[m], Z[m]
                U[m], U°[m] = U°[m], U[m]
            end 
        end #for m

        if DRAW
             I1 = 1 + k % 40
             I2 = I1 + 9
             p = FramedPlot()
             setattr(p, "xrange", (I1, I2))
             setattr(p, "yrange", (-3,20))

             for j in 1:d
                 add(p, Points(ttd[I1:I2], xxd[j, I1:I2+1][:], "color","red"))
                 for i in I1:I2
                     add(p, Curve(Y[i].tt[:],Y[i].yy[j, :][:] , "color","black", "linewidth", 0.5))
                 end
             end
             display(p)
        end
        if false
             C1 = 3
             C2 = 4
             if k % 4 == 0; C1 = 1; C2 = 4; end
             I1 = 1 + k % 40
             I2 = I1 + 9
             p = FramedPlot()
             const Rs = [(-1, 15),(-1, 30),(-5, 15),(-5, 15)]
             setattr(p, "xrange", Rs[C1])
             setattr(p, "yrange", Rs[C2])
             add(p, Points(xxd[C1, I1:I2+1][:], xxd[C2, I1:I2+1][:], "color","red"))
             for i in I1:I2
                 add(p, Curve(Y[i].yy[C1, :][:],Y[i].yy[C2, :][:] , "color","black", "linewidth", 0.5))
             end
             display(p)
        end
		
		VERBOSE > 0 && k > 1 && print(" ", round(100*bbsumall/k/M,1),"(min ", round(100*minimum(bbsum)/k,1))
		VERBOSE > 0 && k > 1 && @printf(" ll %.1f)",round(maximum(ll),1))    
		
		P = Model(th)

        posi = false
        while !posi
            posi = true
		    th°[:] = th
	        try
    	        P = Model(th)		
			    #### compute conjugate proposal
	            mu[:] = 0.
	            Si[:] = 0.
	
	            for m in 1:M
	                t = Y[m].tt
	                yy = Y[m].yy
	                for i in 1:length(t)-1
	                    phii = phi(t[i], yy[:,i], P)
	                    si = sigma(t[i], yy[:,i], P)
	                    aiph= (si*si')\phii
	                    dy = yy[:,i+1] .- yy[:,i]
	                    Si[:] = Si +  phii'*aiph * (t[i+1]-t[i])
	                end
	            end #for m
	            # sampling conditional posterior
	            WW = Si .+ diagm(xi)
	            WL = chol(WW)
	            de = WL\randn(NthI)
	            de_ += de[1]/de[2]
	            print(round(de_/k,2))
	            th° = th .+ delta*de
	            th° = min(1100.,max(0.00091, th°))
	        
	            mu[:] = 0.
	            Si[:] = 0.
			    P° = Model(th°)
	            for m in 1:M
	                t = Y[m].tt
	                yy = Y[m].yy
	                for i in 1:length(t)-1
	                    phii = phi(t[i], yy[:,i], P°)
	                    si = sigma(t[i], yy[:,i], P°)
    	                aiphii= (si*si')\phii
	                    dy = yy[:,i+1] .- yy[:,i]
	                    Si[:] = Si +  phii'*aiphii * (t[i+1]-t[i])
	                end
	            end #for m
	            WW° = Si .+ diagm(xi)
	            WL° = chol(WW°)
	            qq = prod(diag(WL°))/prod(diag(WL))*exp( -0.5/delta^2*norm(WL°*(th-th°))^2 + 0.5/delta^2*norm(WL*(th-th°))^2)
	            print(" q ", round(qq,2))        
    		catch e
	            print("\n\n\t\t")
	            println(e)
	            th°[:] = th
	            qq = 1.0
		    end
        
			P° = Model(th°)
		    VERBOSE > 0 && print(" th° ", repr(round(th°,3)))
		
            if PRECOMP
                for i = 1:n-1
                    Phim°[i] = Phims(TTUall[i], Tall, Bof(P°))
                end 
            end
      
            #### update sigma and theta
            laccth = 0.
        

            for m in 1:M
                u = xxd[:, m]
                v = xxd[:, m+1]
            
                Tmin = ttd[m]
                Tmax = ttd[m+1]
                T = Tmax-Tmin
     
                if proptype == :brown || proptype == :brownbrown
                    Pt = MvAffPro(muof(P), sigma(Tmax,v, P))
                    Pt° = MvAffPro(muof(P°),  sigma(Tmax,v, P°))
                elseif proptype == :lin || proptype == :linlin
                    Pt = MvLinPro(Bof(P),betaof(P), sigma(Tmax,v,P))
                    Pt° = MvLinPro(Bof(P°),betaof(P°), sigma(Tmax,v, P°))
                end

                # compute acceptance probability of change in theta with equivalent Innovations
              
                ### bridge Y° (U) using original W (Z)
                if TC && PRECOMP
                    u0 = uofx(0., u, T, v, Pt°, Phim°[1])
                    eulerU!(U°[m], u0, Z[m], Tmin, T, v, Pt°, P°, Phim°)
    NEEDX &&        XofU!(Y°[m], U°[m], Tmin, T, v, Pt°)
                    m1 = +lp(Tmin, u, Tmax, v, Pt°) 
                    m2 = -lp(Tmin, u, Tmax, v, Pt) 
                    m3 = +llikeliU(U°[m], Tmin, T, v, Pt°, P°, Phim°)
                    m4 = -llikeliU(U[m], Tmin, T, v, Pt, P, Phim)
                    laccth += m1 + m2 + m3 + m4
                elseif TC
                    u0 = uofx(0., u, T, v, Pt°)
                    eulerU!(U°[m], u0, Z[m], Tmin, T, v, Pt°, P°)
    NEEDX &&        XofU!(Y°[m], U°[m], Tmin, T, v, Pt°)
                    m1 = +lp(Tmin, u, Tmax, v, Pt°) 
                    m2 = -lp(Tmin, u, Tmax, v, Pt) 
                    m3 = +llikeliU(U°[m], Tmin, T, v, Pt°, P°)
                    m4 = -llikeliU(U[m], Tmin, T, v, Pt, P)
                    laccth += m1 + m2 + m3 + m4
                else
                    guidedeuler!(Y°[m], u, W[m], Tmax, v, Pt°, P°)
                    m1 = +lp(Tmin, u, Tmax, v, Pt°) 
                    m2 = -lp(Tmin, u, Tmax, v, Pt) 
                    m3 = +llikeliXcirc(Y°[m], Pt°, P°)
                    m4 = -llikeliXcirc(Y[m], Pt, P)
                    laccth += m1 + m2 + m3 + m4
                end
                if REJNEG && any(Y°[m].yy .<= 0)
                        posi = false 
                        print("\$")
                        rejectedth += 1
                end

            end #for m
        
        end #while !posi

        ll2max = max(laccth, ll2max)
        accth = min(1., exp(laccth) * pth(th°)/pth(th)*qq)   
        aacc[k] = false
        aaccp[k] = accth
        VERBOSE > 0 && print(" Pacc $(repr(round(accth,3))) avg ", round(100aaccsum/k,2))  
        

        if rand() <= accth
            aacc[k] = true
            aaccsum += 1

            th[:] = th°
            for m = 1:M
                # W and Z do not change
                Y[m], Y°[m] = Y°[m], Y[m]
                U[m], U°[m] = U°[m], U[m]
            end 

        end
		VERBOSE > 0 && @printf(" ll %.1f %s", round(laccth,1), aacc[k] ? "acc" : "")             

		thetas[:, k] = th
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
            println(f,"estimate: $(repr(mean(thetas,2)[:]))")
            println(f,"bridge acc ",  repr(round(100*mean(1.0bb),2)), "%, avg prob  ", repr(round(100*mean(bbp),2)) )
            println(f, "llmax ", maximum(llmax), "rejected neg. path", rejected)
            println(f,"sigma acc ",  round(100*sum(aacc)/K,2), "% avg prob ", round(100*sum(aaccp)/K, 2), "% " )
        end
    end

    println("estimate: $(repr(mean(thetas,2)[:]))")
    println("bridge acc ",  round(100*mean(1.0bb,2),2), "%, avg prob  ", repr(round(100*mean(bbp, 2),2)) )
    println("sigma acc ",  round(100*sum(aacc)/K,2), "% avg prob ", round(100*sum(aaccp)/K, 2), "% " )
    return thetas, accs
end


function runexamples(k=100000)
    load("autoreg50fo", PAR)
    for n in [20, 10, 50], t in [:lin], tc in [true]
        println("$n $t $tc")
        iter(k, n; proptype=t, DRAW=false, TC=tc, PRECOMP=true, delta=0.8, REJNEG=true)
    end    
end


end #module
                       
