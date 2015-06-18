using SDE
using Cubature
import SDE.b, SDE.sigma, SDE.a

srand(9)
if !isdefined(:MvTest)

type MvTest <: MvPro
    B
    beta
    Sigma
    A
    d::Int   
end
end

b(t,x, P::MvTest) =   exp(-0.2*t)*P.B*x + P.beta
sigma(t,x, P::MvTest) = exp(-0.1*(t))*P.Sigma
a(t,x,  P::MvTest) = exp(-0.2*(t))*P.A

d = 2

const Q95 = sqrt(2.)*erfinv(0.95)
function mc(M)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  round([m, Q95*ste], round(Int,2-log(10,Q95*ste)))
end

@Test.test 5.0 == soft(tofs(5., 1., 10.,), 1.,10.)

T =  1.7
tmax = T
tmin = 0.
 
B0 = [ -2.21225   -0.49590;   0.631753  -1.34795]
beta0 = [0.5,-0.5]
sigma0 = .1* [1.5*[1, 1]  1.5*[1,-1]]
ph0(t,s ) =  5.*exp(-0.2*s)-5.*exp(-0.2*t)


begin 



    u = [1.,0.]
    v = [0.5, -0.2]

    #taken slightly off 
    B = 1.2exp(-0.2*T)*B0 
 
    N = 100
 

    P = MvTest(B0, beta0,sigma0, sigma0*sigma0', 2)
    Pt = MvLinPro(B, beta0, sigma(T, v, P))



    tt = linspace(0., T, N)
    dt = (tt[N]-tt[1])/(N-1)

    ss = linspace(0.,T, N)
    ds = (ss[N]-ss[1])/(N-1)

    ttofss = map(s -> tofs(s, 0, T),ss)

    M = 500

    Yll = zeros(M) 
    Y2ll = zeros(M) 
    Ull = zeros(M) 
    Ull2 = zeros(M) 

    Yt0 = zeros(2, M)
    Y2t0 = zeros(2, M)
    YofUt0 = zeros(2, M)
    Us0 = zeros(2, M)
 
    nt = floor(Int,5N/10)


    t0 =  dt*(nt-1)
    W = SDE.sample(tt, Wiener(d))
    @time for i in 1:M
	    resample!(W, Wiener(d))
	    Y = guidedeuler(u, W, T, v, Pt, P)
	    ll = llikeliXcirc(Y, Pt, P)
	    Yt0[:,i] = Y.yy[:,nt]
	    Yll[i] = ll
    end	

    println("X° at ", round(t0,3), "($nt): ", "[", repr(mc(Yt0[1,:])), ", ",repr(mc(Yt0[2,:])), "]"	)


    x1 = zeros(2, M)
    s0 = soft(t0, 0, T)
    ns = 1 + floor(Int,s0/ds)
    s0 = ds*(ns-1)


    u0 = uofx(0.0, u,  T, v, Pt)
    W = SDE.sample(ss, Wiener(d))
    U = eulerU(u0, W, 0., T, v, Pt, P)
    srand(9)
    @time for i in 1:M
 	    resample!(W, Wiener(d))
	    eulerU!(U, u0, W, 0., T, v, Pt, P)
	    ll =  llikeliU(U, 0., T, v, Pt, P)
	    Us0[:,i] =  U.yy[:,ns]
	    YofUt0[:,i] =  xofu(s0, U.yy[:,ns],  T, v, Pt)
	    Ull[i] = ll
    end
    srand(9)
    #Profile.init(10^7, 0.0003)
    #@profile  begin
    @time begin
    Phim = [ Phims(ss[i], T, Pt.B) for i in 1:N-1] #precompute Phim
    U2 = eulerU(u0, W, 0., T, v, Pt, P, Phim)
    for i in 1:M
 	    resample!(W, Wiener(d))
	    eulerU!(U2, u0, W, 0., T, v, Pt, P, Phim)
	    ll2 = llikeliU(U2, 0., T, v, Pt, P, Phim)
	    Us0[:,i] =  U2.yy[:,ns]
	    xofu(s0, U2.yy[:,ns],  T, v, Pt, Phim[ns])
	    Ull2[i] = ll2
    end
    end

    println("diff normal/precompiled: ", norm(Ull .- Ull2))

    println("X°(U) at ", round(tofs(s0,0., T),3), "($ns): ",  "[", repr(mc(YofUt0[1,:])), ",", repr(mc(YofUt0[2,:])), "]")




 
    W = SDE.sample(ttofss, Wiener(d))
    @time for i in 1:M
	    resample!(W, Wiener(d))
	    Y2 = guidedeuler(u, W, T, v, Pt, P)
	    ll = llikeliXcirc(Y2, Pt, P)
	    Y2t0[:,i] = Y2.yy[:,ns]
	    Y2ll[i] = ll
    end	

    println("X°° at ", round(ttofss[ns],3), "($ns): ", "[", repr(mc(Y2t0[1,:])), ", ", repr(mc(Y2t0[2,:])), "]"	)
	 			 



    pt = exp(lp(0., u, T, v, Pt))
	p = exp(SDE.varlp(0., u, T, v, ph0, B0, beta0, s->a(s, u, P)))
	println("p ", round(p,5), " ", repr(mc(exp(Yll)*pt)), repr(mc(exp(Ull)*pt)) , repr(mc(exp(Y2ll)*pt)))
	println("~p ", round(pt,5))
	w1 = mc(pt/p*exp(Yll))
	w2 = mc(pt/p*exp(Ull))
	w2b = mc(pt/p*exp(Ull2))
	w3 = mc(pt/p*exp(Y2ll))
	println(" Σ ~p/p*psi(X°) = ", repr(w1), " ≈ 1\n", " Σ ~p/p*psi(U) = ", repr(w2), repr(w2b), " ≈ 1\n", " Σ ~p/p*psi(X°°)= ", repr(w3), " ≈ 1")
	println("lmax ", maximum(exp(Yll)), " ", maximum(exp(Ull)), " ", maximum(exp(Y2ll)))
	

    Test.@test abs(w1[1] - 1.) < 1.96w1[2]
    Test.@test abs(w2[1] - 1.) < 1.96w2[2]
    Test.@test abs(w3[1] - 1.) < 1.96w3[2]

end

