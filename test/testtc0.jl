using SDE
using Cubature
import SDE.b, SDE.sigma, SDE.a
#srand(7)
if !isdefined(:Model)

type Model <: UvPro
    B
    beta
    Sigma
    A
    d::Int   
end
end

b(t,x, P::Model) =   exp(-0.2*t)*P.B*x + P.beta
sigma(t,x, P::Model) = exp(-0.1*(t))*P.Sigma
a(t,x,  P::Model) = exp(-0.2*(t))*P.A

let d = 1

Q95 = sqrt(2.)*erfinv(0.95)
function mc(M)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  round([m, Q95*ste], floor(Int,2-log(10,Q95*ste)))
end

tmax = 0.6
tmin = 0.3
T = tmax - tmin
 
B0 = -2.0
beta0 = 0.0
sigma0 = .75
ph0(t,s ) =  5.*exp(-0.2*s)-5.*exp(-0.2*t)




u = 0.3
v = 0*0.5

B = exp(-0.2*tmax)*B0 
 
N = 1000
 

P = Model(B0, beta0,sigma0, sigma0^2, 1)
Pt = UvLinPro(B,  beta0, sigma(tmax, v, P))
#Pt = UvLinPro(-1E-6, 0.0, sigma(tmax, v, P))
#Pt = UvAffPro(0.0, sigma(tmax, v, P))


tt = linspace(tmin, tmax, N)
dt = (tt[N]-tt[1])/(N-1)

ss = linspace(0.,T, N)
ds = (ss[N]-ss[1])/(N-1)

ttofss = map(s -> tofs(s, tmin, T), ss)

M = 10000

Yll = zeros(M) 
UYll = zeros(M)

Y2ll = zeros(M) 
UY2ll = zeros(M) 

Ull = zeros(M) 
YUll = zeros(M)

Yt0 = zeros(M)
Y2t0 = zeros(M)
YUt0 = zeros(M)
UYt0 = zeros(M)
Us0 = zeros(M)
 
nt = floor(Int,5N/10)


t0 =  tmin + dt*(nt-1)
W = SDE.sample(tt, Wiener())

@time for i in 1:M
	resample!(W, Wiener())
	Y = guidedeuler(u, W, tmax, v, Pt, P)
	UY = UofX(Y, tmin, T, v, Pt) 
	ll = llikeliXcirc(Y, Pt, P)
	ll2 = llikeliU(UY, tmin, T, v, Pt, P)
	Yt0[i] = Y.yy[nt]
#	println(U.tt[ns], s0)
    s0 = UY.tt[nt]
	UYt0[i] = xofu(UY.tt[nt], UY.yy[nt], T, v, Pt) 
	Yll[i] = ll
	UYll[i] = ll2
end	

#       Y               UY                               U                       YU
println("\ntesttc0.jl:")
println("X° - proposal, U(X°) - timechanged X° proposal, U - proposal in U-time, X°(U) - undone timechange, X°° - X° at same times as Us")

println("X° at ", round(t0,3), "($nt): ""[", repr(mc(Yt0)),"]")
println("U(X°) at ", round(tofs(s0, tmin, T),3), "($nt): ""[", repr(mc(UYt0)),"]")

s0 = soft(t0, tmin, T)
ns = 1 + floor(Int,s0/ds)
print(s0)
s0 = ds*(ns-1)
println(" ", s0)


@Test.test abs(1.-SDE.J(T-0.001, T, Pt)*a(tmax, 0., Pt)) < 0.001

u0 = uofx(0., u,  T, v, Pt)
W = SDE.sample(ss, Wiener())
U = eulerU(u0, W, tmin, T, v, Pt, P)
@time for i in 1:M
 	resample!(W, Wiener())
	U = eulerU(u0, W, tmin, T, v, Pt, P)
	YU = XofU(U, tmin, T, v, Pt) 
	llu = llikeliU(U, tmin, T, v, Pt, P)
	llyu = llikeliXcirc(YU, Pt, P)
	Us0[i] = U.yy[ns]
	YUt0[i] = xofu(s0, U.yy[ns],  T, v, Pt)
	Ull[i] = llu
	YUll[i] = llyu
end



println("X°(U) at ", round(tofs(s0,tmin, T),3), "($ns): ",  "[", repr(mc(YUt0[:])),"]")

 
W = SDE.sample(ttofss, Wiener())
@time for i in 1:M
	resample!(W, Wiener())
	Y2 = guidedeuler(u, W, tmax, v, Pt, P)
	Y2ll[i]  = llikeliXcirc(Y2, Pt, P)
	Y2t0[i] = Y2.yy[ns]
	U2 = UofX(Y2, tmin, T, v, Pt) 
	UY2ll[i] = llikeliU(U2, tmin, T, v, Pt, P)
    
end	

println("X°° at ", round(ttofss[ns],3), "($ns): ", "[", repr(mc(Y2t0[:])),  "]"	)
	 			 



    pt = exp(lp(tmin, u, tmax, v, Pt))
	p = exp(SDE.varlp(tmin, u, tmax, v, ph0, B0, beta0, s->a(s, u, P)))
	println("p ", round(p,5), " ", repr(mc(exp(Yll)*pt)), repr(mc(exp(Ull)*pt)) , repr(mc(exp(Y2ll)*pt)))
	println("~p ", round(pt,5))
	w1 = mc(pt/p*exp(Yll))
	w1b = mc(pt/p*exp(UYll))
	
	w2 = mc(pt/p*exp(Ull))
	w2b = mc(pt/p*exp(YUll))
	w3 = mc(pt/p*exp(Y2ll))
    w3b = mc(pt/p*exp(UY2ll))	
    
	println(" Σ ~p/p*psi(X°) = ", repr(w1), " ≈ 1\n",
	        " Σ ~p/p*psi(UX°) = ", repr(w1b), " ≈ 1\n",
            " Σ ~p/p*psi(U) = ", repr(w2), " ≈ 1\n",
            " Σ ~p/p*psi(YU) = ", repr(w2b), " ≈ 1\n", 
            " Σ ~p/p*psi(Xs)= ", repr(w3), " ≈ 1\n",
            " Σ ~p/p*psi(UXs)= ", repr(w3b), " ≈ 1"
            )
	println("lmax ", maximum(exp(Yll)), " ", maximum(exp(Ull)), " ", maximum(exp(Y2ll)))

    Test.@test abs(w1[1] - 1.) < 1.96w1[2]
    Test.@test abs(w2[1] - 1.) < 1.96w2[2]
    Test.@test abs(w3[1] - 1.) < 1.96w3[2]

end

