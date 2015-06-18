#include("../SDE.jl")
using SDE
using Base.Test
using Randm
#include(Pkg.dir("SDE","src", "quad.jl"))
import SDE.eps2

function mc(M::Vector)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  m, ste
end
qu = sqrt(2)*erfinv(2*0.995-1)

macro repeat(ex, n)
    quote
       r1 = $ex
       y = zeros(size(r1)..., $n)
       s = prod(size(r1))
       y[1:s] = r1
       for i in 1:$n-1
        y[i*s+1:(i+1)*s] = $ex
       end
       y
    end
end
      
for d in 1:3
    N = 10000
    xi = 0.3
	global T, t, v, h, P, B, A, s
	println("\ntestlinproc.jl: dimension $d xi $xi")

	if (d== 1) srand(7)
	elseif (d== 2) srand(8)
	else srand(6)
	end
	B = xi.*Randm.randstable(d)
	Si = rand(d,d)
	A = Si*Si'
	beta = 0.5*randn(d)
	
	P = SDE.MvLinPro(B, beta, Si)
	P0 = MvAffPro(beta, Si)
	P00 = MvLinPro(eye(d)*1E-7, beta, Si)
	
	function aa(s)
	 A
	end
	

	T = 1.7
	t = 0.5
	v = zeros(d)



#	println("Test SDE.H, SDE.K")
	

    H1 = SDE.H(t,T, P)
    H2 = SDE.varH(t, T, B, aa)
    K = SDE.K(t, T, P)
    H3 = expm((T-t)*B')*inv(K)*expm((T-t)*B)
    
    H0 = SDE.H(t, T, P0)
    K0 = SDE.K(t, T, P0)
    H00 = SDE.H(t, T, P00)
    @test norm(H0 - H00) < 1E-5*d^2
    

	h = T-t
	phi = expm(h*B)
	phim = expm(-h*B)
	K2 = P.lambda - phi*P.lambda*phi'
	H4 = phi'*inv(P.lambda - phi*P.lambda*phi')*phi
	
	@test norm(H1 - H2) < 1E-10
	@test norm(H1 - H3) < 1E-10
	@test norm(H1 - H4) < 1E-10

	@test norm(K-K2) < d*24*eps()
    

	function varr3(t, x, T, v, P)
	    h = T - t
		binv = inv(P.B)
		phi = expm(h*P.B)
		phim = expm(-h*P.B)
		vt =  phim*( v) + phim * binv * beta -binv*beta 
		SDE.H(t, T, P) * (vt-x )
	
	end

#	println("Test r")

	srand(5)
	x0 = x = randn(d)/d
	v = randn(d)
	r1 = SDE.r(t, x, T, v, P)
	r2 = SDE.varr(t, x, T, v, P)
	@test norm(r1 - r2) < 1E-10

	r3 = varr3(t, x, T, v, P)
	@test norm(r1 - r3) < 1E-10
	
	r10 = SDE.r(t, x, T, v, P0)
	r100 = SDE.r(t, x, T, v, P00)
	@test norm(r10 - r100) < 1E-5*d^2
	

#	println("Test mu")

	mu = SDE.mu(t, x, T, P)
	mu0 = SDE.mu(t, x, T, P0)
	mu00 = SDE.mu(t, x, T, P00)
	mu2 = SDE.varmu(t, x, T, P)
	@test norm(mu - mu2) < 1E-13
	@test norm(mu0 - mu00) < 1E-6


    println("Testing samplep(::LinPro)")

    y = @repeat(SDE.samplep(t, x, T, P), N)
    mcy = mc(sum(y,1)[:])
	println(repr(mcy), " ", sum(mu))
    @test abs(mcy[1]- sum(mu)) < qu*mcy[2]
    muhat = mean(y,2)
	Khat = cov((y.-x)')
	println("K vs. sample covariance ", norm(Khat - K), " < ", 0.04*d)
	@test norm(Khat-K) < 0.04*d

    println("Testing samplep(::AffPro)")
    	
    y = @repeat(SDE.samplep(t, x, T, P0), N)
    mcy = mc(sum(y,1)[:])
	println(repr(mcy), " ", sum(mu0))
    @test abs(mcy[1]- sum(mu0)) < qu*mcy[2]
    muhat = mean(y,2)
	Khat = cov((y.-x)')
	println("K vs. sample covariance ", norm(Khat - K0),  " < ", 0.04*d)
	@test norm(Khat-K0) < 0.04*d
	
	
	#checks if transition density p integrates to 1 by monte carlo integration 
	println("Testing lp")

	function testp(h, x0, N, B, beta, Si)
	    local P, P0
		mu =  B*x0+beta
		P = MvLinPro(B, beta, Si)
		P0 = MvAffPro(mu, Si)
    	# monte carlo integration of the transition density p with proposals according to p0
		w1 = zeros(N)
		for n in 1:N
			y = SDE.samplep(0., x0, h, P0 )
			w1[n] =  exp(SDE.lp(0., x0, h,  y, P)-SDE.lp(0., x0, h, y, P0))
		end
		# monte carlo integration of the transition density p0 with proposals according to p
	    w2 = zeros(N)
		for n in 1:N
			y = SDE.samplep(0., x0, h, P )
			w2[n] =  exp(SDE.lp(0., x0, h,  y, P0)-SDE.lp(0., x0, h, y, P))
		end

		return mc(w1), mc(w2), w1, w2
	end	


	srand(5)
	n1, n2, w1, w2 = testp(0.05, x, N, B, beta, Si)
    print("Testing if transition density p integrates to 1: ")
	println(repr(n1), repr(n2))
	@test abs(n1[1]-1.) < qu*n1[2]
	@test abs(n2[1]-1.) < qu*n2[2]

	@test abs(SDE.lp(t, x,T, v, P) - SDE.varlp(t, x,T, v, (t,s) -> t-s, B, beta, s->A)) < 1E-10


    tmin = t/2
    tmax = T
    T = tmax-tmin
    
	s = SDE.soft(t, tmin, T)
	@test norm(t - SDE.tofs(s, tmin, T)) < 1E-13
	@test norm(t - (tmax - (T - s)^2/T)) < 1E-13

	v1 = SDE.Vs(s, T, v, P0)
	v2 = SDE.V(t, tmax, v, P0)
	@test norm(v1-v2) < d^2*1E-13
	
	v1 = SDE.Vs(s, T, v, P)
	v2 = SDE.V(t, tmax, v, P)
 	@test norm(v1-v2) < d^2*1E-13
	
    v1 = SDE.dotVs(s, T, v, P0)
    v2 = SDE.dotVs(s, T, v, P00)
#    println ("dotv ", repr(v1) ,repr(v2))
    @test norm(v1 - v2) < 1E-5

	us =  SDE.uofx(s, x, T, v,  P)
	xt =  SDE.xofu(s, us, T, v,  P)
	@test norm(x - xt) < 1E-10
	
	
	r1 =  SDE.r(t, x, tmax, v, P)
 	r2 =  SDE.J(s, T, P) *us*T/(T-s)
	@test norm(r1 .- r2) < 1E-10

end

