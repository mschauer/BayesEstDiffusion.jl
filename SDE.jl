#%  .. currentmodule:: SDE
#%    

#%  SDE
#%  ------- 
#%
#%  
#%  

# http://sdejl.readthedocs.org/en/latest/
module SDE
using Cubature
#include("Schauder.jl")
#using Debug
#using Randm
#using Distributions
#using NumericExtensions

import Base.length
import Base.vec
import Base: copy

export b, sigma, a, H, r, p, Bstar, Bcirc, Bsharp, euler, euler!, guidedeuler, guidedeuler!,  llikeliXcirc, samplep, lp, linexact, linll
 
export CTPro, CTPath, UvPath, MvPath, MvPro, UvPro
export UvLinPro, UvAffPro, MvWiener, MvLinPro, MvAffPro, Wiener, Diffusion, UvTime
export diff1, resample!, sample, samplebridge, setv!

export soft, tofs, txofsu, uofx, xofu, XofU, XofU!, UofX, UofX!, eulerU, eulerU!, llikeliU, MvLinProInhomog, Phims

#%  Miscellaneous
#%  ~~~~~~~~~~~~~

#%  .. function:: syl(a, b, c)
#%               
#%      Solves the Sylvester equation ``AX + XB = C``, where ``C`` is symmetric and 
#%      ``A`` and ``-B`` have no common eigenvalues using (inefficient)
#%      algebraic approach via the Kronecker product, see http://en.wikipedia.org/wiki/Sylvester_equation
#%  
function issquare(a::Matrix)
     size(a,2) == size(a,1)
end

vec(x::Number) = x
range(x) = (min(x), max(x))
extendr(R1, mar) = (R1[1] -mar*(R1[2]-R1[1]),R1[2] + mar*(R1[2]-R1[1]))
hrange(x) = extendr(range(x), 1/7)

function syl(a, b, c)
    !(issquare(a) && issquare(b) && issquare(c)) && error("Arguments not square matrices.")

    k = kron(eye(a), a) + kron(b', eye(b))
    xvec=k\vec(c)
    reshape(xvec,size(c))
end

lyap(bt, c) = syl(bt', bt, c)
lyap(b::Float64, mina::Float64) = 0.5mina/b

#%  
#%  Stochastic Processes
#%  ~~~~~~~~~~~~~~~~~~~~
#%

 
diff1(a) = [i > 1 ? a[i] - a[i-1] : a[1] for i=1:length(a) ]
diff0(a) = [i > 1 ? a[i] - a[i-1] : zero(a[1]) for i=1:length(a) ]
eps2 = sqrt(eps())


cumsum!(v::AbstractVector) = Base.cumsum_pairwise(v, v,  zero(v[1]), 1, length(v))

 
function randn!{T,N}(X::StridedArray{T, N})
    n = length(X)
    for i in 1:n
        X[i] = randn()
    end
end


abstract CTPro{Rank} #rank = dim(W) + dim(t) = dim(W) + 1
 

typealias MvPro CTPro{2}
typealias UvPro CTPro{1}

immutable CTPath{Rank}
    tt :: Array{Float64,1}
    yy :: Array{Float64,Rank}
    CTPath(tt, yy) = new(tt, yy)
end       

typealias MvPath CTPath{2}
typealias UvPath CTPath{1}
length(X::CTPath) = length(X.tt)  

getindex(V::UvPath, I) = (V.tt[I], V.yy[I])
getindex(V::MvPath, I) = (V.tt[I], V.yy[:,I])
endof(V::CTPath) = endof(V.tt)

function setindex!(V::MvPath,y, I) 
    V.tt[I],V.yy[:,I] = y 
end
function sub(V::MvPath, I) 
    SubVecProcPath(sub(V.tt,I), sub(V.yy, 1:size(V.yy,1), I)) 
end

  
function setv!(X::MvPath, v)
    X.yy[:, end] = v
    X
end

MvPath(tt::Array{Float64,1}, n::Integer) = MvPath(tt, zeros(n, length(tt)))
UvPath(tt::Array{Float64,1}) = UvPath(tt, zeros(length(tt)))

type Wiener{Rank}  <: CTPro{Rank}
    #dims::NTuple{Rank-1, Int}  #would be nice to do that.
    dims
end

Wiener() =  Wiener{1}(())
Wiener(d::Int) =  Wiener{2}((d,))
Wiener(d1::Int, d2::Int) =  Wiener{3}((d1,d2))

Wiener{Dim}(dims::NTuple{Dim,Int}) = Wiener{Dim+1}(dims)


function resample!(W::CTPath, P::Wiener)
    assert(P.dims == size(W.yy)[1:end-1])
    sz = prod(P.dims)::Int
    for i = 2:length(W.tt)
        rootdt = sqrt(W.tt[i]-W.tt[i-1])
        for j = 1:sz
            W.yy[sz*(i-1) + j] = W.yy[sz*(i-2) + j] + rootdt*randn()
        end
    end
end


function resamplebridge!(W::CTPath, T, v, P::Wiener)
    #white noise
    assert(P.dims == size(W.yy)[1:end-1])
    sz = prod(P.dims)
    TT = T - W.tt[1]

    wtotal = zeros(sz)
    for i = 2:length(W.tt) 
        rootdt = sqrt(W.tt[i]-W.tt[i-1])
        for j = 1:sz
            wtotal[j] +=  W.yy[sz*(i-1) + j] = rootdt*randn()
        end
    end

    # noise between tt[end] and T
    rootdt = sqrt(T-W.tt[end])
    for j = 1:sz
            wtotal[j] +=  rootdt*randn() + (W.yy[j] - v[j])
    end

    # cumsum
    for i = 2:length(W.tt)
        dt =  (W.tt[i]-W.tt[i-1]) /TT
        for j = 1:sz
            W.yy[sz*(i-1) + j] = W.yy[sz*(i-2) + j]  +  W.yy[sz*(i-1) + j] - wtotal[j]*dt
        end
    end
   
end

function samplebridge(tt, u, T, v, P::Wiener)
    yy = zeros(P.dims..., length(tt))
    yy[1:prod(P.dims)] = u
   
    W = CTPath{length(P.dims)+1}(copy(tt), yy)
    resamplebridge!(W, T, v, P)
    W
end


function sample(tt, P::Wiener)
    yy = zeros(P.dims..., length(tt))
    W = CTPath{length(P.dims)+1}(copy(tt), yy)
    resample!(W, P)
    W
end

typealias MvWiener Wiener{2}

type MvAffPro <: MvPro
    mu::Vector{Float64}
    Sigma::Matrix{Float64}
    A::Matrix{Float64}
    Gamma::Matrix{Float64}
    detGamma::Float64
    d::Int    
    dr::Int
    
    function MvAffPro(mu, Sigma) 
         d = length(mu)
         size(Sigma, 1) == d || throw(ArgumentError("The dimensions of mu and Sigma are inconsistent."))
         A = Sigma*Sigma'
         Gamma = inv(A)
         new(mu, Sigma, A, Gamma , det(Gamma), d, size(Sigma, 2))
    end
end


type UvLinPro <: UvPro
    B::Float64
    beta::Float64
    betabyB::Float64
    Sigma::Float64
    lambda::Float64
    d::Int
    function UvLinPro(B, beta, Sigma) 
        (norm(B) > eps2) || throw(ArgumentError("norm(B) < $eps2")) 
        new(B, beta, B\beta, Sigma, -0.5*(Sigma*Sigma)/B, 1)
    end
end

type MvLinPro <: MvPro
    B::Matrix{Float64}
    beta::Vector{Float64}
    betabyB::Vector{Float64}
    Sigma::Matrix{Float64}
    A::Matrix{Float64}
    lambda::Matrix{Float64}
    d::Int
    dp :: Int    
    function MvLinPro(B::Matrix{Float64}, beta::Vector{Float64}, Sigma) 
        d = length(beta)
        size(B,2) == size(B,1) == d || throw(ArgumentError("The dimensions of beta and B are inconsistent."))
        size(Sigma,1) == d || throw(ArgumentError("The dimensions of beta and Sigma are inconsistent."))
        dp = size(Sigma,2)
        (norm(B) > eps2) || throw(ArgumentError("norm(B) < $eps2, use MvAffPro")) 
        A = Sigma*Sigma'
        lambda = lyap(B', -A)
     
    
        new(B, beta, B\beta, Sigma, A, lambda, d, dp)
    end
end


type MvLinProInhomog <: MvPro
    B::Matrix{Float64}
    beta::Vector{Float64}
    betabyB::Vector{Float64}
    Sigma::Matrix{Float64}
    A::Matrix{Float64}
    lambda::Matrix{Float64}
    d::Int
    dp :: Int    
    ph
    function MvLinProInhomog(B::Matrix{Float64}, beta::Vector{Float64}, ph,  Sigma) 
        d = length(beta)
        size(B,2) == size(B,1) == d || throw(ArgumentError("The dimensions of beta and B are inconsistent."))
        size(Sigma,1) == d || throw(ArgumentError("The dimensions of beta and Sigma are inconsistent."))
        dp = size(Sigma,2)
        (norm(B) > eps2) || throw(ArgumentError("norm(B) < $eps2, use MvAffPro")) 
        A = Sigma*Sigma'
        lambda = lyap(B', -A)
     
    
        new(B, beta, B\beta, Sigma, A, lambda, d, dp, ph)
    end
end


type Diffusion{Rank} <: CTPro{Rank}
    b 
    sigma
    dims::Tuple    
end


type UvAffPro <: UvPro
    mu::Float64
    Sigma::Float64
end

type UvTime <: UvPro
    beta
    Beta
    Sigma::Float64
end


typealias AffPro Union(UvAffPro, MvAffPro)

typealias LinPro Union(UvLinPro, MvLinPro)



function a(s, x, P::CTPro)
    sigma(s,x, P)*sigma(s,x, P)'
end


function b(s, x, P::Diffusion)
    P.b(s,x)
end

function sigma(s, x, P::Diffusion)
    P.sigma(s,x)
end



function b(s, x, P::LinPro)
    P.B*x + P.beta
end

function b(s, x, P::AffPro)
    P.mu
end

function b(s, x, P::UvTime)
    P.beta(s)
end

function a(s, x, P::Union(MvLinPro, MvAffPro))
    P.A
end

function a(s, x, P::Union(UvLinPro, UvAffPro, UvTime))
    P.Sigma*P.Sigma'
end


function sigma(s, x, P::LinPro)
    P.Sigma
end

function sigma(s, x, P::AffPro)
    P.Sigma
end


function gamma(P::Union(UvAffPro, UvTime))
    inv(P.Sigma*P.Sigma)    
end 
function gamma(P::MvAffPro)
    P.Gamma
end 


#%  .. function:: mu(t, x, T, P)
#%           
#%      Expectation :math:`E_(t,x)(X_{T})`
#%      
function mu(t, x, T, P::LinPro)
    phi = expm((T-t)*P.B)
    phi*(x + P.betabyB) - P.betabyB
end    

function mu(t, x, T, P::AffPro)
    x + (T-t) * P.mu
end    

function mu(t, x, T, P::UvTime)
    x + (P.Beta(T) - P.Beta(t))
end    

#%  .. function:: K(t, T, P)
#%           
#%      Covariance matrix :math:`Cov(X_{T}-x_t)`
#%      

function K(t, T, P::LinPro)
    phi = expm((T-t)*P.B)
    P.lambda - phi*P.lambda*phi'
end

function K(t, T, P::MvAffPro)
     (T-t)*P.A
end
function K(t, T, P::Union(UvAffPro, UvTime))
     (T-t)*P.sigma^2
end




#%  .. function:: H(t, T, P)
#%           
#%      Negative Hessian of :math:`\log p(t,x; T, v)` as a function of ``x``.
#%      

function Hinv(t, T, P::LinPro)
    phim = expm(-(T-t)*P.B)
    (phim*P.lambda*phim'-P.lambda)
end

Hinv(t, T, P::MvAffPro) = K(t, T, P::MvAffPro)

function Hx(t, T, P::LinPro, x)
     phim = expm(-(T-t)*P.B)
    (phim*P.lambda*phim'-P.lambda)\x
end

function H(t, T, P::LinPro)
     phim = expm(-(T-t)*P.B)
     inv(phim*P.lambda*phim'-P.lambda)
end

function H(t, T, P::AffPro)
    gamma(P)/(T-t)
end
function Hx(t, T, P::AffPro, x)
    gamma(P)/(T-t)*x
end

function H(t, T, P::UvTime)
    gamma(P)/(T-t)
end
function Hx(t, T, P::UvTime, x)
    gamma(P)/(T-t)*x
end


# cholesky factor of H^{-1}, note that x'inv(K)*x =  norm(chol(K, :L)\x)^2

function L(t,T, P::MvPro)
    chol(Hinv(t, T, P), Val{:L})
end


# technical function

function V(t, T, v, P::LinPro)
    phim = expm(-(T-t)*P.B)
    phim*(v + P.betabyB) - P.betabyB  
end

function V(t, T, v, P::MvAffPro)
    return v - (T-t)*P.mu
end

function V(t, T, v, P::UvAffPro)
    return v - (T-t)*P.mu
end

function V(t, T, v, P::UvTime)
    return v - (P.Beta(T)-P.Beta(t))
end



#%  .. function:: r(t, x, T, v, P)
#%           
#%      Returns :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` where
#%      ``p`` is the transition density of the process ``P``.
#%  

function r(t, x, T, v, P)
    Hx(t, T, P, V(t, T, v, P)-x)
end


#%  .. function:: bstar(t, x, T, v, P::MvPro)
#%           
#%      Returns the drift function of a vector linear process bridge which end at time T in point v.
#%      

function bstar(t, x, T, v, P::CTPro)
    b(t, x,  P) + a(t, x, P) * r(t, x, T, v, P)
end    

#%  .. function:: bcirc(t, x, T, v, Pt::Union(MvLinPro, MvAffPro), P::MvPro)
#%           
#%      Drift for guided proposal derived from a vector linear process bridge which end at time T in point v.
#%      

function bcirc(t, x, T, v, Pt::CTPro, P::CTPro)
    b(t,x, P) + a(t,x, P) * r(t, x, T, v, Pt)
end    




#%  .. function:: lp(t, x, T, y, P)
#%           
#%      Returns :math:`log p(t,x; T, y)`, the log transition density of the process ``P``
#%  
function lp(t, x, T, y, P::MvLinPro)
    z = x - V(t, T, y, P)
    l = L(t, T, P)
    -0.5*P.d*log(2pi) - log(prod(diag(chol(K(t, T, P))))) - 0.5*norm(l\z)^2
end
function lp(t, x, T, y, P::UvLinPro)
    z = x - V(t, T, y, P)
    -0.5log(2pi*K(t, T, P)) - 0.5*norm(z)^2*H(t, T, P) 
end



#%  .. function:: samplep(t, x, T, P) 
#%           
#%      Samples from the transition density of the process ``P``.
#%  

#function samplep(t, x, T, P::MvLinPro) 
#    phi = expm((T-t)*P.B)
#    mu = phi*(x + P.betabyB) - P.betabyB 
#    k = P.lambda - phi*P.lambda*phi'
#    l = chol(k)

#    z = randn(length(x))
#    mu + l*z
#end

function samplep(t, x, T, P::MvLinPro) 
    
    m = mu(t, x, T, P)
    k = K(t, T, P)
    l = chol(k,Val{:L})

    z = randn(P.d)
    m + l*z
end


function samplep(t, x, T, P::MvAffPro) 
        z = randn(length(x))
        return x + sqrt(T-t)*P.Sigma*z + (T-t)*P.mu
end
function samplep(t, x, T, P::UvAffPro) 
        z = randn()
        return x + sqrt(T-t)*P.Sigma*z + (T-t)*P.mu
end

#%  .. function:: exact(u, tt, P)
#%           
#%      Simulate process ``P`` starting in `u` on a discrete grid `tt` from its transition probability.
#%  

function exact(u, tt, P::MvPro)
    M = length(tt)
    dt = diff(tt)
    xx = zeros(length(u), M)
    xx[:,1] = u
    for i in 1 : M-1
         xx[:,i+1] = samplep(dt[i], xx[:,i], P) 
    end
    MvPath(tt, xx)
end

#%  .. function:: ll(X,P)
#%           
#%      Compute log likelihood evaluated in `B`, `beta` and Lyapunov matrix `lambda`
#%      for a observed linear process on a discrete grid `dt` from its transition density.
#%  

function ll(X, P::MvPro)
    M = size(X.xx)[end]
    ll = 0.0
    for i in 1 : M-1
        ll += lp(X.tt[i+1]-X.tt[i], X.xx[:,i], X.xx[:,i+1], P) 
    end
    ll
end



#%  .. function:: lp(s, x, t, y, P)
#%           
#%      Returns :math:`log p(t,x; T, y)`, the log transition density
#%  
function lp(s, x, t, y, P::MvAffPro)
      -1/2*P.d*log(2pi*(t-s)) + 0.5*log(det(P.Gamma))  -dot(0.5*(y-x-(t-s)*P.mu), P.Gamma*(y-x-(t-s)*P.mu)/(t-s))
     
end
function lp(s, x, t, y, P::UvAffPro)
      -1/2*log(2pi*(t-s)) - log(abs(P.Sigma))  - 0.5*(y-x-(t-s)*P.mu)*(y-x-(t-s)*P.mu)/(P.Sigma*P.Sigma*(t-s))
end

function varlp(t, x, T, y, ph, B, beta, a)
    z = (x -  varV(t,T, y, ph, B, t -> beta))
    Q = varQ(t, T, ph,  B, a )
    l = chol(Q, Val{:L})
    K =  expm(ph(T,t)*B)*Q*expm(ph(T,t)*B)'
    -1/2*length(x)*log(2pi) -log(prod(diag(chol(K)))) - 0.5*norm(l\z)^2
end




#%  .. function:: euler(u, W::CTPath, P::CTPro)
#%  
#%      Multivariate euler scheme for ``U``, starting in ``u`` using the same time grid as the underlying Wiener process ``W``.
#%      
copy(W::UvPath) = UvPath(copy(W.tt), copy(W.yy))
euler(u, W::UvPath, P::UvPro) = euler!(UvPath(copy(W.tt), copy(W.yy)),u, W, P)
 
function euler!(Y::UvPath, u, W::UvPath, P::UvPro)
    
    N = length(W)
    N != length(Y) && error("Y and W differ in length.")
  
    ww = W.yy
    tt = Y.tt  
    yy = Y.yy
    tt[:] = W.tt
  
    y = u
        
    for i in 1:N-1
        yy[i] = y
        y = y +  b(tt[i],y, P)*(tt[i+1]-tt[i]) + sigma(tt[i],y, P)*(ww[i+1]-ww[i])
    end
    yy[N] = y
    Y
end
guidedeuler(u, W::UvPath, T, v, Pt::UvPro,  P::UvPro) = guidedeuler!(UvPath(copy(W.tt), copy(W.yy)), u, W, T, v, Pt,  P)

function guidedeuler!(Y::UvPath, u, W::UvPath, T, v, Pt::UvPro,  P::UvPro)

    N = length(W)
    N != length(Y) && error("Y and W differ in length.")
  
    ww = W.yy
    tt = Y.tt  
    yy = Y.yy
    tt[:] = W.tt
  
    y = u
        
    for i in 1:N-1
        yy[i] = y
        y = y +  bcirc(tt[i], y, T, v, Pt, P)*(tt[i+1]-tt[i]) + sigma(tt[i],y, P)*(ww[i+1]-ww[i])
    end
    yy[N] = v
    Y
end


function euler(u, W::MvPath, P::MvPro)
    ww = W.yy
    tt = copy(W.tt)

    N = length(tt)
   
    yy = zeros(size(u)..., N)

    y = copy(u)
        
    for i in 1:N-1
        yy[:,i] = y
        y[:] = y .+  b(tt[i],y, P)*(tt[i+1]-tt[i]) .+ sigma(tt[i],y, P)*(ww[:, i+1]-ww[:, i])
    end
    yy[:,N] = y
    MvPath(tt,yy)
end

function guidedeuler!(Y::MvPath, u, W::MvPath, T, v, Pt::MvPro,  P::MvPro)
    ww = W.yy
    tt = Y.tt
    tt[:] = W.tt

    N = length(tt)
   
    yy = Y.yy

    y = copy(u)
        
    for i in 1:N-1
        yy[:,i] = y
        y[:] = y .+  bcirc(tt[i], y, T, v, Pt, P)*(tt[i+1]-tt[i]) .+ sigma(tt[i],y, P)*(ww[:, i+1]-ww[:, i])
    end
    yy[:,N] = v
    Y
end
guidedeuler(u, W::MvPath, T, v, Pt::MvPro,  P::MvPro) = guidedeuler!(MvPath(copy(W.tt), copy(W.yy)), u, W, T, v, Pt,  P)


#%  .. function:: llikeliXcirc(t, T, Xcirc, b, a,  B, beta, lambda)
#%           
#%      Loglikelihood (log weights) of Xcirc with respect to Xstar.
#%  
#%          t, T -- timespan
#%          Xcirc -- bridge proposal (drift Bcirc and diffusion coefficient sigma) 
#%          b, sigma -- diffusion coefficient sigma target
#%          B, beta -- drift b(x) = Bx + beta of Xtilde
#%          lambda -- solution of the lyapunov equation for Xtilde
#%      


function llikeliXcirc(Xcirc::MvPath, Pt::MvPro, P::MvPro)
    tt = Xcirc.tt
    xx = Xcirc.yy
 
    N = length(tt)
    T = tt[N]    
    v = xx[:, N]
    
    som = 0.
    x = similar(v)
    for i in 1:N-1 #skip last value, summing over n-1 elements
      s = tt[i]
      x[:] = xx[:, i]
      R = r(s, x, T, v, Pt) 
      som += (dot(b(s,x, P) - b(s,x, Pt), R) - 0.5 *trace((a(s,x, P) - a(s,x, Pt)) *(H(s, T, Pt) - (R*R')))) * (tt[i+1]-tt[i])
    end
    
    som
end

function llikeliXcirc(Xcirc::UvPath, Pt::UvPro, P::UvPro)
    tt = Xcirc.tt
    xx = Xcirc.yy

    N = length(tt)
    T = tt[N]    
    v = xx[N]
    
    som = 0.
    for i in 1:N-1 #skip last value, summing over n-1 elements
      s = tt[i]
      x = xx[i]
      R = r(s, x, T, v, Pt) 
      som += ((b(s,x, P) - b(s,x, Pt))*R - 0.5 *((a(s,x, P) - a(s,x, Pt)) *(H(s, T, Pt) - (R*R))) ) * (tt[i+1]-tt[i])
    end
    
    som
end


################################################################

#%  .. function:: tofs(s, tmin, T)
#%                soft(t, tmin, T)
#%  
#%      Time change mapping s in [0, T=t_2 - t_1] (U-time) to t in [t_1, t_2] (X-time), and inverse.
#%      

#t = tmin + s*(2 - s/T) = tmin + T - T(1 - s/T)^2 = tmax - (T - s)^2/T
tofs(s, tmin, T) = tmin .+ s.*(2. .- s/T) 
soft(t, tmin, T) = T-sqrt(T*(T + tmin - t))



#%  .. function:: Vs (s, T, v, B, beta)
#%                dotVs (s, T, v, B, beta)
#%  
#%      Time changed V and time changed time derivative of V for generation of U
#%      


function Vs (s, T, v, P::LinPro, phim = expm(-P.B*T*(1. - s/T)^2))
    phim*( v .+ P.betabyB) .-  P.betabyB
end
function dotVs (s, T, v, P::LinPro, phim = expm(-P.B*T*(1. - s/T)^2))
    phim*( P.B*v .+ P.beta) 
end


#  return v - (tmax-t)*P.mu
#  t = tmax - (T - s)^2/T
#  v - (T-s)^2/T*P.mu

function Vs (s, T, v, P::UvAffPro)
    return v - T*(1. - s/T)^2*P.mu
end
function Vs (s, T, v, P::MvAffPro)
    return v - T*(1. - s/T)^2*P.mu
end

function dotVs (s, T, v,  P::MvAffPro)
    P.mu
end
function dotVs (s, T, v,  P::UvAffPro)
    P.mu
end


#%  .. function:: XofU(UU, tmin, T, v, P) 
#%    
#%      U is the scaled and time changed process 
#%      
#%          U(s)= exp(s/2.)*(v(s) - X(tofs(s))) 
#%      
#%      XofU transforms entire process U sampled at time points ss to X at tt.
#%      
    
# 
xofu(s, u, T, v,  P::CTPro) = Vs(s, T, v, P) .- (T-s)*u
xofu(s, u, T, v,  P::CTPro, phim) = Vs(s, T, v, P, phim) .- (T-s)*u

#careful here, s is in U-time
uofx(s, x, T, v,  P::CTPro)  = (Vs(s, T, v, P) .- x)/(T-s)
uofx(s, x, T, v,  P::CTPro, phim)  = (Vs(s, T, v, P, phim) .- x)/(T-s)

txofsu(s,u, tmin, T, v, P::CTPro) = (tofs(s, tmin, T), xofu(s, u, T, v, P))
txofsu(s,u, tmin, T, v, P::MvLinPro, phim) = 
    (tofs(s, tmin, T), Vs(s, T, v, P, phim) .- (T-s)*u)


function XofU{L}(U, tmin, T, v, P::CTPro{L}) 
    X = CTPath{L}(copy(U.tt), copy(U.yy))
    XofU!(X, U, tmin, T, v, P) 
end


function XofU!(X, U, tmin, T, v, P::MvPro) 
    ss = U.tt
    U = U.yy
    for i in 1:length(ss)
        s = ss[i]
        u = U[:, i]
        X.tt[i] = tofs(s, tmin, T)
        X.yy[:, i] = xofu(s, u, T, v, P)
    end
    X
end


function XofU!(X, U, tmin, T, v, P::UvPro) 
    ss = U.tt
    U = U.yy
    for i in 1:length(ss)
        s = ss[i]
        u = U[i]
        X.tt[i] = tofs(s, tmin, T)
        X.yy[i] = xofu(s, u, T, v, P)
    end
    X
end

function UofX{L}(X, tmin, T, v, P::CTPro{L}) 
    U = CTPath{L}(copy(X.tt), copy(X.yy))
    UofX!(U, X, tmin, T, v, P) 
end

function UofX!(U, X, tmin, T, v, P::UvPro) 
    tt = X.tt
    yy = X.yy
    for i in 1:length(tt)
        t = tt[i]
        x = yy[i]
        s = soft(t, tmin, T)
        U.tt[i] = s
        U.yy[i] = uofx(s, x, T, v, P)
    end
    if norm(tmin + T - tt[end]) < eps(T)
        U.yy[end] = 0.0 
    end
    U
end

function UofX!(U, X, tmin, T, v, P::MvPro) 
    tt = X.tt
    yy = X.yy
    for i in 1:length(tt)
        t = tt[i]
        x = yy[:, i]
        U.tt[i] = s = soft(t, tmin, T)
        U.yy[:, i] = uofx(s, x, T, v,P)
    end
    U
end


#helper functions

Phims(s, T, B) = expm(-T*(1. - s/T)^2*B)

function Ju(s,T, P::LinPro, x, phim = expm(-T*(1. - s/T)^2*P.B))
     sl = P.lambda*T/(T-s)^2
    ( phim*sl*phim'-sl)\x
#    LAPACK.sysv!('L', phim*sl*phim'-sl, copy(x))[1]
#     LAPACK.posv!('L', phim*sl*phim'-sl, copy(x))[2]
end
function Ju!(s,T, P::LinPro, j, x, phim = expm(-T*(1. - s/T)^2*P.B))
     j[:] = phim*P.lambda*phim'-P.lambda
     j[:] *= T/(T-s)^2
     j\x
#    LAPACK.sysv!('L', phim*sl*phim'-sl, copy(x))[1]
#     LAPACK.posv!('L', j, x)[2]
end


function Ju(s,T, P::AffPro, x, _ = 0)
    gamma(P)*x
end    

function J(s,T, P::LinPro, phim = expm(-T*(1. - s/T)^2*P.B))
    sl = P.lambda*T/(T-s)^2
    inv( phim*sl*phim'-sl)
end

function J(s,T, P::AffPro, _ = 0)
    gamma(P)
end

function bU(s, u, tmin, T, v, Pt::Union(LinPro, AffPro), P)
    t, x = txofsu(s, u, tmin, T, v, Pt)
    2./T*dotVs(s,T,v, Pt) - 2/T*b(t, x, P) + 1./(T-s)*(u-2.*a(t, x, P)*Ju(s, T, Pt, u) )
end


function llikeliU(U, tmin, T, v, Pt::Union(MvLinPro, MvAffPro), P)
    ss = U.tt
    uu = U.yy
    
    N = length(ss)
    som = 0. 
    for i in 1:N-1
        s = ss[i]
        u = uu[:, i]
        j = J(s, T, Pt)
        ju = j*u
        t, x = txofsu(s, u, tmin, T, v, Pt)

        z1 = 2*dot(b(t, x, P)  - b(t,x, Pt),ju)
        z2 = -1./(T-s)*trace((a(t,x, P) - a(t,x, Pt)) *( j - T*ju*ju' ))
        som += (z1 + z2)*(ss[i+1]-ss[i])
    end
    som
 
end

#  sum((chol(hinv)'\si).^2) = trace(h*a)
#  si'h*x = ((chol(hinv)'\si)'*(chol(hinv)'\x))
      

function llikeliU(U::MvPath, tmin, T, v, Pt::Union(MvLinPro, MvAffPro), P, Phim)
    ss = U.tt
    uu = U.yy
    u = zeros(Pt.d)
    ju = zeros(Pt.d)
    j = zeros(Pt.d, Pt.d)
    ad = zeros(Pt.d, Pt.d)
    jad = zeros(Pt.d, Pt.d)
    N = length(ss)
    som = 0. 
    @inbounds for i in 1:N-1
        s = ss[i]
        u[:] = uu[:, i]
        
        t, x = txofsu(s, u, tmin, T, v, Pt, Phim[i])
        ad[:] = a(t,x, P) - a(t,x, Pt)
        jad[:] = ad
        Ju!(s, T, Pt, j, jad, Phim[i])
        ju[:] = u
        LAPACK.potrs!('L', j, ju)
        
#        jlsi = jli\sigma(t,x, P)
#        jlsit = jli\sigma(t,x, Pt)
#        jlu = jli\u
#        ju = jli'\jlu
        

        z1 = 2*dot(b(t, x, P)  - b(t,x, Pt),ju)
        z2 = -1./(T-s)*(trace( jad) - T*dot(ju,ad*ju))
#        z2 = -1./(T-s)*(vecnorm(jlsi,2)^2-vecnorm(jlsit,2)^2 - T*(norm(jlsi'*jlu)^2 - norm(jlsit'*jlu)^2))
        som += (z1 + z2)*(ss[i+1]-ss[i])
    end
    som
 
end

function llikeliU(U, tmin, T, v, Pt::Union(UvLinPro, UvAffPro), P)
    ss = U.tt
    uu = U.yy
    
    N = length(ss)
    som = 0. 
    for i in 1:N-1
        s = ss[i]
        u = uu[i]
        j = J(s, T, Pt)
        ju = j*u
        t, x = txofsu(s, u, tmin, T, v, Pt)

        z1 = 2.*(b(t, x, P)  - b(t,x, Pt))*ju
        z2 = -1./(T-s)*((a(t,x, P) - a(t,x, Pt)) *( j - T*ju*ju))
        som += (z1 + z2)*(ss[i+1]-ss[i])
    end
    som
 
end

function eulerU!(U::MvPath, ustart, W::MvPath, tmin, T, v, Pt::Union(MvLinPro, MvAffPro),  P::MvPro)
    ww = W.yy
    U.tt[:] = W.tt
    ss, uu = U.tt, U.yy

    N = length(ss)

    u = copy(ustart)
        
    for i in 1:N-1
        uu[:,i] = u
        t, x = txofsu(ss[i], u, tmin, T, v, Pt)
        bU = 2/T*dotVs(ss[i],T,v, Pt) - 2/T*b(t, x, P) +   1/(T-ss[i])*(u - 2.*a(t, x, P)*Ju(ss[i], T, Pt, u) )
        sigmaU = -sqrt(2.0/(T*(T-ss[i])))*sigma(t, x, P)
        u[:] = u .+  bU*(ss[i+1]-ss[i]) .+ sigmaU*(ww[:, i+1]-ww[:, i])
    end
    uu[:,N] = 0u
    U
end

function eulerU!(U::MvPath, ustart, W::MvPath, tmin, T, v, Pt::Union(MvLinPro, MvAffPro),  P::MvPro, Phim)
    ww = W.yy
    U.tt[:] = W.tt
    ss, uu = U.tt, U.yy

    N = length(ss)

    u = copy(ustart)
        
    @inbounds for i in 1:N-1
        uu[:,i] = u
        t, x = txofsu(ss[i], u, tmin, T, v, Pt, Phim[i])
        bU = 2/T*dotVs(ss[i],T,v, Pt, Phim[i]) - 2/T*b(t, x, P) +   1/(T-ss[i])*(u - 2.*a(t, x, P)*Ju(ss[i], T, Pt, u, Phim[i]) )
        sigmaU = -sqrt(2.0/(T*(T-ss[i])))*sigma(t, x, P)
        u[:] = u .+  bU*(ss[i+1]-ss[i]) .+ sigmaU*(ww[:, i+1]-ww[:, i])
    end
    uu[:,N] = 0u
    U
end

function eulerU!(U::UvPath, ustart, W::UvPath, tmin, T, v, Pt::Union(UvLinPro, UvAffPro),  P::UvPro)
    ww = W.yy
    U.tt[:] = W.tt
    ss, uu = U.tt, U.yy

    N = length(ss)

    u = ustart
        
    for i in 1:N-1
        uu[i] = u
        s = ss[i]
        t, x = txofsu(s, u, tmin, T, v, Pt)
        bU = 2/T*dotVs(s,T,v, Pt) - 2/T*b(t, x, P) +   1/(T-s)*(u - 2.*a(t, x, P)*Ju(s, T, Pt, u) )
        sigmaU = -sqrt(2.0/(T*(T-s)))*sigma(t, x, P)
        u += bU*(ss[i+1]-s) + sigmaU*(ww[i+1]-ww[i])
    end
    uu[N] = 0*u
    U
end

function eulerU(ustart, W::MvPath, tmin, T, v, Pt::Union(MvLinPro, MvAffPro),  P::MvPro)
    eulerU!(MvPath(copy(W.tt), zeros(Pt.d, length(W.tt))), ustart, W, tmin, T, v, Pt,  P)
end
function eulerU(ustart, W::MvPath, tmin, T, v, Pt::Union(MvLinPro, MvAffPro),  P::MvPro, Phim)
    eulerU!(MvPath(copy(W.tt), zeros(Pt.d, length(W.tt))), ustart, W, tmin, T, v, Pt,  P, Phim)
end

function eulerU(ustart, W::UvPath, tmin, T, v, Pt::Union(UvLinPro, UvAffPro),  P::UvPro)
    eulerU!(UvPath(copy(W.tt), zeros(length(W.tt))), ustart, W, tmin, T, v, Pt,  P)
end

#%  .. function:: stable(Y, d, ep)
#%           
#%      Return real stable `d`-dim matrix with real eigenvalues smaller than `-ep` parametrized with a vector of length `d*d`, 
#%  
#%  
#%      For maximum likelihood estimation we need to search the maximum over all stable matrices.
#%      These are matrices with eigenvalues with strictly negative real parts.
#%      We obtain a dxd stable matrix as difference of a antisymmetric matrix and a positive definite matrix.
#%  


function stable(Y, d, ep)

    # convert first d*(d+1)/2 values of Y into upper triangular matrix
    # positive definite matrix
    x = zeros(d,d)
    k = 1
    for i in 1:d
        for j in i:d
        x[i,j] = Y[k]
        k = k + 1
        end
    end
    # convert next d*(d+1)/2 -d values of Y into anti symmetric matrix
    y = zeros(d,d)
    for i in 1:d
        for j  in i+1:d
        y[i,j] = Y[k]
        y[j,i] = -y[i, j]
        k = k + 1
        end
    end
    assert(k -1 == d*d == length(Y))
    
    # return stable matrix as a sum of a antisymmetric and a positive definite matrix
    y - x'*x - ep*eye(2) 
end


###
 

#obtaining r via quadrature
function varmu(t, x, T, P)
    function f(s, y)
        y[:] = expm(-s*P.B)*P.beta
    end
    integral = Cubature.hquadrature(length(P.beta), f, 0, T-t; reltol=1E-15, abstol=1E-15, maxevals=05000)[1]
    expm((T-t)*P.B)*(x + integral)
end    


function varr(t, x, T, v, P)
    mu = varmu(t, x, T, P)
    expm((T-t)*P.B')*inv(K(t, T, P))*(v - mu) 
end


# inhomogenous cases

function varH(t, T, B, a::Function)
    d = size(B,1)
    function f(s, y)
        y[:] = vec(expm(-(s-t)*B)*a(s)*expm(-(s-t)*B)')
    end
    Q = reshape(Cubature.hquadrature(d*d, f, t, T; reltol=1E-15, abstol=1E-15, maxevals=5000)[1], d,d) 
    inv(Q)
end    


function varQ(s, T, ph, B, a::Function)
    d = size(B,1)
    
    function f(tau, y)
        y[:] = vec(expm(ph(s,tau)*B)*a(tau)*expm(ph(s,tau)*B)')
    end
    Q = reshape(Cubature.hquadrature(d*d, f, s, T; reltol=1E-15, abstol=1E-15, maxevals=05000)[1], d,d)
    Q
end
 
function varV(s, T, v, ph, B, beta)
    function f(tau, y)
        y[:] = expm(ph(s,tau)*B)*beta(tau)
    end
    expm(-ph(T,s)*B)*v - Cubature.hquadrature(length(v), f, s, T; reltol=1E-15, abstol=1E-15, maxevals=05000)[1]
    
end

    


end
