const libsigma = Pkg.dir("SDE", "deps", "libsigma")


#void fe_mu(double *mu, double *y, int N, int L, int K)
function fe_mu_c(y::Vector{Float64}, L, K)

	n = 2^(L-1)
	mu = zeros(Float64,2n-1 + K)
	N = length(y)
	product = ccall( (:fe_mu_at, libsigma),
                  Void,
                  (Ptr{Float64}, Ptr{Float64}, Int32, Int32),
                  mu, y, N, L)
	return mu;
end


#void fe_Sigma_at(double *S, double *y, int N, double dt, int L)
function fe_Sigma_c(y::Vector{Float64}, dt::Float64, L)
	n = 2^(L-1)
	N = length(y)
	S = zeros(Float64,(2n-1)^2)
	product = ccall( (:fe_Sigma_at,  libsigma),
                  Void,
                  (Ptr{Float64}, Ptr{Float64}, Int32, Float64, Int32),
                  S, y, N, dt, L)
	return reshape(S, 2n-1 , 2n-1)
	
end
function fe_SigmaB1_c(y::Vector{Float64}, dt::Float64, L)
	n = 2^(L-1)
	N = length(y)
	S = zeros(Float64,(2n)^2)
	product = ccall( (:fe_SigmaB1_at,  libsigma),
                  Void,
                  (Ptr{Float64}, Ptr{Float64}, Int32, Float64, Int32),
                  S, y, N, dt, L)
	return reshape(S, 2n , 2n)
	
end

fe_mu = fe_mu_c
fe_SigmaB1 = fe_SigmaB1_c
fe_Sigma = fe_Sigma_c
