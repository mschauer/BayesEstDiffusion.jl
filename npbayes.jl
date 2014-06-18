export ex1, fex1, fex2, bayes_drift, test1, levelK, finger_pm, visualize_posterior

#%  Introduction
#%  ~~~~~~~~~~~~
#%  
#%
#%  The procedure is as follows.
#%  Consider the diffusion process :math:`(x_t\colon 0 \le t \le T)` given by
#%  
#%  	:math:`dx_t = b(x_t) dt + dw_t`
#%  
#%  
#%  where the drift ``b`` is expressed as linear combination 
#%  
#%  	:math:`f(x) = \sum_{i =1}^n c_i \phi_i(x)`
#%  
#%  (see :ref:`modschauder`) and 
#%  prior distribution on the coefficients 
#%  
#%  	:math:`c_i \sim N(0,\xi_i)`
#%  
#%  Then the posterior distribution of :math:`b` given observations :math:`x_t` is given by
#%  
#%  	:math:`c_i | x_s \sim N(W^{-1}\mu, W^{-1})`
#%  	:math:`W = \Sigma + (\operatorname{diag}(\xi))^{-1}`, 
#%  
#%  with the nxn-matrix 
#%  
#%  	:math:`\Sigma_{ij} = \int_0^T \phi_i(x_t)\phi_j(x_t) dt`
#%  
#%  and the n-vector
#%  
#%  	:math:`\mu_i = \int_0^T \phi_i(x_t) d x_t`.
#%  
#%  
#%  Using the recursion detailed in :ref:`modschauder`, one rather computes
#%  
#%  	:math:`\Sigma^\prime_{ij} = \int_0^T \psi_i(x_t)\psi_j(x_t) dt`
#%  
#%  and the n-vector
#%  
#%  	:math:`\mu^\prime_i = \int_0^T \psi_i(x_t) d x_t`
#%  
#%  and uses ``pickup_mu!(mu)`` and ``pickup_Sigma!(Sigma)`` to obtain :math:`\mu` and :math:`\Sigma`. 
#%  

#%  
#%  Optional additional basis functions
#%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#%  
#%  One can extend the basis by additional functions, implemented are variants. ``B1`` includes a constant, ``B2`` two linear functions
#%  
#%  ``B1``
#%  	:math:`\phi_1 \dots \phi_n, c`
#%  
#%  ``B2``
#%  	:math:`\phi_1 \dots \phi_n, \max(1-x, 0), \max(x, 0)`
#%  
#%  To compute ``mu``, use
#%  
#%  .. code-block:: matlab
#%  
#%  	mu = pickup_mu!(fe_mu(y,L, 0))
#%  	mu = fe_muB1(mu, y);
#%  
#%  or
#%  
#%  .. code-block:: matlab
#%  
#%  	mu = pickup_mu!(fe_mu(y,L, 0))
#%  	mu = fe_muB2(mu, y);
#%  
#%  



#%  Reference 
#%  ~~~~~~~~~

#%  
#%  Functions taking ``y` without parameter [a,b] expect ``y`` to be shifted into the intervall ``[0,1]``.
#%  



#%  .. function:: pickup_mu!(mu)
#%               
#%  	computes mu from mu'
#%  	


function pickup_mu!(mu)
	L = level(mu)
	n = 2^L-1
	for l in 1:L-1
		for i in 2^l:2^l:n
	 		 mu[i] *= 2
			 mu[i] += mu[i-2^(l-1)] + mu[i+2^(l-1)] 
		end
		 
	end
	mu
end

#%  .. function:: drop_mu!(mu)
#%               
#%  	Computes mu' from mu.
#%  	



function drop_mu!(mu)
	L = level(mu)
	n = 2^L-1
	for l in L-1:-1:1
		for i in 2^l:2^l:n
	 		 
			mu[i] -= mu[i-2^(l-1)] + mu[i+2^(l-1)] 
			mu[i] /= 2
		end
		 
	end
	mu
end

#%  .. function:: pickup_Sigma!(Sigma)
#%               
#%  	Transforms Sigma' into Sigma.
#%  	


function pickup_Sigma!(Sigma)
	L = level(Sigma)
	n = 2^L-1
	for l in 1:L-1
		for i in 2^l:2^l:n
		 Sigma[i, :] *= 2
		 Sigma[:, i] *= 2		 
		 Sigma[i, :] += Sigma[i-2^(l-1), :] + Sigma[i+2^(l-1), :] 
		 Sigma[:, i] += Sigma[:, i-2^(l-1)] + Sigma[:, i+2^(l-1)] 
		end
		 
	end
	Sigma
end


#%  .. function:: drop_Sigma!(Sigma)
#%               
#%  	Transforms Sigma into Sigma'.
#%  	

function drop_Sigma!(Sigma)
	L = level(Sigma)
	n = 2^L-1


	for l in L-1:-1:1
		for i in 2^l:2^l:n
		 Sigma[i, :] -= Sigma[i-2^(l-1), :] + Sigma[i+2^(l-1), :] 
		 Sigma[:, i] -= Sigma[:, i-2^(l-1)] + Sigma[:, i+2^(l-1)] 
		 Sigma[i, :] /= 2
		 Sigma[:, i] /= 2
		end
		 
	end
	Sigma
end


#use after drop_Sigma
function drop_SigmaB2!(Sigma)
	(L,K) = levelK(Sigma)
	assert(K==2)
	n = 2^L-1
	for i in 1:n
		 Sigma[i, :] -= Sigma[end, :]*i/(n+1) + Sigma[end-1, :]*(n+1-i)/(n+1) 
		 Sigma[:, i] -= Sigma[:, end]*i/(n+1)  + Sigma[:, end-1]*(n+1-i)/(n+1) 

	end
	Sigma
end

#use after drop_Sigma
function drop_SigmaB1!(Sigma)
	(L,K) = levelK(Sigma)
	assert(K==1)
	n = 2^L-1
	for i in 1:n
		 Sigma[i, :] -= Sigma[end, :]
		 Sigma[:, i] -= Sigma[:, end]

	end
	Sigma
end


compose(f, g) = x -> f(g(x))

drop_mu = compose(drop_mu!,copy)
drop_Sigma = compose(drop_Sigma!,copy)


#%  .. function:: fe_mu(y, L, K)
#%               
#%  	Computes mu' from the observations `y` using ``2^L-1`` basis elements
#%  	and returns a vector with ``K`` trailing zeros (in case one ones to customize
#%  	the basis.


function fe_mu_jl(y, L, K)
	n = 2^(L-1) #number of even element/elements in lowest level!
	mu = zeros(2n-1 + K)

	dy = [diff(y),0]
	for i in 1:n - 1
		mu[2i-1] = dot(hat(y*n .- (i - 1)),dy)
		mu[2i] = dot(hat(y*n .- (i -.5)),dy)
	end
	mu[2n-1] = dot(hat(y*n.- (n-1)),dy)
	mu
end




#%  .. function:: fe_muB1(mu, y)
#%               
#%  	Append :math:`\mu_{n+1} = \int_0^T \phi_{n+1} d x_t` with :math:`\phi_{n+1} = 1`.
#%  	


function fe_muB1(mu, y)
	[mu, y[end] - y[1]]
end

#%  .. function:: fe_muB2(mu, y)
#%               
#%  	Append :math:`\mu_{n+1} = \int_0^T \phi_{n+1} d x_t` with :math:`\phi_{n+1} =  \max(1-x, 0)`
#%  	and :math:`\mu_{n+2} = \int_0^T \phi_{n+2} d x_t` with :math:`\phi_{n+2} =  \max(x, 0)`
#% 


function fe_muB2(mu,y)
	dy = [diff(y),0]
	[mu, dot(max(1 .- y, 0), dy), dot(max(y, 0), dy)]
end


#%  .. function:: fe_Sigma(y, dt, L)
#%               
#%  	Computes the matrix Sigma' from the observations `y` uniformly spaced at distance ``dt``
#%	using ``2^L-1`` basis elements.
#%




#int(chol(finger_matrix(pickup_Sigma!(fe_Sigma(rand(10000), 0, 1, 5)))) .!= 0)
#memory intensive, using dot
#dont be confused, even wavelets correspond to odd indices as julia indices start at 1
function fe_Sigma_dot(y, dt::Float64, L)
	n = 2^(L-1)
	S = zeros(2n-1, 2n-1)
	yn = y*n
	evnext = hat(yn)
	for i in 1 : n-1
		ev = evnext
		od = hat(yn .- (i - 0.5))
		evnext = hat(yn .-i)
		S[2i-1, 2i-1] = dot(ev,ev)*dt
		S[2i, 2i] = dot(od,od)*dt
		S[2i-1, 2i] = S[2i, 2i-1] = dot(od, ev)*dt
		S[2i, 2i+1] = S[2i+1, 2i] = dot(evnext,od)*dt
	end
	S[2n-1, 2n-1] = dot(evnext, evnext)*dt
	S
end


#wavelets with a constant
function fe_SigmaB1_dot(y, dt::Float64, L)
	K = 1
	n = 2^(L-1)
	S = zeros(2n-1 + K, 2n-1+K)
	yn = y*n
	evnext = hat(yn)
	for i in 1:n - 1
		ev = evnext
		od = hat(yn .- (i - 0.5))
		evnext = hat(yn .- i)
		
		S[2i-1, 2i-1] = dot(ev,ev)*dt
		S[2i, 2i] = dot(od,od)*dt
		S[2i-1, 2i] = S[2i, 2i-1] = dot(od, ev)*dt
		S[2i, 2i+1] = S[2i+1, 2i] = dot(evnext,od)*dt
		S[end, 2i-1] = S[2i-1, end] = sum(ev)*dt
		S[2i, end] = S[end, 2i] = sum(od)*dt
	end
	S[2n-1, 2n-1] = dot(evnext, evnext)*dt
	S[end, 2n-1] = S[2n-1, end] = sum(evnext)*dt
	S[end, end]  = length(y)*dt
	S
end

function fe_SigmaB2_dot(y, dt::Float64, L)
	K = 2
	n = 2^(L-1)
	S = zeros(2n-1 + K, 2n-1+K)
	yn = y*n
	evnext = hat(yn)
	k1 = max(1-y, 0)
	k2 = max(y, 0)

	for i in 1:n - 1
		ev = evnext
		od = hat(yn .- (i - 0.5))
		evnext = hat(yn .- i)
		
		S[2i-1, 2i-1] = dot(ev,ev)*dt
		S[2i, 2i] = dot(od,od)*dt
		S[2i-1, 2i] = S[2i, 2i-1] = dot(od, ev)*dt
		S[2i, 2i+1] = S[2i+1, 2i] = dot(evnext,od)*dt
		S[end, 2i-1] = S[2i-1, end] = dot(ev,k2)*dt
		S[2i, end] = S[end, 2i] = dot(od, k2)*dt
		S[end-1, 2i-1] = S[2i-1, end-1] = dot(ev,k1)*dt
		S[2i, end-1] = S[end-1, 2i] = dot(od, k1)*dt

	end
	S[2n-1, 2n-1] = dot(evnext, evnext)*dt
	S[end, 2n-1] = S[2n-1, end] = dot(k2,evnext)*dt
	S[end, end]  = dot(k2,k2)*dt
	S[end-1, 2n-1] = S[2n-1, end-1] = dot(k1,evnext)*dt
	S[end-1, end-1] = dot(k1,k1)*dt
	S[end, end-1] = S[end-1, end] = dot(k1,k2)*dt
	S
end


#variant accessing y once, fastest if implemented in c

function fe_Sigma_at(y, dt::Float64, L)
	n = 2^(L-1)
	N = length(y)

	S = zeros(2n-1, 2n-1)
	

	for t = 1:N
		yn = y[t]*n
		i = clamp(ceil(yn), 1, n)
		j = clamp(ceil(yn .- 0.5), 1, n-1)
		S[2i-1, 2i-1] += ((hat(yn .- (i - 1)).^2))*dt
		S[2j, 2j] += ((hat(yn .- (j - 0.5)).^2))*dt
		S[2i-1, 2j] = S[2j, 2i-1] += (hat(yn .- (i - 1)).*hat(yn .- (j - 0.5)))*dt
	end
	S
end




fe_mu = fe_mu_jl
fe_Sigma = fe_Sigma_dot
fe_SigmaB1 = fe_SigmaB1_dot
fe_SigmaB2 = fe_SigmaB2_dot

#can libsigma be used?
if find_library(["libsigma"], [Pkg.dir("SDE","deps")]) != ""
	include("sigma.jl")
 
end






# fe_Sigma([0, 0.25, 0.5, 0.625,0], 2, 0,1, 2)
#3x3 Float64 Array:
# 2.0  0.0  0.0
# 0.0  2.5  0.5
# 0.0  0.5  0.5


#expect [xiA, ..., xiB, xi1, ..., xiL] where xiA ... xiB are variances for the trailing elements
#and xi is the variance on level i with 2^(i-1) elements
function expand_levelwise_coeff(xi,K)
	L = length(xi)-K
	rev = 2^L
	xifull = zeros(2^L-1 + K)
	for i in 1:L
		xifull[rev-(2^i-1):rev-2^(i-1)] = xi[K+i]
	end
	xifull[end-K+1:end] = xi[1: K]
	p = finger_pm(L,K)
	pm = permutationmatrix(p)

	ipermute!(xifull, p)
end
	

#%  .. function:: bayes_drift(x, dt, a, b, L, xirem, beta, B)
#%               
#%  	Performs estimation of drift on observations ``x`` in [a,b] spaced at distance ``dt``
#%	using the Schauder basis of level ``L`` and level wise coefficients decaying at rate ``beta``.
#%	A Brownian motion like prior is obtained for beta= 0.5. The ``K`` remaining optional 
#%  	basiselements have variance ``xirem``. 
#%  
#%  	The result is returned as ``[designp coeff se]`` where ``coeff`` are coefficients of finite elements with maximum at the designpoints ``designp`` and standard error ``se``.
#%  
#%  
#%  	
#%	Observations outside [a,b] may influence the result through ``phi_{n+1}, ..., phi_{n+K}``

function bayes_drift(x, dt, a, b, L, xirem, beta, B)
	
	n = 2^L -1
	
	# shift process
	y = (x-a)/(b-a)

	# number of additional elements
	K = 0 
	
	if (B == "B1"); K = 1; 
	elseif (B == "B2"); K = 2;
	end


	xil = [xirem, 1.*2.^(-2.*[L-1:-1:0]    -beta*[1:L])]
	
	
	println("preparations"); tic()
	xi = expand_levelwise_coeff(xil, K)
 	toq()
	println("Xi ", xi')


	println("mu"); tic()
#	mu = pickedup_mu(y, L)
	mu = pickup_mu!(fe_mu(y,L, 0))
	toq()
	
	if(B=="B1"); mu = fe_muB1(mu, y);
	elseif(B=="B2"); mu = fe_muB2(mu, y);
	end

	
	println("mu ", mu)

	println("Sigma pickup"); tic()

	if B=="B1"
		feS = fe_SigmaB1(y,dt, L)
	elseif B == "B2"
		feS = fe_SigmaB2(y,dt, L)
	else 
		feS = fe_Sigma(y,dt, L)
	end

	println("Sigma f.e.", feS)
	toc()
	tic()
	sigma = pickup_Sigma!(feS)
	toc()
	println("Sigma schauder ", sigma)

	println("chol, coeff")
	tic()
	M = chol(sigma + diagm(xi.^(-1)))
	toc()

	tic()
	coeff = M\(M'\mu)
	toc()

	ci2 =  inv(M)*inv(M') 


 	drop_Sigma!(ci2) #check
	drop!(coeff)
	coeff .*= (b-a)

	if B=="B1"
		drop_SigmaB1!(ci2)
		se = diag(ci2).^(.5)
		coeff = dropB1(coeff)
		designp = [(1:n+1)*(b-a)/(n+1)+a]
	elseif B=="B2"
		drop_SigmaB2!(ci2)
		se = diag(ci2).^(.5)
		se = se[[n+1, 1:n, n+2]]
		coeff = dropB2(coeff)
		designp = [(0:n+1)*(b-a)/(n+1)+a]
	else
		se = diag(ci2).^(.5)
		se = 0.5(se + [se[2],se[1:end-1]])
		designp = [(1:n+K)*(b-a)/(n+1)+a]
	end
	


	println("coeff f.e.", coeff)


	[designp coeff se*(b-a)]
end




 
