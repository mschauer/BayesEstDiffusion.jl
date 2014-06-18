module Schauder
export level, levelK, vectoroflevels, number, finger_pm, finger_permute, permutationmatrix, pickup!, pickup, drop!, drop, dropB1, dropB2, hat, fe_transf, fe_transfB2, fe_transfB1

#%  .. currentmodule:: Schauder
#%    
#%  .. _modschauder:
#%
#%  Module Schauder
#%  ---------------
#%
#%  Introduction
#%  ~~~~~~~~~~~~
#%
#%  In the following ``hat(x)`` is the piecewise linear function taking values
#%  values ``(0,0), (0.5,1), (1,0)`` on the interval ``[0,1]`` and ``0`` elsewhere.
#%  
#%  The Schauder basis of level ``L > 0`` in the interval ``[a,b]`` can be defined recursively
#%  from ``n = 2^L-1`` classical finite elements :math:`\psi_i(x)` on the grid
#%  
#%  	``a + (1:n)/(n+1)*(b-a)``.
#%  
#%  
#%  Assume that f is expressed as linear combination 
#%  
#%  	:math:`f(x) = \sum_{i =1}^n c_i \psi_i(x)`
#%  
#%  with 
#%  
#%  	:math:`\psi_{2j-1}(x) = hat(nx-j + 1)` 	for :math:`j = 1 \dots 2^{L-1}` 
#%  
#%  and 
#%  
#%  	:math:`\psi_{2j}(x) = hat(nx-j + 1/2)` 	for :math:`j = 1 \dots 2^{L-1}-1`
#%
#%  Note that these coefficients are easy to find for the finite element basis, just
#%  
#%  .. code-block:: matlab
#%  
#%  	function fe_transf(f, a,b, L)
#%  	    n = 2^L-1 	
#%  	    return map(f, a + (1:n)/(n+1)*(b-a))
#%  	end
#%
	
#%    
#%  Then the coefficients of the same function with respect to the Schauder basis 
#%  
#%  	:math:`f(x) = \sum_{i = 1}^n c_i \phi_i(x)`
#%  
#%  where for ``L = 2``
#%  
#%  	:math:`\phi_2(x) = 2 hat(x)`
#%  
#%  	:math:`\phi_1(x) = \psi_1(x) = hat(2x)`
#%  
#%  	:math:`\phi_3(x) = \psi_3(x) = hat(2x-1)`
#%  
#%  can be computed directly, but also using the recursion
#%  
#%  	:math:`\phi_2(x) = \psi_1(x) + 2\psi_2(x) + \psi_2(x)`
#%  
#%  This can be implemented inplace (see ``pickup()``), and is used throughout. 
#%  

#%  
#% .. code-block:: matlab
#%  
#%  	for l in 1:L-1
#%  		b = sub(c, 2^l:2^l:n)
#%  		b[:] *= 0.5
#%  		a = sub(c, 2^(l-1):2^l:n-2^(l-1))
#%  		a[:] -= b[:]
#%  		a = sub(c, 2^l+2^(l-1):2^l:n)
#%  		a[:] -= b[1:length(a)]
#%  	end
#%  


#%  Reference 
#%  ~~~~~~~~~
#%  
#%  
#%  .. function:: pickup!(x)
#%  
#%  	Inplace computation of 2^L-1 Schauder-Faber coefficients from 
#%  	``2^L-1`` overlapping finite-element coefficients ``x``.
#%  
#%  	-- inverse of ``Schauder.drop``
#%  
#%  	-- L = level(xj)
#%  
#%  .. function:: drop!(x)
#%  
#%  	Inplace computation of 2^L-1 finite element coefficients from 
#%  	2^L-1 Faber schauder coefficients ``x``.
#%  
#%  	-- inverse of ``Schauder.pickup``
#%  
#%  .. function:: finger_permute(x)
#%  
#%  	Reorders vector ``x`` or matrix ``A`` according to the reordering
#%  	of the elements of a Faber-Schauder-basis from
#%  	left to right, from bottom to top.
#%  
#%  .. function:: finger_pm(L, K)
#%  
#%  	Returns the permuation used in ``finger_permute``.
#%  	Memoized reordering of faber schauder elements from low level to high level. The last K elements/rows are left untouched.
#%  

order(a,b) = a > b ? (b,a): (a,b)


function ilogbi(n::Integer)
	Base.exponent(float(n))
end


#%  .. function:: level(x)
#%  
#%  	Gives the no. of levels of the biggest Schauder basis with less then length(x) elements.
#%  		level(x) = ilogb(size(x,1)+1)
#%  
level(x) = ilogbi(size(x,1)+1)

#%  .. function:: level(x, K)
#%  
#%  	Gives the no. of levels ``l`` of the biggest Schauder basis with less then length(x) elements
#%  	and the number of additional elements ``n-2^l+1``.
#%  
function levelK(x) 
	n = size(x,1)
	l = ilogbi(n+1)
	(l, n-2^l+1)
end

#%  .. function:: vectoroflevels(L, K)
#%               
#%  	Gives a vector with the level of the hierarchical elements.
#%  	


function vectoroflevels(L, K)
	n = 2^L
	o = zeros(Int32, n-1+K)
	i = 0
	for l in 0:L-1
	   for j in 2^l:2^l:n-1
		o[j] = l
	   end
	end
	# additional elements (levels) at the end
	for l in n^2:n^2-1+K
		o[l] = l
	end
	o
end

function number(L, K)
	n = 2^L
	o = zeros(Int32, n-1+K)
	i = 0
	for l in 0:L-1
	   for j in 2^l:2^(l+1):n-1
		i+=1
		o[j] = i
	   end
	end
	for l in n:n-1+K
		i+=1
		o[l] = i
	end
	o
end
order_tbl = nothing

#memoized reordering of faber schauder elements from low level to high level
function finger_pm(L, K)
	global order_tbl
	if order_tbl == nothing
		order_tbl = Dict()
	end
	

	if haskey(order_tbl, (L,K))
               pm = order_tbl[(L,K)]
	else
		pm = sortperm(number(L,K))
		assign(order_tbl,pm, (L,K))

	end
	pm
end

# permutationmatrix(sortperm(number(L))
#  pm * pickup_Sigma!(fe_sigma(rand(100), 0, 1, 3),3) * pm'
function permutationmatrix(p)
	d = length(p)
	P = zeros(eltype(p), d,d)
	for i in 1:d
		P[i, p[i]] = 1
	end
	P
end

function finger_permute(x::Vector)
	(L,K) = levelK(x)
	x[finger_pm(L,K)] #undo with ipermute!(x, pm)
end

function finger_permute(Sigma::Matrix)
	(L,K) = levelK(Sigma) 
	pm = permutationmatrix(finger_pm(L,K))
	pm * Sigma * pm'
end




#leaves K elements untouched
function pickup!(xj)
	L = level(xj)
	n = 2^L-1

	for l in 1:L-1
		b = sub(xj, 2^l:2^l:n)
		b[:] *= 0.5
		a = sub(xj, 2^(l-1):2^l:n-2^(l-1))
		a[:] -= b[:]
		a = sub(xj, 2^l+2^(l-1):2^l:n)
		a[:] -= b[1:length(a)]
	end
	xj
end



#leaves K elements untouched
function drop!(xj)
	L = level(xj)
	n = 2^L - 1
 
	for l in L-1:-1:1
		b = sub(xj, 2^l:2^l:n)
		a = sub(xj, 2^(l-1):2^l:n-2^(l-1))
		a[:] += b[:]
		a = sub(xj, 2^l+2^(l-1):2^l:n)
		a[:] += b[1:length(a)]
		b[:] *= 2
		
	end
	xj
end	

compose(f, g) = x -> f(g(x))
drop = compose(drop!,copy)
pickup = compose(pickup!,copy)


function dropB1(yj)
	xj = copy(yj)
	(L, K) = levelK(xj)
	assert(K==1)
	n = 2^L-1
	xj[1:n] .+= xj[end]
	xj
end

function dropB2(yj)
	xj = copy(yj)
	(L, K) = levelK(xj)
	assert(K==2)
	n = 2^L-1

	for i in 1:n
		xj[i] += xj[end]*i/(n+1) + xj[end-1]*(n+1-i)/(n+1) 
	end
	xj = xj[[n+1, 1:n, n+2]] #permute
end

#%  .. function:: hat(x) 
#%               
#%  	Hat function. Piecewise linear functions with values (-inf,0), (0,0),(0.5,1), (1,0), (inf,0).
#%  	-- ``x`` vector or number
function hat(x) 
	max(1. .- abs(2.*x .- 1.),0.)
end




function fe_transf(f, a,b, L)
	n::Float64 = 2^L-1 	
	return map!(f, a .+ [1.:n]/(n+1.)*(b-a))
end

# versions for different basis elements

# phi1(x) = x
# phi2(x) = 1-x
function fe_transfB2(f, a,b, L)
	n::Float64 = 2^L-1 	
		return	[map!(f, a .+ [1.:n]/(n+1.)*(b-a)) .- [1.:n]/(n+1.)*f(a) .- [n:-1.:1.]/(n+1.)*f(b) ,f(a) , f(b)]
end
# phi1(x) = c
# periodic boundary condition
function fe_transfB1(f, a,b, L)
	n::Float64 = 2^L-1 	
	return	[map!(f, a .+ [1:n]/(n+1)*(b-a)) .- 0.5(f(a) + f(b)), 0.5(f(a) + f(b))]
end

include("npbayes.jl")


end

