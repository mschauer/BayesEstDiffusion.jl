using Distributions
using SDE
using Base.Test
srand(8)
# sum([ norm(mean([ sample(MvWiener(d),  linspace(0.,2.,5)).X[1,end] for i in 1:n])) for j in 1:1000] .< r*sqrt(2/n)) 
#tests with alpha 0.01


n = 1000
d= 2

#0.99 quantile vor mean of sum of squares of d-dimensional normals
#r = sqrt(quantile(Distributions.Chisq(d), 0.99))
r = sqrt(9.211) # for d = 2

#  call 'SDE.sample(t,MvWiener(d)).X[:,i]' n times
function Wn(d, t, n, i)
 wn = zeros(d, n)
 for j in 1:n
   wn[:, j] = SDE.sample(t, Wiener(d)).yy[:,i]
 end
 wn
end
  
@test size(SDE.sample(linspace(3.,4.,2), Wiener(d)).yy) == (d,2)
 
# the mean and variance of a Brownian motion at t=2 is 0 and 2
@test norm(mean([ SDE.sample(linspace(0.,2.,5),Wiener(d)).yy[:,end] for i in 1:n])) < r*sqrt(2/n) #scale with std(W_2)/sqrt(n)
#chiupper =  quantile(Distributions.Chisq(n), 0.995) #upper 0.005 percentile  
#chilower = quantile(Distributions.Chisq(n), 0.005) #lower 0.005 percentile  

chiupper = 1119 #n = 1000, 0.995
chilower = 888.5

@test chiupper > n*var( Wn(1, linspace(0.,2.,5), n, 5))/2 > chilower


# check that W(2) has the right covariance matrix
@test (d-> norm(cov( Wn(d, linspace(0.,2.,5), n, 5)') - diagm(2ones(d)))*sqrt(n)/sqrt(d))(50) < 5. # did not figure that out exactly, should fail less then 99 %
@test (d-> norm(cov( Wn(d,[0.,0.1, 0.3, 0.5, 2.], n, 5)') - diagm(2ones(d)))*sqrt(n)/sqrt(d))(50) < 5. # did not figure that out exactly, should fail less then 99 %




function Bn(d, u, v, t, n, i)
 wn = zeros(d, n)
 for j in 1:n
   wn[:, j] = SDE.samplebridge(t, u, t[end], v, Wiener(d)).yy[:,i]
 end
 wn
end
 

# a "deterministic" bridge with only start end endpoint
@test (global t1 = norm(SDE.samplebridge( linspace(1.,4.,2), [1.,2.], 4., [5.,7.], Wiener(2)).yy - [1. 5.; 2. 7.])) < eps(10.)

# Check that the bridge has the right mean and Variance

# Covariance of a Brownian bridge from t_1 to t_2 (t_2-t)(s-t_1)/(t_2-t_1)
# here (3-0.5)*(0.5-0)/(3) = 0.4166666666666667

#cq = quantile(Distributions.Chisq(2), 0.99)
qchisq = 9.21034037197618
@test norm(mean(Bn(2,[1.,2.], [5.,7.],[1.,1.1, 1.3, 1.5, 3.], n, 4),2) -  [1.,2.] - ([5.,7.] - [1.,2.])*.5/2) < sqrt(qchisq)*sqrt(0.416/n)

 # (3-0.5)*(0.5-0.1)/(2.9) = 0.3448275862068966
@test  chilower < n*var(Bn(1,[1.], [5.],[0.1, 0.3, 0.5, 3.], n, 3),2 )[1]/0.345  <  chiupper


