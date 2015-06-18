include(Pkg.dir("SDE","src", "Diffusion.jl"))

using Diffusion
srand(8)

@test brown1(3,2,1) == [3.0]
@test length(brown1(3,4,2)) == 2

#test brown1(0,2,5)[end] sim N(0, 2)

#tests with alpha 0.01
r = 2.576
n = 1000

varmu(x, mu) = dot(x-mu, x-mu)/length(x)
var0(x) = dot(x, x)/length(x)

#testing mean and variance of Brownian motion sampled at few points
@test abs(mean([brown1(0,2,5)[end] for i in 1:n])) < r*sqrt(2/n)
chiupper = 1118.95 #upper 0.005 percentile for n = 1000
chilower = 888.56 #lower 0.005 percentile for n = 1000
@test chiupper >1000.*var([brown1(0,2,5)[end] for i in 1:1000])/2 > chilower

@test brown(5., 1, 2, 1) == [5. 5.]'
@test abs(mean(mean([brown([5.,2.], 2, 2, 5)[:,end] for i in 1:n/2]))-3.5) <  r*sqrt(2/n)


#quadratic variance is proportional to time plus sampling bias
qu = [quvar(brown1(0,2,1000)) for j in 1: 1000]
s2 = var(qu)
@test abs(mean(qu)-2)/sqrt(s2)*sqrt(1000) < 3 

@test diff(ito([1., 2., 3.])) == [1., 2., 3.]


# int0^T w_t dw_t = w_T^2/2 - T/2
@test abs((b -> (ito(ydx(b, diff(b)))[end] - (0.5b[end]^2 - 1)))(brown1(0, 2, 10000))) < 0.1
@test abs((dw -> ito((ydx(ito(dw), dw)))[end] - (0.5ito(dw)[end]^2 - 1) )(dW1(2., 10000))) < 0.1
@test abs((dw -> ito((ydx(ito(dw), dw)))[1,end] - (0.5ito(dw)[1,end]^2 - 1) )(dW(2.,1, 10000))) < 0.1
@test abs((dw -> ito(ito(dw), dw)[end] - (0.5ito(dw)[end]^2 - 1) )(dW1(2., 10000))) < 0.1

# E int0^T w_t dw_t = E W_T^2/2 - T/2 = 0
# takes too long to test...
# mean([(dw -> ito((ydx(ito(dw), dw)))[end])(dW1(2., 10000)) for i in 1:100000]) approx 0




#test (roughly) the quadratic variance of Euler approximation
qu = [quvar(euler(0, 1, (t,x)-> -5x, (t,x) -> 1, 1/1000., dW1(1.,1000))) for j in 1: 1000]
s2 = var(qu)
@test abs(mean(qu)-1)/sqrt(s2)*sqrt(1000) < 5.
 

#  B(t) = (1-t) W(t/(1-t)). 

@test abs(mean([bb(0,1,1,3)[2] for i in 1:1000])-0.5)/sqrt(0.25)*sqrt(1000) < r
@test abs(mean([bb(0,1,1,5)[3] for i in 1:1000])-0.5)/sqrt(0.25)*sqrt(1000) < r
@test chilower < 1000*abs(var0([bb(0,1,1,3)[2] for i in 1:1000]-0.5))/.25 < chiupper
@test chilower < 1000*abs(var0([bb(0,1,1,5)[3] for i in 1:1000]-0.5))/.25 < chiupper

@test_approx_eq sum(dWcond(2, 5, 10)) 2

# Covariance of a Brownian bridge from t_1 to t_2 (t_2-t)(s-t_1)/(t_2-t_1)
# here (2-1)(1)/2 = 1/2

@test chilower < 1000*abs(var0([bb(0,1,2,5)[3] for i in 1:1000]-0.5))*2 < chiupper
@test chilower < 1000*abs(var0([cumsum0(dWcond(2, 2, 5))[3] for i in 1:1000]-1))*2. < chiupper

a1 = sum([abs(2 - quvar(ito(aug(dW1(2., 100),1/50,10)))) for j in 1:1000])
a2 = sum([abs(2 - quvar(ito(aug(dW1(2., 100),1/51,10)))) for j in 1:1000])
a3 = sum([abs(2 - quvar(ito(aug(dW1(2., 100),1/49,10)))) for j in 1:1000])

@test a1 < min(a2, a3)
