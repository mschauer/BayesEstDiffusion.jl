BayesEstDiffusion.jl
====================

This repository contains code accompanying the paper

Frank van der Meulen, Moritz Schauer: Bayesian estimation of  discretely observed  multi-dimensional diffusion processes using guided proposals, http://arxiv.org/abs/1406.4704

Abstract: Estimation of parameters of a diffusion based on discrete time observations poses a difficult problem due to the lack of a closed form expression for the likelihood. From a Bayesian
computational perspective it can be casted as a missing data problem where the diffusion bridges in between discrete-time observations are missing. Next, the computational problem can
be dealt with using a Markov-chain Monte-Carlo method known as  data-augmentation.

However, if unknown parameters appear in the diffusion coefficient, direct implementation of data-augmentation results in a Markov chain that is reducible. Furthermore,
data-augmentation requires  efficient sampling of diffusion bridges, which can be difficult, especially in the multidimensional case.

We present a general framework to deal with with these problems that does not rely on discretisation.  The construction generalises previous approaches and sheds light on the
assumptions necessary to make these approaches work. We illustrate our methods  using guided proposals for sampling diffusion bridges. These are Markov processes obtained by adding a
guiding term to the drift of the diffusion.  In a number of examples we give  general guidelines on the construction of these proposals. We introduce a time changing and scaling of
the guided proposal process for stable numerical implementation. Two numerical  examples demonstrate the performance of our methods.


SDE.jl
======
Code to simulate multivariate diffusion and diffusion bridges. 

SDE.jl - contains the module SDE.jl
misc.jl - contains some helper functions.
test - directory with test for SDE.jl


Inno.jl
=======
Implementation of the innovation scheme for the 1-dimensional example, example 6.1.
                             
ProkarNC.jl 
===========
Implementation of the innovation scheme for the prokaryotic auto regulation network 6.2 
using algorithm 1.

ProkarC.jl
==========
Implementation of the innovation scheme for the prokaryotic auto regulation network 6.2 
using algorithm 2.


Addional files
==============
LICENSE
README.md

autoreg50fo.csv - observation for the autoregulation network


Several programs for plotting the pictures
plotatan.R
plotparC.R
plotpar.jl
plotparNC.R
plotparobs.jl

summary* - create a summary of all simulations



