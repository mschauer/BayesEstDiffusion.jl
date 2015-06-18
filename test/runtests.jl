using Base.Test
include("../SDE.jl")
using SDE
#is called with arg "travis" from travis
test_type = length(ARGS) == 1 ? ARGS[1] : ""

include("testsde.jl")
include("testlinproc.jl")
include("testtc0.jl")
include("testtc.jl")


