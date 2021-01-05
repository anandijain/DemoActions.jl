""" 
Learning repo for testing various CI and benchmarking tools.

Mainly benchmarking my naive ForwardDiff implementation and Dormand-Prince solver on models I find interesting.

Working through `Mathematical Physiology` by James Keener and `Biomolecular Feedback Systems` by Murray.
"""
module DemoActions

using StaticArrays

include("mysolvers.jl")
include("models.jl")
include("fwd.jl")

# solver interfaces
export MyProbInplace, mysolve, mysolve_inplace

# solver algorithms
export eulers_inplace_f!, 
    midpoint!, 
    midpoint_inplace_f!, 
    rk4!, 
    rk4_inplace_f!, 
    rk4_inplace_f_prealloc_k!, 
    dp5, 
    dp5!, 
    dp5_inplace!

# mini model-zoo
export lotka!, transcriptional_component!, bistable_gene_circuit!, negative_autoregulation!, phosphorylation_cycle!

end
