using DemoActions
using Test

@testset "DemoActions.jl" begin
    @test 2*5 == 10
end

# """ 
# The `mysolvers_inplace` test set verifies that the higher order solvers 
# have lower error when comparing against the reference DifferentialEquations.jl solve. 

# The following algorithms will be tested:

#     * eulers_inplace_f!
#     * midpoint_inplace_f!
#     * rk4_inplace_f!
#     * dp5! - this should error currently
   
   
# On the following models (found in src/models.jl): 
#     * lotka!


#todo test that the inplace vs allocing return same answers
# """
@testset "mysolvers" begin

    test_algs = [eulers_inplace_f!, midpoint_inplace_f!, rk4_inplace_f!, dp5_inplace!]
    test_models = [lotka!] 
    u0 = [1., 1.]
    p = [1.5, 1.0, 3.0, 1.0]
    tspan = (0., 10.)
    dt = 0.25

    for alg in test_algs 
        @testset "$(alg)_inplace" begin
            for f in test_models 
                @testset "$f" begin
                    @time sol = mysolve_inplace(alg, f, u0, tspan, dt; p=p)
                    @show sol
                    @test !isnothing(sol)
                end 
            end
        end
    end
end