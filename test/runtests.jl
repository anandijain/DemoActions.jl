using DemoActions
using Test

@testset "DemoActions.jl" begin
    @test 2*5 == 10
end

""" 
The `mysolvers_inplace` test set verifies that the higher order solvers 
have lower error when comparing against the reference DifferentialEquations.jl solve. 

The following algorithms will be tested:

    * eulers_inplace_f!
    * midpoint_inplace_f!
    * rk4_inplace_f!
    * dp5! - this should error currently
   
   
On the following models (found in src/models.jl): 
    * lotka!
"""
@testset "mysolvers" begin

    test_algs = [eulers_inplace_f!, midpoint_inplace_f!, rk4_inplace_f!, dp5!]
    test_models = [lotka!] 

    for alg in test_algs 
    @testset "$(alg)" begin
        for f in test_models 
            @testset "$f" begin
                @test !isnothing(mysolve_inplace(alg, f, u0, tspan, dt; p=p))
            end 
        end
    end
    # @testset "allocating" begin
        
    #     @test 3*4 == 12    
    # end

end



        # lv_euler = mysolve_inplace(eulers_inplace_f!, lotka!, u0, tspan, dt; p=p)
        # lv_mid = mysolve_inplace(midpoint_inplace_f!, lotka!, u0, tspan, dt; p=p)
        # lv_rk4_inp = mysolve_inplace(rk4_inplace_f!, lotka!, u0, tspan, dt; p=p)
        # # fix to use mysolve_inplace
        # lv_dp5 = mysolve(dp5!, lotka, u0, tspan, dt; p=p)
        
        # # lv_prob = ODEProblem(lotka!, u0, tspan, p)
        # # lv_de = solve(lv_prob, saveat=ts)
        # # mysols = [lv_euler, lv_mid, lv_rk4_inp, lv_dp5]    
        # # myerrs = sum.(abs2, reference_soln .- [lv_eul_err, lv_mid_err, lv_rk4_err, lv_dp5_err])
        
        # # tests that error is decreasing as solver order is increasing
        # @test issorted(myerrs, rev=true)