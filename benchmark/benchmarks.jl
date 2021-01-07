using DemoActions
using BenchmarkTools

SUITE = BenchmarkGroup()
SUITE["inplace"] = BenchmarkGroup()


test_algs = [eulers_inplace_f!, midpoint_inplace_f!, rk4_inplace_f!, dp5_inplace!]
test_models = [lotka!] 
u0 = [1., 1.]
p = [1.5, 1.0, 3.0, 1.0]
tspan = (0., 10.)
dt = 0.25
# prob = MyProbInplace(test_algs[1], lotka!, u0, tspan, dt, p)

for alg in test_algs
    for f in test_models
        prob = MyProbInplace(alg, f, u0, tspan, dt, p)
        SUITE["inplace"][string(alg) * "_" * string(f), u0] = @benchmarkable mysolve_inplace($prob)
    end
end


run(SUITE)