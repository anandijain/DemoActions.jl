"""
	mysolve(alg, f, u0, tspan, dt; p=nothing)

Allocating solve. 
The difference between `mysolve` and `mysolve_inplace` is 
`mysolve` expects `f` to return an array, not mutate one.
"""
function mysolve(alg, f, u0, tspan, dt; p=nothing)
    m = Matrix{eltype(u0)}(undef, length(u0), Int((tspan[2] - tspan[1]) / dt))
    alg(m, f, u0, tspan, dt, p)
end

"""
	mysolve_inplace(alg!, f!, u0, tspan, dt; p=nothing)

Non-allocating solve.
The difference between `mysolve` and `mysolve_inplace` is 
`mysolve` expects `f` to return an array, not mutate one.
"""
function mysolve_inplace(alg!, f!, u0, tspan, dt; p=nothing)
    m = Matrix{eltype(u0)}(undef, length(u0), Int((tspan[2] - tspan[1]) / dt))
    alg!(m, f!, u0, tspan, dt, p)
end	

"""
	eulers_inplace_f!(m, f!, u0, tspan, dt, p)

Non-allocating explicit first order method. Fixed dt.
"""
function eulers_inplace_f!(m, f!, u0, tspan, dt, p)
    du = similar(u0)
    m[:, 1] = 	u0
    for i in 2:size(m, 2)
        t = (i - 1) * dt + tspan[1]
        f!(du, m[:, i - 1],	 p, t)
        m[:, i] = m[:, i - 1] + dt * du
    end
    m
end	

"""
	midpoint!(m, f, u0, tspan, dt, p)

Explicit second order method. Fixed dt.
Allocates for intermediate `k` values.
This means `f` does not return `nothing`.
"""
function midpoint!(m, f, u0, tspan, dt, p)
# du = similar(u0)
    m[:, 1] = u0
    for i in 2:size(m, 2)
        t = (i - 1) * dt + tspan[1]
        k1 = f(m[:, i - 1], p, t)
        k2 = f(m[:, i - 1] + (dt / 2) * k1, p, t + dt / 2)
        m[:, i] = m[:, i - 1] + dt * k2
    end
    m
end

"""
	midpoint_inplace_f!(m, f, u0, tspan, dt, p)

Explicit second order method. Fixed dt.

TODO: have du allocated beforehand, so that this function
has zero allocations.
"""
function midpoint_inplace_f!(m, f!, u0, tspan, dt, p)
    du = similar(u0)
    m[:, 1] = u0
    for i in 2:size(m, 2)
        t = (i - 1) * dt + tspan[1]
        f!(du, m[:, i - 1], p, t)
        f!(du, m[:, i - 1] + (dt / 2) * du, p, t + dt / 2)
        m[:, i] = m[:, i - 1] + dt * du
    end
    m
end

"""
	rk4!(m, f, u0, tspan, dt, p)

Explicit fourth-order Runge-Kutta method. Fixed dt.
`m` is mutated inplace, but `f` is not an in-place function.
"""
function rk4!(m, f, u0, tspan, dt, p)
    m[:, 1] = u0
    for i in 2:size(m, 2)
        u = m[:, i - 1]
        t = (i - 1) * dt + tspan[1]

        k1 = f(u, p, t)
        k2 = f(u + (dt / 2) * k1, p, t + dt / 2)
        k3 = f(u + (dt / 2) * k2, p, t + dt / 2)
        k4 = f(u + dt * k3, p, t + dt)

        m[:, i] = u + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    end
    m
end

"""
	rk4_inplace_f!(m, f!, u0, tspan, dt, p)

Explicit fourth-order Runge-Kutta method. Fixed dt.
`m` is mutated inplace and `f` is an in-place function.
"""
function rk4_inplace_f!(m, f!, u0, tspan, dt, p)
	
	du = similar(u0)
	m[:, 1] = u0
	# k = Matrix{eltype(u0)}(undef, length(u0), 3)
	
	for i in 2:size(m, 2) # what is the better way to write this, prealloc ks and overwrite?
		u = m[:, i-1]
		t = (i-1)*dt + tspan[1]
		f!(du, u, p, t)
		k1 = copy(du)
		f!(du, u + (dt/2)*k1, p, t + dt/2)
		k2 = copy(du)
		f!(du, u + (dt/2)*k2, p, t + dt/2)
		k3 = copy(du)
		f!(du, u + dt*k3, p, t + dt)

		m[:, i] = u + dt*(k1 + 2*k2 + 2*k3 + du)/6
	end
	m
end

"""
	rk4_inplace_f_prealloc_k!(m, f!, u0, tspan, dt, p)

Explicit fourth-order Runge-Kutta method. Fixed dt.
`m` is mutated inplace and `f` is an in-place function.

The difference with the other `rk4` methods is that 
the `k` values are allocated before the loop.

I want to see which of these is fastest.
"""
function rk4_inplace_f_prealloc_k!(m, f!, u0, tspan, dt, p)

    du = similar(u0)
    m[:, 1] = u0
    k = Matrix{eltype(u0)}(undef, length(u0), 3)

    for i in 2:size(m, 2)
        u = m[:, i - 1]
        t = (i - 1) * dt + tspan[1]
        f!(du, u, p, t)
        k[:, 1] = copy(du)
        f!(du, u + (dt / 2) * k[:, 1], p, t + dt / 2)
        k[:, 2] = copy(du)
        f!(du, u + (dt / 2) * k[:, 2], p, t + dt / 2)
        k[:, 3] = copy(du)
        f!(du, u + dt * k[:, 3], p, t + dt)

        m[:, i] = u + dt * (k[:, 1] + 2 * k[:, 2] + 2 * k[:, 3] + du) / 6
    end
    m
end

"Dummy struct for the Dormand-Prince tableau"
struct DP5Cache
    a
    b
    c
end

"""
	dp5_cache() 

allocates the DP5 tableau.
"""
function dp5_cache() 
    dp_a_1 = [0]
    dp_a_2 = [1 / 5]
    dp_a_3 = [3 / 40, 9 / 40]
    dp_a_4 = [44 / 45, -56 / 15, 32 / 9]
    dp_a_5 = [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729]
    dp_a_6 = [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656]
    dp_a_7 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]
    a = [dp_a_1, dp_a_2, dp_a_3, dp_a_4, dp_a_5, dp_a_6, dp_a_7]
    b = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0]
    c = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1]
    DP5Cache(a, b, c)
end

"""
	dp5(f, u0, tspan, dt, p)
	
Uses `dp5!` to solve f and store in `m`.
"""
function dp5(f, u0, tspan, dt, p)
    m = Matrix{eltype(u0)}(undef, length(u0), Int((tspan[2] - tspan[1]) / dt))
    dp5!(m, f, u0, tspan, dt, p)
end

"""
	dp5!(m, f, u0, tspan, dt, p)
	
Explicit Dormand-Prince method. Fixed dt.
Non-allocating.
"""
function dp5!(m, f, u0, tspan, dt, p)
    cache = dp5_cache() # allocates, inefficient
    a, b, c = cache.a, cache.b, cache.c
    S = length(c) # number of stages
    N = length(u0)
    m[:, 1] = u0
    K = Matrix{eltype(u0)}(undef, N, S) 
    for i in 2:size(m, 2)
        u = m[:, i - 1]
        t = (i - 1) * dt + tspan[1]
        K[:, 1] = f(u, p, t)
        for j in 2:S # calculate K vals
            tmp = sum([a[j][k] * K[:, k] for k in 1:j - 1]) 
            K[:, j] = f(u + dt * tmp, p, t + dt * c[j])
        end
        m[:, i] = m[:, i - 1] + dt * sum([b[k] .* K[:, k] for k in 1:S])
    end
    m
end