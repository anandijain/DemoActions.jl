# Canonical test models 

"""
    lotka!(du, u, p, t)

lotka-volterra, inplace 
"""
function lotka!(du, u, p, t)
	x, y = u
	a, b, c, d = p
	du[1] = a*x - b*x*y
	du[2] = -c*y + d*x*y
	nothing
end

# Murray Biomolecular Feedback Systems http://www.cds.caltech.edu/~murray/BFSwiki/index.php/Main_Page
# Ch 3 examples

"Example 3.1 (Transcriptional component)"
function transcriptional_component!(du, u, p, t)
    ctrl = 5. # idk
    α, Κ, δ, κ, γ, n = p
    du[1] = α / (1 + (ctrl / Κ)^n) - δ * u[1]
    du[2] = κ * u[1] - γ * u[2]
    nothing
end

"Example 3.2 (Bistable gene circuit) - Murray"
function bistable_gene_circuit!(du, u, p, t)
    a, b, c, d, e, f = p
    du[1] = a / (1 + (u[2] / b)^c) - d * u[1]
    du[2] = e / (1 + (du[1] / f)^c) - d * u[2]
    nothing 
end

"Example 3.4 (Negative autoregulation) - Murray"
function negative_autoregulation!(du, u, p, t)
    du[1] = 1 / (1 + u[2]) - u[1]
    du[2] = u[1] - u[2]
    nothing 
end

""" Example 3.5 (Phosphorylation cycle) - Murray 

Z + X -> Z + X*
Y + X* -> Y + X
X + X* = Xtot
Ytot = c

Notice that in this example there is a fixed amount of Xtot and Ytot.
"""
function phosphorylation_cycle!(out, du, u, p, t)
    Z, X, X2, Y = u
    Xtot, Ytot = p
    out[1] = - 0.04u[1]              + 1e4 * u[2] * u[3] - du[1]
    out[2] = + 0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2]
    out[3] = X + X2 - Xtot
end

