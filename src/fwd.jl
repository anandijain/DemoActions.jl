struct Dual{T}
    x::T # primal ?
    der::T # derivative
end

# """ a and b are functions
# $ c(x) = a(x) + b(x) $
# $ dc/dx = da/dx + db/dx $
# """

Base.:+(a::Dual, b::Dual) = Dual(a.x + b.x, a.der + b.der)

Base.:+(a::Dual, b::Number) = Dual(a.x + b, a.der)
Base.:+(a::Number, b::Dual) = b + a


Base.:-(a::Dual, b::Dual) = Dual(a.x - b.x, a.der - b.der)

Base.:-(a::Dual, b::Number) = Dual(a.x - b, a.der)
Base.:-(a::Number, b::Dual) = Dual(a - b.x, -b.der)

Base.:-(a::Dual) = Dual(-a.x, -a.der)
# Base.:inv(a::Dual) = Dual(1/a.x, -a.der)

Base.:*(a::Dual, b::Dual) = Dual(a.x * b.x, a.der * b.x + a.x * b.der) # product rule

Base.:*(a::Dual, b::Number) = Dual(a.x * b, a.der * b)
Base.:*(a::Number, b::Dual) = b * a 

Base.:/(a::Dual, b::Dual) = Dual(a.x / b.x, (a.der * b.x - a.x * b.der) / (b.x ^ 2))
Base.:/(a::Dual, b::Number) = a * inv(b)
Base.:/(a::Number, b::Dual) = Dual(a / b.x, ( -a * b.der ) / ( b.x ^ 2 ) ) # 

Base.:^(a::Dual, b::Integer) = Base.power_by_squaring(a, b) 
Base.:^(a::Number, b::Dual) = Dual(Float64(a^b.x), a^b.x * log(a))

Base.:≈(a::Dual, b::Dual) = a.x ≈ b.x && a.der ≈ b.der

derivative(f, x) = f(Dual(x, one(x))).der

# primatives
Base.sin(a::Dual) = Dual(sin(a.x), a.der * cos(a.x))
Base.cos(a::Dual) = Dual(cos(a.x), -a.der * sin(a.x))
Base.tan(a::Dual) = sin(a)/cos(a)

Base.csc(a::Dual) = Dual(-csc(a.x)*cot(a.x), -a.der * sin(a.x))
Base.sec(a::Dual) = Dual(tan(a.x), -a.der * sin(a.x))
Base.cot(a::Dual) = 1/tan(a)

Base.exp(a::Dual) = Dual(exp(a.x), exp(a.x) * a.der) 
Base.log(a::Dual) = Dual(log(a.x), a.der/a.x)

# used with ∇ and N-dimensional scalar valued functions
# struct MultiDual{N, T}
#     x::T
#     ders::SVector{N, T}
# end

# Base.:+(a::MultiDual, b::MultiDual) = MultiDual(a.x + b.x, a.ders .+ b.ders)
# Base.:*(a::MultiDual, b::MultiDual) = MultiDual(a.x * b.x, a.ders .* b.x .+ a.x .* b.ders)

# ∇(f, x) = f(MultiDual(x, ones(length(x)))).ders

# struct MultiDual2{M, N, T}
#     x::SVector{M, T}
#     J::SMatrix{M, N, T}
# end