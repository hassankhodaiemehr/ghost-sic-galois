using HomotopyContinuation, LinearAlgebra, Random
d = parse(Int, ARGS[1]); Random.seed!(parse(Int, get(ARGS,2,"1")))
ω = cis(2π/d)
X = [mod(i-j-1,d)==0 ? 1.0+0im : 0im for i in 1:d, j in 1:d]   # shift: X e_j = e_{j+1}
Z = Diagonal([ω^k for k in 0:d-1])
Dp(a,b) = X^a * Matrix(Z)^b
@var x[1:d] y[1:d]
ps = [(a,b) for a in 0:d-1 for b in 0:d-1 if (a,b)!=(0,0)]
s = sum(y[i]*x[i] for i in 1:d)
eqs = [ (d+1)*(transpose(y)*Dp(a,b)*x)*(transpose(y)*Dp(a,b)'*x) - s^2 for (a,b) in ps ]
# random square subsystem in P^{d-1} x P^{d-1}
m = 2d-2
R = randn(ComplexF64, m, length(eqs))
F = System(R*eqs; variables=[x;y], variable_groups=[x,y])
res = solve(F; show_progress=false, threading=true)
sols = solutions(res; only_nonsingular=false)
println("d=$d paths=$(nresults(res)) finite sols=$(length(sols))")
# filter by full system, and exclude s=0
full(v) = begin xv=v[1:d]; yv=v[d+1:end]; xv/=norm(xv); yv/=norm(yv)
  sv = transpose(yv)*xv
  abs(sv) < 1e-6 && return Inf
  maximum(abs((d+1)*(transpose(yv)*Dp(a,b)*xv)*(transpose(yv)*Dp(a,b)'*xv)/sv^2 - 1) for (a,b) in ps) end
good = [v for v in sols if full(v) < 1e-6]
# real structure: [conj(y)] == [x]
isreal_(v) = begin xv=v[1:d]/norm(v[1:d]); yv=conj(v[d+1:end]); yv/=norm(yv); abs(abs(dot(xv,yv))-1) < 1e-6 end
nreal = count(isreal_, good)
println("solutions of full complexified system: $(length(good));  sigma-fixed (genuine SIC fiducials): $nreal;  non-real (ghost-type): $(length(good)-nreal)")
sing = count(r->is_singular(r), results(res)); println("singular endpoints: $sing")
