using HomotopyContinuation, LinearAlgebra, Random
N = parse(Int, ARGS[1]); a = parse(Int, get(ARGS,2,"1")); Random.seed!(11)
z2N = cis(pi/N); met(x) = z2N^(a*mod(x,N)*(mod(x,N)+N)); bic(x,y) = met(x+y)/(met(x)*met(y))
half = [u for u in 1:N-1 if u < mod(-u,N)]; selfinv = [u for u in 1:N-1 if u == mod(-u,N)]
@var e[1:length(half)] s[1:length(selfinv)] e0
function Ev(x)
  x = mod(x,N); x==0 && return e0
  i = findfirst(==(x),half); i!==nothing && return e[i]
  j = findfirst(==(x),selfinv); j!==nothing && return s[j]
  i = findfirst(==(mod(-x,N)),half); return (1/met(x))/e[i]   # E(-u)=<u>^{-1}/E(u)
end
vars = [e; s; e0]
# clear denominators: multiply each equation by prod e
D = prod(e)
Chat = -(1/sqrt(N))*sum(Ev(x) for x in 0:N-1)   # C = -Ehat(0)
eqs = Expression[]
for x in 0:N-1, y in 0:N-1
  ex = bic(x,y)*Ev(x)*Ev(y) - Chat - (1/sqrt(N))*sum(Ev(y-z)*Ev(z)*met(z)*Ev(x-z) for z in 0:N-1)
  push!(eqs, expand(ex*D^3))
end
push!(eqs, e0^2 - sqrt(N)*e0 - 1)
for u in selfinv; push!(eqs, s[findfirst(==(u),selfinv)]^2*met(u) - 1); end
R = randn(ComplexF64, length(vars), length(eqs))
res = solve(System(R*eqs; variables=vars); show_progress=false, threading=true)
full(v) = begin
  sub = Dict(zip(vars,v)); vals=[ComplexF64(to_number(evaluate(Ev(x),vars=>v))) for x in 0:N-1]
  C = -sum(vals)/sqrt(N)
  r = maximum(abs(bic(x,y)*vals[mod(x,N)+1]*vals[mod(y,N)+1] - C - sum(vals[mod(y-z,N)+1]*vals[mod(z,N)+1]*met(z)*vals[mod(x-z,N)+1] for z in 0:N-1)/sqrt(N)) for x in 0:N-1, y in 0:N-1)
  (r, vals) end
uniq = Vector{Vector{ComplexF64}}()
for v in solutions(res; only_nonsingular=false)
  any(abs.(v[1:length(half)]) .< 1e-8) && continue
  r, vals = full(v); r < 1e-8 || continue
  any(norm(vals-t)<1e-6 for t in uniq) || push!(uniq, vals)
end
println("N=$N a=$a paths=$(nresults(res)) solutions=$(length(uniq))")
open("fqd_N$(N)_a$(a).txt","w") do io
  for v in uniq; println(io, join(string.(real.(v))," ")," | ",join(string.(imag.(v))," "), " | 0 | 0"); end
end
