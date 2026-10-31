using HomotopyContinuation, LinearAlgebra, Random, JSON
D = JSON.parsefile(ARGS[1]); tag = ARGS[2]; Random.seed!(3)
N = D["N"]; neg = D["neg"] .+ 1; add = [r .+ 1 for r in D["add"]]
met = [m[1]*cis(2pi*m[2]/(2N)) for m in D["met"]]
bic = [cis(2pi*D["bic"][i][j]/N) for i in 1:N, j in 1:N]
zero = findfirst(i->add[i][i]==i && all(add[i][j]==j for j in 1:N), 1:N)
sub(x,y) = add[x][neg[y]]
others = [i for i in 1:N if i!=zero]
half = Int[]; selfi = Int[]
for i in others
  if neg[i]==i; push!(selfi,i) elseif !(neg[i] in half); push!(half,i) end
end
@var e[1:length(half)] s[1:length(selfi)] e0
function Ev(x)
  x==zero && return e0
  i = findfirst(==(x),half); i!==nothing && return e[i]
  j = findfirst(==(x),selfi); j!==nothing && return s[j]
  i = findfirst(==(neg[x]),half); return (1/met[x])/e[i]
end
vars=[e;s;e0]; Dn = prod(e)
Chat = -(1/sqrt(N))*sum(Ev(x) for x in 1:N)
eqs = Expression[]
for x in 1:N, y in 1:N
  push!(eqs, expand((bic[x,y]*Ev(x)*Ev(y) - Chat - (1/sqrt(N))*sum(Ev(sub(y,z))*Ev(z)*met[z]*Ev(sub(x,z)) for z in 1:N))*Dn^3))
end
push!(eqs, e0^2 - sqrt(N)*e0 - 1)
for (j,u) in enumerate(selfi); push!(eqs, s[j]^2*met[u] - 1); end
R = randn(ComplexF64, length(vars), length(eqs))
res = solve(System(R*eqs; variables=vars); show_progress=false, threading=true)
vals(v) = [ComplexF64(to_number(evaluate(Ev(x), vars=>v))) for x in 1:N]
function resid(E)
  C = -sum(E)/sqrt(N)
  maximum(abs(bic[x,y]*E[x]*E[y] - C - sum(E[sub(y,z)]*E[z]*met[z]*E[sub(x,z)] for z in 1:N)/sqrt(N)) for x in 1:N, y in 1:N)
end
U = Vector{Vector{ComplexF64}}()
for v in solutions(res; only_nonsingular=false)
  any(abs.(v[1:length(half)]) .< 1e-8) && continue
  E = vals(v); resid(E) < 1e-8 || continue
  any(norm(E-t)<1e-6 for t in U) || push!(U,E)
end
println("$tag N=$N paths=$(nresults(res)) solutions=$(length(U))")
for (k,E) in enumerate(U)
  realt = all(abs(conj(E[x]) - met[x]*E[x])<1e-7 for x in 1:N)
  unit = all(abs(conj(E[x]) - met[x]*E[neg[x]])<1e-7 for x in 1:N)
  even = all(abs(E[x]-E[neg[x]])<1e-7 for x in 1:N)
  r = findfirst(t->norm(t-[E[neg[x]] for x in 1:N])<1e-6, U)
  println("  sol$k E0=$(round(real(E[zero]),digits=4)) real=$realt unitary=$unit even=$even refl->sol$r")
end
open("fqdG_$(tag).txt","w") do io; for E in U; println(io, join(string.(real.(E))," "), " | ", join(string.(imag.(E))," ")); end; end
