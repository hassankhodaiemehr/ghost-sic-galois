using HomotopyContinuation, LinearAlgebra, Random, JSON
idx = JSON.parsefile("scan/index.json")
function analyze(D, k)
  Random.seed!(5)
  N = D["N"]; neg = D["neg"] .+ 1; add = [r .+ 1 for r in D["add"]]
  met = [m[1]*cis(2pi*m[2]/(2N)) for m in D["met"]]
  bic = [cis(2pi*D["bic"][i][j]/N) for i in 1:N, j in 1:N]
  zero = findfirst(i->all(add[i][j]==j for j in 1:N), 1:N)
  sub(x,y) = add[x][neg[y]]
  F = [bic[x,neg[y]]/sqrt(N) for y in 1:N, x in 1:N]
  W = Diagonal(1 ./ met) * F
  G0 = sum(1 ./ met)/sqrt(N)
  resid(E) = (C = -sum(E)/sqrt(N); maximum(abs(bic[x,y]*E[x]*E[y] - C - sum(E[sub(y,z)]*E[z]*met[z]*E[sub(x,z)] for z in 1:N)/sqrt(N)) for x in 1:N, y in 1:N))
  U = Vector{Vector{ComplexF64}}()
  pairs = unique([minmax(u,neg[u]) for u in 1:N if u!=zero])
  ev = eigen(W)
  for k in 0:2, sgn in (1,-1)
    lam = G0^(1/3)*cis(2pi*k/3)
    B = ev.vectors[:, abs.(ev.values .- lam) .< 1e-6]; m = size(B,2); m == 0 && continue
    e0 = (sqrt(N) + sgn*sqrt(N+4))/2
    row = B[zero:zero,:]; cp = row \ [e0]; Nl = nullspace(row); r = size(Nl,2)
    if r == 0; E = ComplexF64.(B*cp); resid(E)<1e-8 && push!(U,E); continue; end
    @var t[1:r]
    Ef = B*(cp .+ Nl*t)
    eqs = [Ef[u]*Ef[v]*met[u] - 1 for (u,v) in pairs]
    R = randn(ComplexF64, r, length(eqs))
    res = solve(System(R*eqs; variables=t); show_progress=false)
    for v in solutions(res; only_nonsingular=false)
      E = B*(cp .+ Nl*v); resid(E) < 1e-7 || continue
      any(norm(E-w)<1e-6 for w in U) || push!(U,E)
    end
  end
  stats = Dict{String,Int}()
  for E in U
    realt = all(abs(conj(E[x]) - met[x]*E[x])<1e-6 for x in 1:N)
    unit = all(abs(conj(E[x]) - met[x]*E[neg[x]])<1e-6 for x in 1:N)
    even = all(abs(E[x]-E[neg[x]])<1e-6 for x in 1:N)
    key = (real(E[zero])>0 ? "+" : "-") * (realt ? "R" : "") * (unit ? "U" : "") * (!(realt||unit) ? "X" : "") * (even ? "e" : "")
    stats[key] = get(stats,key,0)+1
  end
  open("sol_$(k).txt","w") do io; println(io, zero); for E in U; real(E[zero])>0 || continue; println(io, join(string.(real.(E))," "), " | ", join(string.(imag.(E))," ")); end; end
  stats
end
for it in idx
  k = it["k"]; D = JSON.parsefile("scan/g$(k).json")
  t = @elapsed s = analyze(D, k)
  println("N=$(it["N"]) cycle=$(it["cycle"]) ", join(["$k:$v" for (k,v) in sort(collect(s))], " "), "  ($(round(t,digits=1))s)")
  flush(stdout)
end
