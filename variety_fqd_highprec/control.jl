using HomotopyContinuation, LinearAlgebra, Random, JSON
# classify a solution x (vector over G with zero index z, negation neg)
function classify(x, pairs)
  tre = all(abs(imag(x[u]+x[w]))<1e-6 for (u,w) in pairs)
  re = all(abs(imag(x[u]))<1e-6 && abs(imag(x[w]))<1e-6 for (u,w) in pairs)
  un = all(abs(abs(x[u])-1)<1e-6 for (u,w) in pairs)
  (tre, re||un)
end
# Control A: random P- and conj-invariant subspace, square system
function controlA(N, mp, mm, e, seed)
  Random.seed!(seed)
  pairs = [(u+1, N-u+1) for u in 1:div(N-1,2)]
  ev = [ (v = zeros(N); v[1]=randn(); for (u,w) in pairs; r=randn(); v[u]=r; v[w]=r; end; v) for _ in 1:mp]
  od = [ (v = zeros(N); for (u,w) in pairs; r=randn(); v[u]=r; v[w]=-r; end; v) for _ in 1:mm]
  B = hcat(ev..., od...); row = B[1:1,:]; cp = row \ [e]; Nl = nullspace(row); r = size(Nl,2)
  @assert r == length(pairs)
  @var t[1:r]; x = B*(cp .+ Nl*t)
  res = solve(System([x[u]*x[w]-1 for (u,w) in pairs]; variables=t); show_progress=false)
  sols = [B*(cp .+ Nl*v) for v in solutions(res)]
  c = [classify(s,pairs) for s in sols]
  (length(sols), count(first,c), count(last,c))
end
# Control B: true Weil eigenspace (x-coords), but only a random square subset of pairing equations
function controlB(D, seed)
  Random.seed!(seed)
  N = D["N"]; neg = D["neg"] .+ 1; add=[r .+ 1 for r in D["add"]]
  met = [m[1]*cis(2pi*m[2]/(2N)) for m in D["met"]]; bic = [cis(2pi*D["bic"][i][j]/N) for i in 1:N, j in 1:N]
  z = findfirst(i->all(add[i][j]==j for j in 1:N), 1:N)
  F = [bic[x,neg[y]]/sqrt(N) for y in 1:N, x in 1:N]; W = Diagonal(1 ./ met)*F; G0 = sum(1 ./ met)/sqrt(N)
  sq = [sqrt(met[u]) for u in 1:N]; for u in 1:N; if neg[u]<u; sq[u]=sq[neg[u]]; end; end
  pairs = unique([minmax(u,neg[u]) for u in 1:N if u!=z && neg[u]!=u])
  ev = eigen(W); out = Tuple{Int,Int,Int}[]
  for k in 0:2
    lam = G0^(1/3)*cis(2pi*k/3); B = ev.vectors[:, abs.(ev.values .- lam) .< 1e-6]
    X = Diagonal(sq)*B; m = size(X,2); m < 2 && continue
    e0 = (sqrt(N)+sqrt(N+4))/2
    row = X[z:z,:]; cp = row \ [e0]; Nl = nullspace(row); r = size(Nl,2)
    r > length(pairs) && continue
    sel = randperm(length(pairs))[1:r]
    @var t[1:r]; x = X*(cp .+ Nl*t)
    res = solve(System([x[pairs[i][1]]*x[pairs[i][2]]-1 for i in sel]; variables=t); show_progress=false)
    sols = [X*(cp .+ Nl*v) for v in solutions(res)]
    c = [classify(s,pairs) for s in sols]
    push!(out, (length(sols), count(first,c), count(last,c)))
  end
  out
end
println("== Control A: random invariant subspace, square system (solutions, t-real, real-or-unimodular)")
for (N,mp,mm) in [(7,2,2),(9,3,2),(11,3,3),(13,4,3)], seed in 1:3
  println("N=$N dims=($mp,$mm) seed=$seed: ", controlA(N,mp,mm,2.5,seed))
end
println("== Control B: Weil eigenspace, random square subset of pairing equations")
idx = JSON.parsefile("cyc/index.json")
for it in idx
  (it["N"] in (11,13,17,19)) || continue
  occursin(r"a=(1|2|3)$", it["cycle"]) || continue
  D = JSON.parsefile("cyc/g$(it["k"]).json")
  for seed in 1:2
    println("N=$(it["N"]) $(it["cycle"]) seed=$seed: ", controlB(D,seed))
  end
end
