using HomotopyContinuation, LinearAlgebra, Random, JSON
idx = JSON.parsefile("cyc/index.json")
function analyze(D)
  Random.seed!(5)
  N = D["N"]; neg = D["neg"] .+ 1; add = [r .+ 1 for r in D["add"]]
  met = [m[1]*cis(2pi*m[2]/(2N)) for m in D["met"]]
  bic = [cis(2pi*D["bic"][i][j]/N) for i in 1:N, j in 1:N]
  zero = findfirst(i->all(add[i][j]==j for j in 1:N), 1:N)
  sub(x,y) = add[x][neg[y]]
  F = [bic[x,neg[y]]/sqrt(N) for y in 1:N, x in 1:N]
  W = Diagonal(1 ./ met) * F
  G0 = sum(1 ./ met)/sqrt(N)
  pres(E) = maximum(abs(E[u]*E[neg[u]]*met[u]-1) for u in 1:N if u!=zero)
  resid(E) = (C = -sum(E)/sqrt(N); maximum(abs(bic[x,y]*E[x]*E[y] - C - sum(E[sub(y,z)]*E[z]*met[z]*E[sub(x,z)] for z in 1:N)/sqrt(N)) for x in 1:N, y in 1:N))
  # consistent square roots of met with sqrt(met(u)) = sqrt(met(-u))
  sq = [sqrt(met[u]) for u in 1:N]
  for u in 1:N; if neg[u] < u; sq[u] = sq[neg[u]]; end; end
  U = Vector{Vector{ComplexF64}}()
  pairs = unique([minmax(u,neg[u]) for u in 1:N if u!=zero])
  ev = eigen(W); dims=Int[]
  for k in 0:2, sgn in (1,-1)
    lam = G0^(1/3)*cis(2pi*k/3)
    B = ev.vectors[:, abs.(ev.values .- lam) .< 1e-6]; m = size(B,2); sgn==1 && push!(dims,m); m == 0 && continue
    e0 = (sqrt(N) + sgn*sqrt(N+4))/2
    row = B[zero:zero,:]; cp = row \ [e0]; Nl = nullspace(row); r = size(Nl,2)
    if r == 0; E = ComplexF64.(B*cp); pres(E)<1e-8 && push!(U,E); continue; end
    @var t[1:r]
    Ef = B*(cp .+ Nl*t)
    eqs = [Ef[u]*Ef[v]*met[u] - 1 for (u,v) in pairs]
    R = randn(ComplexF64, r, length(eqs))
    res = solve(System(R*eqs; variables=t); show_progress=false)
    for v in solutions(res; only_nonsingular=false)
      E = B*(cp .+ Nl*v); pres(E) < 1e-7 || continue
      any(norm(E-w)<1e-6 for w in U) || push!(U,E)
    end
  end
  npair = length(U); npent = 0; treal = 0; tbad = 0; pentbad=0; mixed=0
  for E in U
    ok = resid(E) < 1e-7; npent += ok
    others = [u for u in 1:N if u!=zero]
    tv = [sq[u]*(E[u]+E[neg[u]]) for u in others]
    isr = all(abs(imag(x))<1e-6 for x in tv)
    treal += isr
    if !isr && ok; pentbad += 1; end
    rt = all(abs(conj(E[x]) - met[x]*E[x])<1e-6 for x in 1:N); ut = all(abs(conj(E[x]) - met[x]*E[neg[x]])<1e-6 for x in 1:N); mixed += !(rt||ut)
  end
  (dims=dims, npair=npair, npent=npent, treal=treal, pent_nonreal=pentbad, mixed=mixed)
end
for it in idx
  nothing=0
  k = it["k"]; D = JSON.parsefile("cyc/g$(k).json")
  s = analyze(D)
  println("N=$(it["N"]) cycle=$(it["cycle"]) eigdims=$(s.dims) pairing-solutions=$(s.npair) pentagon-solutions=$(s.npent) t-real=$(s.treal) pentagon-but-t-nonreal=$(s.pent_nonreal) mixed-type(pairing)=$(s.mixed)")
  flush(stdout)
end
