using HomotopyContinuation, LinearAlgebra, Random
# x-coordinates: x in X (conj- and P-invariant), x_0 = e, x_u x_{-u} = 1 ; test whether x real or |x|=1 (x̄ = Px)
function test(N, mp, mm, e; seed=1)
  Random.seed!(seed)
  pairs = [(u, N-u) for u in 1:div(N-1,2)]           # N odd, indices 0..N-1 -> 1-based +1
  # random real even vectors (support on 0 and pairs symmetric) and odd vectors
  ev = [ (v = zeros(N); v[1]=randn(); for (u,w) in pairs; r=randn(); v[u+1]=r; v[w+1]=r; end; v) for _ in 1:mp]
  od = [ (v = zeros(N); for (u,w) in pairs; r=randn(); v[u+1]=r; v[w+1]=-r; end; v) for _ in 1:mm]
  B = hcat(ev..., od...)
  row = B[1:1,:]; cp = row \ [e]; Nl = nullspace(row); r = size(Nl,2)
  @var t[1:r]
  x = B*(cp .+ Nl*t)
  eqs = [x[u+1]*x[w+1] - 1 for (u,w) in pairs]
  R = randn(ComplexF64, r, length(eqs))
  res = solve(System(R*eqs; variables=t); show_progress=false)
  tot=0; tr=0
  for v in solutions(res; only_nonsingular=false)
    xv = B*(cp .+ Nl*v)
    maximum(abs(xv[u+1]*xv[w+1]-1) for (u,w) in pairs) < 1e-7 || continue
    tot += 1
    isreal = norm(imag.(xv)) < 1e-6
    isunit = all(abs(conj(xv[u+1]) - xv[w+1])<1e-6 for (u,w) in pairs)
    tr += (isreal || isunit)
  end
  (tot, tr)
end
for (N,mp,mm) in [(7,1,1),(7,2,1),(9,2,1),(11,2,2),(13,3,2)], seed in 1:4
  println("random X: N=$N dim+=$mp dim-=$mm seed=$seed  (solutions, totally-real-or-unitary) = ", test(N,mp,mm,2.5; seed=seed))
end
