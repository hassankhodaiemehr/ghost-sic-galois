import json,sys,numpy as np
D=json.load(open(sys.argv[1])); N=D['N']
m=lambda i: D['met'][i][0]*np.exp(2j*np.pi*D['met'][i][1]/(2*N))
b=lambda i,j: np.exp(2j*np.pi*D['bic'][i][j]/N)
err=max(abs(m(D['add'][i][j])/(m(i)*m(j))-b(i,j)) for i in range(N) for j in range(N))
print('N',N,'bicharacter consistency err',err, 'even metric',max(abs(m(i)-m(D['neg'][i])) for i in range(N)))
