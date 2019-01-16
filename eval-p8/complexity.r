N = (8192) * (8192)
M = 540 * 4060770

W = 150
I.major = 5
I.clean = 35000

S = 150
I.cd = 10
J=12


ld.N = log(N, base=2)
ld.2N = log(2*N, base=2)
clean = I.major * 2 * (M + W*(2*N*ld.2N+2*N) + N * ld.N)
clean = clean + I.major * (I.clean * 2 * N)

cd = S * 3 * M + I.cd * (S * 8*M + J*2*M)

clean/cd