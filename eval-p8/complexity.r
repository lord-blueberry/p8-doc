N = (8192) * (8192)
M = 540 * 4060770

W = 128
I.major = 10
I.clean = 35000

S = 4000
I.cd = 1
J=5


ld.N = log(N, base=2)
ld.2N = log(2*N, base=2)
clean = I.major * 2 * (M + W*(2*N*ld.2N+2*N) + N * ld.N)
clean = clean + I.major * (I.clean * 2 * N)

CD <- function(S,M, I.cd, J) {
  return(S * 7 * M + I.cd * (S * 4*M + J*2*M) + J*(M + 2*N*ld.2N))
}

cd_cycles <- c(1, 5, 10, 15)
results <- c()
for (i in cd_cycles) {
  res <-clean/CD(S, M, i, J)
  results <- c(results, toString(format(round(res, 2), nsmall = 2)))
}


print(paste(round(clean/CD(1000,M, cd_cycles, J), 2), collapse=" & "))
print(paste(round(clean/CD(2000,M, cd_cycles, J), 2), collapse=" & "))
print(paste(round(clean/CD(3000,M, cd_cycles, J), 2), collapse=" & "))
print(paste(round(clean/CD(4000,M, cd_cycles, J), 2), collapse=" & "))

N = (8192*4) * (8192*4)
ld.N = log(N, base=2)
ld.2N = log(2*N, base=2)
clean = I.major * 2 * (M + W*(2*N*ld.2N+2*N) + N * ld.N)
clean = clean + I.major * (I.clean * 2 * N)
print(paste(round(clean/CD(1000,M, cd_cycles, J), 2), collapse=" & "))
print(paste(round(clean/CD(2000,M, cd_cycles, J), 2), collapse=" & "))
print(paste(round(clean/CD(3000,M, cd_cycles, J), 2), collapse=" & "))
print(paste(round(clean/CD(4000,M, cd_cycles, J), 2), collapse=" & "))

print(paste(round(clean/CD(10000,M, cd_cycles, J), 2), collapse=" & "))
print(paste(round(clean/CD(15000,M, cd_cycles, J), 2), collapse=" & "))
print(paste(round(clean/CD(20000,M, cd_cycles, J), 2), collapse=" & "))