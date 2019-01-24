N = (2048) * (2048) 
M = 75 * 4060770

W = 32
I.major = 6
I.clean = 35000

S = 250
I.cd = 10
J=8


ld.N = log(N, base=2)
ld.2N = log(2*N, base=2)
clean = I.major * 2 * (M + W*(2*N*ld.2N+2*N) + N * ld.N)
clean = clean + (I.clean * 2 * N)

CD <- function(S,M, I.cd, J) {
  return(S * 7 * M + I.cd * (S * 4*M + J*2*M) + J*(M + 2*N*ld.2N))
}

print("memory")
M*250*128/8/1024/1024/1024
print("...")

CD(250,M, 10, J)/clean

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