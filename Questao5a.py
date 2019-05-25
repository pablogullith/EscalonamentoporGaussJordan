#Questão para fazer em casa
#Autor: Pablo Gullith
# Bibliotecas
from numpy import loadtxt, zeros

matriz = loadtxt("matriz.txt", float)
vetor = loadtxt("vetor.txt", float)
pivo = 0
aux = 0
#Sistema por gauss_jordan
def gauss_jordan(matriz, vetor):

    tam = len(vetor)
    
    for i in range(0, tam):
        pivo = matriz[i][i]
        matriz[i][:] = matriz[i][:] / pivo
        vetor[i] = vetor[i] / pivo
       
        for j in range(i+1, tam):
            
            aux = matriz[j][i]            
            matriz[j][:] = matriz[j][:] - aux*matriz[i][:]           
            vetor[j] = vetor[j] - aux*vetor[i]
    
    for i in range(tam-1, -1, -1):
        
        for j in range(i-1, -1, -1):
           
            aux = matriz[j][i]
            
            matriz[j][:] = matriz[j][:] - aux*matriz[i][:]
            
            vetor[j] = vetor[j] - aux*vetor[i]

    return matriz, vetor

#Sistema por gauss
A = loadtxt('matriz.txt',float)
v = loadtxt('vetor.txt',float)
x = zeros(len(v),float)

m = len(v)

#Eliminação Gaussiana
for i in range(m):
    p = A[i,i]
    for j in range(i,m):
        A[i,j]=A[i,j]/p
    v[i]=v[i]/p
    for k in range(i+1,m):
        q = -A[k,i]
        for l in range(i,m):
            A[k,l]=A[k,l]+q*A[i,l]
        v[k]=v[k]+q*v[i]

#Substituição retrocedida
for i in range(m-1,-1,-1):
    x[i] = v[i]
    for j in range(i+1,m):
        x[i] = x[i] - A[i,j]*x[j]

matriz_solucao, vetor_solucao = gauss_jordan(matriz, vetor)

#Laço para tornar -0 em 0
for i in range(0, len(vetor_solucao)):
    for j in range(0, len(vetor_solucao)):
        if(matriz_solucao[i][j] == -0):
            matriz_solucao[i][j] = 0

print("SOLUÇÃO DO SISTEMA POR GAUSS")
print("----------------------------")
print("MATRIZ ESCALONADA POR GAUSS:")
print(A)
print("")
print("SOLUÇÃO POR GAUSS:")
print(x)
print("")
print("SOLUÇÃO DE SISTEMA POR GAUSS-JORDAN")
print("-------------------------------------")
print("MATRIZ ESCALONADA POR GAUSS-JORDAN:")
print(matriz_solucao)
print("")
print("SOLUÇÃO POR GAUSS-JORDAN:")
print(vetor_solucao)


