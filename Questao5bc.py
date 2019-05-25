#Questao para casa 
from numpy import loadtxt, zeros
from time import monotonic

matriz = loadtxt("matriz3.txt", float)
vetor = loadtxt("vetor3.txt", float)
pivo = 0
aux = 0


#Gauss-Jordan
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


t1 = monotonic()

matriz_solucao, vetor_solucao = gauss_jordan(matriz, vetor)

t2 = monotonic()


#Gauss

t3 = monotonic()

A = loadtxt('matriz3.txt',float)
v = loadtxt('vetor3.txt',float)
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

t4 = monotonic()

#Laço
for i in range(0, len(vetor_solucao)):
    for j in range(0, len(vetor_solucao)):
        if(matriz_solucao[i][j] == -0):
            matriz_solucao[i][j] = 0
            

#print("SOLUÇÃO DO SISTEMA POR GAUSS")
#print("----------------------------")
#print("Matriz diagonal:")
#print(A)
#print("")
#print("Vetor solução:")
#print(x)
#print("")
#print("SOLUÇÃO DE SISTEMA POR GAUSS-JORDAN")
#print("-------------------------------------")
#print("Matriz diagonal:")
#print(matriz_solucao)
#print("")
#print("Vetor solução:")
#print(vetor_solucao)
#print("")
print ("TEMPO GASTO EM GAUSS JORDAN", t2-t1)
print("")
print ("TEMPO GASTO EM GAUSS", t4-t3)
print("")
print("RAZÃO ENTRE GAUSS-JORDAN/GAUSS", (t2-t1)/(t4-t3))

"""O algoritmo de gauss-jordan rodou mais rápido"""

""" No método de gauss-jordan temos como resultado final da matriz escalonada
uma matriz identidade, logo não dá pra termos a matriz diagonal com os autovalores"""