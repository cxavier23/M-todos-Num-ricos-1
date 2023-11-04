##author: Christopher Xavier Sanchez Duran
### //UAMC

import math as M
import numpy as NP 
from numpy import linalg as LA #calcular espectro y transponer

NP.set_printoptions(precision = 4) #impresion con decimales
Nmax = 30 #max de iteraciones

#Todo lo relacionado con el metodo de Cholesky

def MetodoCholesky(A, b):

    if is_PosDef(A): 
        L = CholeskyFactoring(A)
        y = forward_sub(L, b)
        x = backward_sub(NP.transpose(L), y)

        print("\n\nL :\n", NP.array(L))
        print("\n\nU :\n", NP.transpose(L))

        print("\n\nSolución:\n\tX = ", NP.array(x))

def CholeskyFactoring(A): #usando el algoritmo de Cholesky factoriza a A
    
    n = len(A)
    L = []

    for i in range(n):

        L.append([])

        for j in range(i + 1):
            S = 0

            for k in range(j):
                S += L[i][k] * L[j][k]

            if i > j:
                L[i].append((A[i][j] - S) / L[j][j])
            else:  # if i == j
                L[i].append(M.sqrt(A[j][j] - S))

        for j in range(i + 1, n):
            L[i].append(0)

    return L

def forward_sub(L, b): #Sustitucion hacia adelante
    
    n = len(b)
    y = []

    for i in range(n):
        S = 0
        for k in range(i):
            S += L[i][k] * y[k]

        y.append((b[i] - S) / L[i][i])

    return y

def backward_sub(U, y): #Sustitucion hacia atras

    n = len(y)
    x = []

    for i in range(n):
        S = 0
        for k in range(i):
            S += U[n - 1 - i][n - 1 - k] * x[k]

        x.append((y[n - 1 - i] - S) / U[n - 1 - i][n - 1 - i])

    x.reverse() #invierte el vector porque los elementos se insertaron al reves

    return x



#Todo lo relacionado con el metodo de Jacobi

def MetodoJacobi(A, b, X0, tol): #x(k) = inv(D)*(b-(L+U)x(k-1))
    
    DUL = split_Matrix_Jacobi(A) #DUL[0]=Inv(D), DUL[1]=U, DUL[2]=L
    Db = MtzVecMult(DUL[0], b)
    DLU = MMult(DUL[0], MSum(DUL[1], DUL[2]))

    err = 2 * tol #asegura que se ejecute la primera iteracion
    i = 1
    
    print("\nIteracion\t\t\tAproximación\t\t\t\tError\n")
    print(0, "\t\t", NP.array(X0)) #imprime aprox inicial
  
    while err > tol and i < Nmax:

      X = VRest(Db, MtzVecMult(DLU, X0)) #calcula nueva iteracion
      err = error_Rel(X, X0) #calcula error relativo
      X0 = X
      print(i, "\t\t", NP.array(X), "\t\t\t\t{:.6f}".format(err))
      i += 1

    if i == Nmax: #si no hubo convergencia
      print("\n\n\nNo se alcanzo la convergencia en ", Nmax, " iteraciones.\n\n\n")
      return

    print("\n\nSolución aproximada:\n\t X = ", NP.array(X))
    print("\nCon un error relativo:\n\t Err = {:.6f}".format(err))

def split_Matrix_Jacobi(A): #Parte la matriz del sistema para el metodo

    D_inv = []  # D es la inversa de la diagonal de A
    U = []
    L = []
    n = len(A)

    for i in range(n):

        D_inv.append([])
        U.append([])
        L.append([])

        for j in range(n):

            if i == j:
                D_inv[i].append(1 / A[i][j])
                U[i].append(0)
                L[i].append(0)

            elif i > j:
                D_inv[i].append(0)
                U[i].append(0)
                L[i].append(A[i][j])

            else:  # if i < j
                D_inv[i].append(0)
                U[i].append(A[i][j])
                L[i].append(0)

    return [D_inv, U, L]

def dominant(A,b): #Encuentra una matriz diagonalmente dominante para el sistema usando permutaciones de renglones

    PA = []
    Pb = []
    n = len(A)

    for i in range(n):
        PA.append([])
        Pb.append(0)
        for j in range(n):
            PA[i].append(0)

    for i in range(n):
      S = sum([abs(a) for a in A[i]]) #Criterio para ser diagonalmente dominante
      for j in range(n):
          if 2*abs(A[i][j]) >= S:
              PA[j] = A[i] #Intercambio de renglones
              Pb[j] = b[i] #Permuta al vector constante tambien

    #print(NP.array(PA))
    #print(NP.array(Pb))

    return [PA,Pb]

def error_Rel(X, X0): #Error relativo con norma infinito

    n = len(X)
    return max([abs(X[i] - X0[i]) for i in range(n)]) / max([abs(X[i]) for i in range(n)])



#Validaciones para matrices

def is_Square(A): #checa si es cuadrada
    if len(A) == len(A[0]):
        return True

    return False

def is_Symmetric(A): #checa si es simetrica
    n = len(A)

    if not is_Square(A):
        return False

    for i in range(1, n):
        for j in range(i):
            if A[i][j] != A[j][i]:
                return False

    return True

def is_PosDef(A): #checa si es positivo-definida

    if not is_Symmetric(A):
        return False

    L = LA.eigvals(A)
    for l in L:
        if l <= 0:
            return False

    return True



#Operaciones auxiliares con vectores y matrices

def MSum(A, B): #Suma de matrices

    m = len(A)
    n = len(A[0])
    C = []

    if len(B) == m and len(B[0]) == n:
        for i in range(m):
            C.append([])
            for j in range(n):
                C[i].append(A[i][j] + B[i][j])

    return C

def VRest(v, u): #Diferencia de vectores

    n = len(v)
    w = []

    if len(u) == n:
      for i in range(n):
        w.append(v[i] - u[i])

      return w

def MMult(A, B): #Producto de Matrices

    n = len(A[0]) #A es mxn, B es nxp

    if len(B) == n: #Checa compatibilidad

        C = []
        m = len(A)
        p = len(B)

        for i in range(m):
            C.append([])
            for j in range(p):
                C[i].append(0)
                for k in range(n):
                    C[i][j] += A[i][k] * B[k][j]

        return C #C es mxp

def MtzVecMult(A, v): #Producto matriz-vector
    
    #Nota: se define aparte de MMult porque MMult usa listas bidimensionales
    n = len(A[0])

    if len(v) == n:

        w = []
        m = len(A)

        for i in range(m):
            w.append(0)
            for k in range(n):
                w[i] += A[i][k] * v[k]

    return w



#Para facilitar la captura de datos

def read_Matrix(): #Lee Matrices

    M = []
    n = int(input("\nNum de variables: "))

    for i in range(n):
        M.append([])
        for j in range(n):
            mij = float(input("\nElemento {},{} : ".format(i + 1, j + 1)))
            M[i].append(mij)

    return M

def read_Vector(n): #Lee vectores
    
    v = []

    for i in range(n):
        v.append(float(input("\nEntrada {} : ".format(i + 1))))

    return v



#Interacciones con el usuario

def bienvenida(): #Inicia el programa

    print("\nEste programa fue creado por:\n",
          "\t- Christopher Xavier Sanchez Duran\n", )

    print("Este programa implementa los metodos de Jacobi y Cholesky para encontrar        la solucion de sistemas de ecuaciones lineales.\n\n\n")

def selectMet(): #Selecciona el metodo (Jacobi o Cholesky)

    while True: #Valida la seleccion del usuario
      print("\nSeleccione el metodo a emplear:\n\n",
            " (1)\tMetodo de Jacobi\n",
            " (2)\tMetodo de Cholesky\n",
            " (3)\tSalir\n")

      met = int(input("\nInserte una opcion:\t"))

      if met == 1 or met == 2 or met == 3: #Para salir del ciclo
        return met

      print("\n\nERROR: Opcion invalida. Intente de nuevo.\n\n")



def menuJacobi(): #Seleccionar sistema a resolver usando Jacobi

    A = []
    b = []
    preg = 0

    while True:

        print("\nSeleccione el sistema a resolver:\n")

        print(" (1)")
        printSistemaJacobi(1)
        print("")

        print(" (2)")
        printSistemaJacobi(2)
        print("")

        print(" (3)\tRegresar al menu anterior.\n\n")

        preg = int(input("Inserte una opcion:\t"))

        if preg == 1 or preg == 2:
            break
        elif preg == 3:
            return

        print("\n\nERROR: Opcion invalida. Intente de nuevo.\n\n")

    if preg == 1:

        A = [ [7.3891, 14.778,1],
              [ 54.598,218.39, 2],
              [ 403.34, 2420.6, 3]]
               
        b = [ 6.7781, 201.2, 797.86]

    else: #if preg == 2:

        A = [ [ 2, 3, 1],
              [ 4, 1,-1],
              [ 1,-1,-5] ]

        b = [ 0, 2, 1]

    print("\n\nIntroduzca una aproximacion inicial:\t")
    X0 = read_Vector(len(b)) #Lee la aproximacion inicial

    tol = float(input("\n\nInserte el valor de la tolerancia:\t"))

    Dom = dominant(A,b) #Intenta dar una matriz diagonal dominante para el sistema
    A = Dom[0]
    b = Dom[1]

    print("\nMetodo:\t Jacobi")
    print("\nSistema:")
    printSistemaJacobi(preg)
    print("Tolerancia:\t", tol)

    MetodoJacobi(A, b, X0, abs(tol)) #abs(tol) ignora el signo de la tolerancia

def printSistemaJacobi(preg): #Imprime los sistemas disponibles para Jacobi

  if preg == 1:
    print("\t-2X\t+5Y\t-1Z\t+0W\t=  2\n",
          "\t+4X\t-2Y\t+0Z\t+0W\t=  0\n"
          "\t+0X\t-1Y\t+4Z\t+2W\t=  3\n",
          "\t+0X\t+0Y\t+2Z\t+3W\t= -2\n")

  elif preg == 2:
    print("\t+2X\t+\t3Y\t+1Z\t=\t  0\n",
          "\t+4X\t+\t1Y\t-1Z\t=\t  2\n",
          "\t+1X\t-\t1Y\t-5Z\t=\t  1\n")



def menuCholesky(): #Seleccionar sistema a resolver usando Cholesky

    A = []
    b = []

    while True:

        print("\nSeleccione el sistema a resolver:\n")

        print(" (1)")
        printSistemaCholesky(1)
        print("")

        print(" (2) \n") #Sistema propuesto por Christopher Xavier Sanchez Duran
        printSistemaCholesky(2)
        print("")

        print(" (3) Introducir un sistema\n\n") #Si se desea introducir un sistema manualmente 

        print(" (4)\tRegresar al menu anterior.\n")

        preg = int(input("Inserte una opcion:\t"))

        if preg > 0 and preg < 4:
            break
        elif preg == 4:
            return

        print("\n\nERROR: Opcion invalida. Intente de nuevo.\n\n")

    if preg == 1:
        A = [ [ 9,-3, 2, 1],
              [-3, 5, 0,-1],
              [ 2, 0, 4, 1],
              [ 1,-1, 1, 7] ]
        b = [7, 1, 3, -2]

    elif preg == 2:
        A = [ [ 4,-1, 0, 2],
              [-1, 4,-1, 0],
              [ 0,-1, 4, 1],
              [ 2, 0, 1, 3] ]
        b = [ 6, 3, 16, 12]

    else: # if preg == 3:
        A = read_Matrix()
        print(NP.array(A))
        b = read_Vector(len(A))
        print(NP.array(b))

    print("\nMetodo:\t Cholesky")
    print("\nSistema:")
    printSistemaCholesky(preg)

    MetodoCholesky(A, b)

def printSistemaCholesky(preg): #Imprime los sistemas disponibles para Cholesky
    
  if preg == 1:
    print("\t+9X\t-3Y\t+2Z\t+1W\t=  7\n",
          "\t-3X\t+5Y\t+0Z\t-1W\t=  1\n",
          "\t+2X\t+0Y\t+4Z\t+1W\t=  3\n",
          "\t+1X\t-1Y\t+1Z\t+7W\t= -2\n")

  elif preg == 2:
    print("\t+4X\t-1Y\t+0Z\t+2W\t=  6\n",
          "\t-1X\t+4Y\t-1Z\t+0W\t=  3\n",
          "\t+0X\t-1Y\t+4Z\t+1W\t=  16\n",
          "\t+2X\t+0Y\t+1Z\t+3W\t=  12\n")




#Ejecucion del programa

bienvenida()

while True: #Se cicla hasta que el usuario elija salir

    met = selectMet()

    if met == 1:
        menuJacobi()

    elif met == 2:
        menuCholesky()

    else:  # if met == 3
        break
    '''
author: Christopher Xavier Sanchez Duran
// Github: https://github.com/cxavier23
// Portafolio: https://cxavier23.github.io/
'''