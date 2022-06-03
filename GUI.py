"""
This file contains numerical methods for solving:
    Systems of linear equations
    Non linear equations roots
    Interpolation problems
Author: Juan Pablo Madrid Peláez - 2022
"""
import PySimpleGUI as psg
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpmath import *
#set the theme for the screen/window
psg.theme('SandyBeach')

def incrementalsearch(funcion,parameters,method = "incremental"):
    x = sp.symbols("x")
    f1 = sp.sympify(funcion)
    f = sp.lambdify(x,f1)
    x = parameters[0]
    x1 = parameters[1]
    delta = parameters[2]
    Found = False
    while x<x1 and not Found:
        mult = f(x)*f(x+delta)
        if mult <= 0:
            Found = True
            if method == "incremental":
                psg.popup("The function has a root in the interval:", "[" + str(x) + "," + str(x+delta) + "]")
        x = x + delta
    if not Found:
        if method == "incremental":
            psg.popup("No roots in the interval")
    return Found
    
#SOLVERS
def biseccion(funcion,parameters):
    # Ingreso datos de entrada para los diferentes métodos a trabajar
    x = sp.symbols("x")
    f1 = sp.sympify(funcion)
    f = sp.lambdify(x,f1)
    a = parameters[0]
    b = parameters[1]
    # guarda valores iniciales
    a0 = a
    b0 = b
    # guarda valores iniciales del error y del número de iteraciones
    tol = parameters[2]  # float(input("Ingrese el valor de la tolerancia: "))
    nmax = 100  # float(input("Ingrese el número máximo de iteraciones: "))
    error = 100
    niter = 0
    has_sol = incrementalsearch(funcion,[a,b,0.0001],"other")
    if not has_sol:
        psg.popup("The interval has no roots")
    else:
    # Método de Bisección
        # evaluo primer valor medio
        m = a + (b - a) / 2
        # Evaluacion de la función en los puntos a, b y m
        fa = f(a)
        fb = f(b)
        fm = f(m)
        results = "# iter\t\t a \t\t f(a) \t\t b \t\t f(b) \t\t m \t\t f(m) \t\t error"
        results = results + '\n' + "{0} \t\t {1:6.4f} \t {2:6.4f} \t {3:6.4f} \t {4:6.4f} \t {5:6.4f} \t {6:6.4f} \t {7:6.4f}".format(niter, a0, fa, b0, fb, m, fm, error)
        # ciclo iterativo
        while error > tol and niter < nmax:
            m = a + (b - a) / 2
            if np.sign(fa) == np.sign(fm):
                a = m
                fa = f(a)
            else:
                b = m
                fb = f(b)
            m = a + (b - a) / 2
            fm = f(m)
            error = abs(b - a)
            niter += 1
            results = results + '\n' + "{0} \t\t {1:6.6f} \t {2:6.6f} \t {3:6.6f} \t {4:6.6f} \t {5:6.6f} \t {6:6.6f} \t {7:6.6f}".format(niter, a, fa, b, fb, m, fm, error)
        psg.popup("La raíz de la función dada en el intervalo [{0:6.4f},{1:6.4f}] es {2:6.7f}".format(a0, b0, m))
        psg.popup("PASO A PASO",results, line_width=300)
    return "True"

def regulaFalsi(funcion,parameters):
    x = sp.symbols("x")
    f1 = sp.sympify(funcion)
    f = sp.lambdify(x,f1)
    a = parameters[0]
    b = parameters[1]
    # guarda valores iniciales
    a0 = a
    b0 = b
    # guarda valores iniciales del error y del número de iteraciones
    tol = parameters[2]
    nmax = 100
    error = 100
    niter = 0
    # Evaluacion de la función en los puntos a, b y m
    fa = f(a)
    fb = f(b)
    # evaluo primer valor medio
    m = a - fa * (b - a) / (fb - fa)

    fm = f(m)
    has_sol = incrementalsearch(funcion,[a,b,0.0001],"other")
    if not has_sol:
        psg.popup("The interval has no roots")
    else:

        print("# iter\t\t a \t\t f(a) \t\t b \t\t f(b) \t\t m \t\t f(m) \t\t error")
        print(
            "{0} \t\t {1:6.4f} \t {2:6.4f} \t {3:6.4f} \t {4:6.4f} \t {5:6.4f} \t {6:6.4f} \t {7:6.4f}".format(
                niter, a0, fa, b0, fb, m, fm, error
            )
        )

        # ciclo iterativo
        while error > tol and niter < nmax:
            m = a - fa * (b - a) / (fb - fa)
            if np.sign(fa) == np.sign(fm):
                a = m
                fa = f(a)
            else:
                b = m
                fb = f(b)

            m = a - fa * (b - a) / (fb - fa)
            fm = f(m)
            error = abs(b - a)
            niter += 1
            print(
                "{0} \t\t {1:6.6f} \t {2:6.6f} \t {3:6.6f} \t {4:6.6f} \t {5:6.6f} \t {6:6.6f} \t {7:6.6f}".format(
                    niter, a, fa, b, fb, m, fm, error
                )
            )

        psg.popup("La raíz de la función dada en el intervalo [{0:6.4f},{1:6.4f}] es {2:6.7f}".format(a0, b0, m))

        print("Error : %s" % ((b0 - a0) / 2 ** niter))

def fixed_point(funcion, parameters):
    #Open second window to insert G
    layout2=[[psg.Text('Insert function G',size=(30, 1), font='Lucida',justification='left')],
            [psg.Input(key='function G')],
            [psg.Button('RUN', font=('Times New Roman',12)),psg.Button('CANCEL', font=('Times New Roman',12))]]
    win2 =psg.Window('Fixed Point',layout2)
    e2,v2=win2.read()
    win2.close()
    equationG = v2['function G']

    x = sp.symbols("x")
    f1 = sp.sympify(funcion)
    f = sp.lambdify(x,f1)
    g1 = sp.sympify(equationG)
    g = sp.lambdify(x,g1)
    p = parameters[0]
    TOL = parameters[1]
    error = 1
    iterations = 0
    while error > TOL:
        p_new = g(p)
        error = abs(p_new - p)
        p = p_new
        iterations += 1
        print(f'p{iterations} = {p: 0.5f}')
    psg.popup(f'Root: {p}\nIterations: {iterations}')

def newtonRaphson(funcion, parameters):
    x = sp.symbols("x")
    f = sp.sympify(funcion)
    df = sp.diff(f)
    f = sp.lambdify(x,f)
    df = sp.lambdify(x,df)
    x0 = parameters[0]
    TOL = parameters[1]
    error = 1
    iterations = 0
    while error > TOL:
        new_x = x0 - f(x0)/df(x0)
        error = abs(new_x - x0)
        x0 = new_x
        iterations += 1
        print(f'x{iterations}: {x0}')
    psg.popup(f"Solution = {x0:.15f},Iterations: {iterations}")

def secant(funcion,parameters):
    x = sp.symbols("x")
    f1 = sp.sympify(funcion)
    f = sp.lambdify(x,f1)
    x1 = parameters[0]
    x2 = parameters[1]
    TOL = parameters[2]
    has_sol = incrementalsearch(funcion,[a,b,0.0001],"other")
    if not has_sol:
        psg.popup("The interval has no roots")
    else:
        error = 1
        iterations = 0
        print(f'x0: {x1}\nx1: {x2}')
        while error > TOL:
            x_new = x2 - (f(x2)*(x2-x1))/(f(x2)-f(x1))
            error = abs(x_new - x2)
            x1 = x2
            x2 = x_new
            iterations += 1
            print(f'x{iterations+1}: {x2}')
        psg.popup(f'Result: {x2}\nIterations: {iterations}')

def eliminacionInf(Ab):
    n = Ab.shape[0]
    for k in range(0, n - 1):
        # print("\n\n Etapa", k + 1, "Eliminación columna:", k + 1, "\n")
        for i in range(k + 1, n):
            multiplicador = Ab[i][k] / Ab[k][k]
            for j in range(k, n + 1):
                Ab[i][j] = Ab[i][j] - multiplicador * Ab[k][j]
        # print("\n", Ab)
    print(Ab)
    return Ab

def eliminacionSup(Ab):
    n = Ab.shape[0]
    for k in range(n - 1, 0, -1):
        for i in range(0, k):
            multiplicador = Ab[i][k] / Ab[k][k]
            for j in range(1, n + 1):
                Ab[i][j] = Ab[i][j] - multiplicador * Ab[k][j]
    print(Ab)
    return Ab

def sustitucionReg(Ab):
    n = Ab.shape[0]
    x = np.zeros(n)
    x[n - 1] = Ab[n - 1][n] / Ab[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        suma = 0
        for j in range(i + 1, n):
            suma = suma + Ab[i][j] * x[j]
        x[i] = (Ab[i][n] - suma) / Ab[i][i]
    [print("x[{0}] = {1:6.4f}".format(i, x[i])) for i in range(n)]
    return x

def sustitucionProg(Ab):
    n = Ab.shape[0]
    x = np.zeros(n)
    x[0] = Ab[0][n]
    for i in range(1, n):
        suma = 0
        for j in range(i - 1, -1, -1):
            suma = suma + Ab[i][j] * x[j]
        x[i] = Ab[i][n] - suma
    [print("x[{0}] = {1:6.4f}".format(i, x[i])) for i in range(n)]
    return x

def pivoteo(Ab):
    n = Ab.shape[0]
    for i in range(n - 1):
        columna = abs(Ab[i:, i])  # i en adelante
        dondemax = np.argmax(columna)
        if dondemax != 0:
            temporal = np.copy(Ab[i, :])
            Ab[i, :] = Ab[dondemax + i, :]
            Ab[dondemax + i, :] = temporal
    print(Ab)

def escalamiento(Ab):
    n = Ab.shape[0]
    A = np.delete(Ab, n, axis=1)
    for i in range(n):
        filamax = np.max(abs(A[i, :]))
        print(filamax)
        Ab[i, :] = Ab[i, :] / filamax
    print(Ab)

def gauss_Simple(Ab):
    Ab = eliminacionInf(Ab)
    x = sustitucionReg(Ab)
    psg.popup("the Solution vector is:",x)
    return x

def gauss_Jordan(Ab):
    Ab = eliminacionInf(Ab)
    Ab = eliminacionSup(Ab)
    n = Ab.shape[0]
    x = np.zeros(n)
    for i in range(n):
        x[i] = Ab[i][n] / Ab[i][i]
    psg.popup("the Solution vector is:", x)

def erroresReondeo(A):
    x = gauss_Simple(Ab)
    bGS = A.dot(x)
    [print("bGS[{0}] = {1:6.20f}".format(i, bGS[i])) for i in range(len(x))]
    # Error relativo %
    [
        print("error[{0}] = {1:6.20f}".format(i, abs(bGS[i] - b[i]) / abs(b[i]) * 100))
        for i in range(len(x))
    ]

def factorizacion(A, b):
    n = A.shape[0]
    U = A
    L = np.identity(n)
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            multiplicador = U[i][k] / U[k][k]
            L[i][k] = multiplicador
            for j in range(k, n):
                U[i][j] = U[i][j] - multiplicador * U[k][j]
    return (L, U)

def factorizacionLU(A, b):
    (L, U) = factorizacion(A, b)
    n = L.shape[0]
    Lb = np.c_[L, b]
    z = np.zeros(n)
    z[0] = Lb[0][n]
    for i in range(1, n):
        suma = 0
        for j in range(i - 1, -1, -1):
            suma = suma + Lb[i][j] * z[j]
        z[i] = Lb[i][n] - suma
    Uz = np.c_[U, z]
    psg.popup("the Solution vector is:",sustitucionReg(Uz))

def nuevoJacobi(A, b, x):
    x1 = np.copy(x)
    x0 = x
    n = x.shape[0]
    for i in range(n):
        suma = 0
        for j in range(n):
            if j != i:
                suma = suma + A[i, j] * x0[j]
        x[i] = b[i] - suma / A[i, i]
    return x

def jacobi(A, b, x0, tol, niter):
    n = A.shape[0]
    D = np.diag(np.diagonal(A))
    L = np.tril(A, k=-1) * -1
    U = np.triu(A, k=1) * -1
    Dinv = np.linalg.inv(D)
    Tj = np.dot(Dinv, L + U)
    Cj = np.dot(Dinv, b)
    dispersion = tol + 1
    contador = 0

    print("Iteracion\t  Resultado\t\t\t\tError")

    while (dispersion > tol) and (contador < niter):
        x1 = np.dot(Tj, x0) + Cj
        dispersion = np.linalg.norm(x1 - x0, ord=2)
        x0 = x1
        contador += 1
        print(
            "{0} \t\t {1:s} \t {2:6.6f}".format(
                contador, np.array2string(x0), dispersion
            )
        )
    if dispersion < tol:
        psg.popup(x1, "es aproximación con una tolerancia de {tol}".format(tol=tol))
    else:
        psg.popup("Fracaso en {niter} iteraciones".format(niter=niter))

def gaussSeidel(A, b, x0, tol, niter):
    n = A.shape[0]
    D = np.diag(np.diagonal(A))
    L = np.tril(A, k=-1) * -1
    U = np.triu(A, k=1) * -1
    Tg = np.dot(np.linalg.inv(D - L), U)
    print(Tg)
    Cg = np.dot(np.linalg.inv(D - L), b)
    dispersion = tol + 1
    contador = 0

    while (dispersion > tol) and (contador < niter):
        x1 = np.dot(Tg, x0) + Cg
        dispersion = np.linalg.norm(x1 - x0, ord=2)
        x0 = x1
        contador += 1
        print("\n\n")
        print(x0)

    if dispersion < tol:
        psg.popup(x1, "es aproximación con una tolerancia de {tol}".format(tol=tol))
    else:
        psg.popup("Fracaso en {niter} iteraciones".format(niter=niter))

def polynomial(p,x0):
    fx0 = p[0]
    x = sp.symbols("x")
    f = p[0]
    i = 1
    for c in p[1:]:
        f = f + c*x**i
        fx0 = fx0 + c*x0**i
        i += 1
    return (f,fx0)

def poly_newton_coefficient(x, y):
    """
    x: list or np array contanining x data points
    y: list or np array contanining y data points
    """

    m = len(x)

    x = np.copy(x)
    a = np.copy(y)
    for k in range(1, m):
        a[k:m] = (a[k:m] - a[k - 1])/(x[k:m] - x[k - 1])
    return a

def newton_polynomial(x_data, y_data, x):
    """
    x_data: data points at x
    y_data: data points at y
    x: evaluation point(s)
    """
    a = poly_newton_coefficient(x_data, y_data)
    n = len(x_data) - 1  # Degree of polynomial
    p = a[n]

    for k in range(1, n + 1):
        p = a[n - k] + (x - x_data[n - k])*p
    psg.popup("The coefs are: ", a, "the value at " + str(x)+ " is " +str(p))
    return p

def vandermonde(xs,ys,x0):
    n = len(ys)
    x = np.array([[i] for i in xs])
    y = np.array([[j] for j in ys])
    A = np.ones((n,n))
    Acol = np.ones((n,1))
    for i in range(1,n):
        Acol = np.multiply(Acol,x)
        A[:,i] = Acol[:,0]
    a = np.linalg.solve(A,y)
    coefs = [i[0] for i in a]
    pol = polynomial(coefs,x0)[0]
    fx0 = polynomial(coefs,x0)[1]
    psg.popup("The coefs are: ", coefs,"The polynomial is:", pol, "Evaluated at " +str(x0)+ " is equal to " + str(fx0))
    return a

def lagrange(x,y,xi):
    # Initialize result
    n = len(x)
    f = [(x[j],y[j]) for j in range(n)]
    result = 0.0
    for i in range(n):
        # Compute individual terms of formula
        term = f[i][1]
        for j in range(n):
            if j != i:
                term = term * (xi - f[j][0]) / (f[i][0] - f[j][0])
        # Add current term to result
        result += term
    psg.popup("The value at "+ str(xi) +" is:", result)
    return result

def interpolate():
    layout1=[[psg.Text('Select the method',size=(20, 1), font='Lucida',justification='left')],
        [psg.Combo(['Divided difference','Vandermonde','Lagrange'],default_value='Vandermonde',key='method')],
        [psg.Text('Insert X values separated by comma',size=(30, 1), font='Lucida',justification='left')],
        [psg.Input(key='Xs')],
        [psg.Text('Insert Y values separated by comma',size=(30, 1), font='Lucida',justification='left')],
        [psg.Input(key='Ys')],
        [psg.Text('Point to evaluate (optional)',size=(30, 1), font='Lucida',justification='left')],
        [psg.Input(key='eval')],
        [psg.Button('RUN', font=('Times New Roman',12)),psg.Button('CANCEL', font=('Times New Roman',12))]]
    #Define Window
    win1 =psg.Window('Interpolation',layout1)
    #Read  values entered by user
    e1,v1=win1.read()
    #close first window
    win1.close()
    #preprocess inputs
    X = list(map(float,v1['Xs'].split(",")))
    Y = list(map(float,v1['Ys'].split(",")))
    x0 = float(v1['eval'])
    #Call selected method
    method = v1['method']
    if method=='Divided difference':
        newton_polynomial(X,Y,x0)
    elif method == 'Vandermonde':
        vandermonde(X,Y,x0)
    elif method == 'Lagrange':
        lagrange(X,Y,x0)

def solve_system_of_equations():
    layout1=[[psg.Text('Select the method',size=(20, 1), font='Lucida',justification='left')],
        [psg.Combo(['Simple Gauss method','Gauss Jordan','LU Factorization','Jacobi','Gauss Seidel'],default_value='Gauss Jordan',key='method')],
        [psg.Text('Insert A: rows separated by semicolon and vals separated by comma',size=(30, 1), font='Lucida',justification='left')],
        [psg.Input(key='A')],
        [psg.Text('Insert b: vals separated by comma',size=(30, 1), font='Lucida',justification='left')],
        [psg.Input(key='b')],
        [psg.Text('Insert tol and niter separated by comma',size=(30, 1), font='Lucida',justification='left')],
        [psg.Input(key='pars', default_text="0.001,100")],
        [psg.Text('Insert initial vals if required separated by comma',size=(30, 1), font='Lucida',justification='left')],
        [psg.Input(key='initial', default_text="0,0,0")],
        [psg.Button('RUN', font=('Times New Roman',12)),psg.Button('CANCEL', font=('Times New Roman',12))]]
    #Define Window
    win1 =psg.Window('Sistem of equations',layout1)
    #Read  values entered by user
    e1,v1=win1.read()
    #close first window
    win1.close()
    #preprocess inputs
    rowsA = v1['A'].split(";")
    entriesA = [list(map(float,r.split(","))) for r in rowsA]
    A = np.array(entriesA)
    b_list = list(map(float,v1['b'].split(",")))
    b = np.array(b_list)
    x0_list = list(map(float,v1['initial'].split(",")))
    x0 = np.array(x0_list)
    parameters = list(map(float,v1['pars'].split(",")))
    tol = parameters[0]
    niter = parameters[1]
    Ab = np.c_[A, b]
    #Call selected method
    method = v1['method']
    if method=='Simple Gauss method':
        gauss_Simple(Ab)
    elif method == 'Gauss Jordan':
        gauss_Jordan(Ab)
    elif method == 'LU Factorization':
        factorizacionLU(A,b)
    elif method == 'Jacobi':
        jacobi(A,b,x0,tol,niter)
    elif method == 'Gauss Seidel':
        gaussSeidel(A,b,x0,tol,niter)

def solve_non_linear_equation():
    layout1=[[psg.Text('Select the method',size=(20, 1), font='Lucida',justification='left')],
            [psg.Combo(['Incremental Search','Bisection','Regula falsi','Fixed Point','Newton Raphson','Secant'],default_value='Bisection',key='method')],
            [psg.Text('Enter parameters',size=(30, 1), font='Lucida',justification='left')],
            [psg.Input(key='parameters')],
            [psg.Text('Insert function',size=(30, 1), font='Lucida',justification='left')],
            [psg.Input(key='function')],
            [psg.Button('RUN', font=('Times New Roman',12)),psg.Button('CANCEL', font=('Times New Roman',12))]]
    #Define Window
    win1 =psg.Window('Solver',layout1)
    #Read  values entered by user
    e1,v1=win1.read()
    #close first window
    win1.close()
    #process inputs
    method = v1['method']
    pars = v1['parameters'].split(",")
    parameters = list(map(float,pars))
    equation = v1['function']

    if method == 'Incremental Search':
        incrementalsearch(equation,parameters)
    if method=="Bisection":
        biseccion(equation,parameters)
    elif method == "Regula falsi":
        regulaFalsi(equation,parameters)
    elif method == "Fixed Point":
        fixed_point(equation,parameters)
    elif method == "Newton Raphson":
        newtonRaphson(equation,parameters)
    elif method == "Secant":
        secant(equation,parameters)
    #display string in a popup
    psg.popup('SUCCESS','Method: ', method, 'Parameters: ', parameters, 'Equation: ', equation, 'OK for run again')

while True:
    layout=[[psg.Text('What will you solve?',size=(20, 1), font='Lucida',justification='left')],
            [psg.Combo(['Find roots of nonlinear equation','System of equations','Interpolate points'],default_value='Find roots of nonlinear equation',key='type')],
            [psg.Button('RUN', font=('Times New Roman',12)),psg.Button('CANCEL', font=('Times New Roman',12))]]
    #Define Window
    win =psg.Window('Solver',layout)
    #Read  values entered by user
    e1,v=win.read()
    #close first window
    win.close()
    #process inputs
    problem = v['type']

    if problem=='Find roots of nonlinear equation':
        solve_non_linear_equation()
    elif problem == 'System of equations':
        solve_system_of_equations()
    elif problem == 'Interpolate points':
        interpolate()
    psg.popup("OK for run again")