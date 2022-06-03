@author Juan Pablo Madrid Pelaez
The program contains all the methods seen in class
Below is the list of all contained methods and the order in which their parameters should be inserted.


METHODS FOR NON LINEAR EQUATIONS (ROOTS)
>incremental search
	params: x0,x1,delta
>Bisection
	params: a,b,tol
>RegulaFalsi
	params: a,b,tol
>Newton-Raphson
	params: x0,tol
>Fixed Point
	params: p0,tol
>secant
	params: x0,x1,tol


METHODS FOR SYSTEM OF EQUATIONS
(rows separated by ; and numbers separated by ,)
>Gauss
	params:A,b 
>LU 	
	params:A,b
>Jacobi
	params:A,b,x0,tol,niter
>Gauss-Seidel
	params:A,b,x0,tol,niter

METHODS FOR INTERPOLATION:
>Divided difference
	params:X values, Y values, x0
>Vandermonde
	params:X values, Y values, x0
>Lagrange
	params:X values, Y values, x0


Please define the functions with a clear syntax to be processed easily.
Some errors will be enhanced in a upcoming version