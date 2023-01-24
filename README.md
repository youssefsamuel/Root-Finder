# Root finder
Implementing a program to find the root of an equation using one of the following methods:
1. Bisection
2. False position
3. Fixed point
4. Newton-Raphson
5. Secant

The root finder program takes as an input the equation, the technique to use and its required parameters (e.g., interval for the bisection method).



## Specifications
The program contains the following features:
1.  An interactive GUI that enables the user to enter equations containing different 
functions such as: {poly, exp, cos, sin}. Reading from files is available as 
well.
2.  Differentiation and Parsing the input.
3.  A way to choose a method to solve the given equation.
4.  A way to enter the precision and the max number of iterations otherwise default 
values are used, Default Max Iterations = 50, Default Epsilon = 0.00001.
5.  The answer for the chosen method indicating the number of iterations, 
execution time, all iterations, approximate root, and precision.
6.  Compute the theoretical bound of the error for the methods.
