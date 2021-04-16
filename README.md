# pic16bproject

### Abstract
In this project, I want to learn about solving PDEs (with a focus on hyperbolic and parabolic equations) using finite difference methods and analyze some particular PDEs with different boundary conditions and initial conditions numerically and visualize the solution. We need to implement the numerical methods in python and visualize the numerical solution.



### Planned Deliverables

Full success : A Python package that is designed to visualize the solution to some PDES with specified parameters, initial and boundary conditions and finite different method. The PDEs include 1-dimensional advection equation, heat equation, wave equation, and some nonlinear PDEs such as Burger's equation and KdV equation. The methods include explicit and implicit methods. The visualization includes 3d visualization and animation. 

Partial success: will not be able to implement all of the planned methods or will not be able to implement the method for the KdV equation.


### Resources Required
The main resources required are textbooks and notes on solving PDEs using finite difference method. Currently, materials I want to use as reference includes:  
Chapter 1, 2, 6 from Finite difference schemes and partial differential equations by John Strikwerda   
Numerical solution of partial differential equations by K. W. Morton, D. F. Mayers    
Lecture notes from http://people.bu.edu/andasari/courses/numericalpython/python.html  
Lecture notes from https://espace.library.uq.edu.au/data/UQ_239427/Lectures_Book.pdf?Expires=1618541682&Key-Pair-Id=APKAJKNBJ4MJBJNC6NLQ&Signature=MT~pUUExCCXmpcg78IGxNd6Mkfbo7k9mAP4-Wo2DLIgiLTYixjwlvvExQ2UJG53vCX54gyXFV0e8njvb4SiVpNj1M7zhZML7l2XVrhkffoT4OWqjb65shSZUJw0g0oLqEUUuTTVzlvznT4GiaL1~ZQP~PpPfblCdj4ylC~6TVjYYEsBIvtkwUBWjY7OMicZFSg-uOVWGsxcFbvVpPfhusIV7kl7VdabC2M03UrzOT29CrcCP0uM3boHUMBwQ~lkqypa7W41Gbytdy61XdRXozcBFH-RczPcQB3rknEUr7DwWtv2BULMQD8qvaq9SFXnTZR1to8bhKJJITV6MtnYmUg__  
KdV Notes https://newtraell.cs.uchicago.edu/files/ms_paper/hmmorgan.pdf  
Burgers equation http://www.bcamath.org/projects/NUMERIWAVES/Burgers_Equation_M_Landajuela.pdf


### Tools/Skills Required:
numpy  
complex visualization:  
matplotlib.animation as animation  
matplotlib colors  
matplotlib.pyplot   
plotly  
for solving matrix: scipy linalg


### Risks: 
What are two things that could potentially stop you from achieving the full deliverable above?   
Even though deriving the different schemes may not be that difficult, analyzing the schemes	to determine if they are useful approximations to the differential equation is complicated. I probably don’t have enough time to thoroughly look into the convergence and stability analysis. 

### Ethics:
Currently I cannot think of any potential biases or harms from the project. I think if we can try to solve a problem both analytically and numerically, we can gain better insights. When we are using analytical approaches, it would be nice if we can verify our result using numerical approaches.

### Tentative Timeline:
Week2  
1D advection equation   
Explicit methods: Lax-Friedrichs scheme and Leapfrog scheme (multistep),  Lax-Wendroff scheme

Week4  
Second order 1-d wave equation
Explicit difference method and Implicit difference method

Week6  
Heat equation  
Forward difference in time and central differences in space (FTCS method) explicit method
Backward-time Central-Space (BTCS method) which is implicit and unconditionally stable

Remaining week  
Inviscid Burger’s equation: Lax-Wendroff   
KdV equation: FTBS for the first two terms, then derive the approximation for $$u_{xxx}$$


