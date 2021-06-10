This Python package is designed to help visualize the solution to some PDEs with specified initial conditions. The PDEs include 1d advection equation, 1d and 2d heat equation, 1d and 2d wave equation, and Burger's equation. The numerical solutions are solved using finite difference methods. My intended audience are students who have a little bit knowledge of PDEs. In an introductory PDE course, it can be confusing when you learn about the analytic solution. In order to get a better understanding of how the solution looks like, they can use this package to visualize the solution to some basic PDEs.   

### Installation:
```python
pip install PdeSimulation==0.0.7
```

### Example Usage
```python
from PdeSimulation import transport_equation as te
from PdeSimulation import wave_equation as we
from PdeSimulation import heat_equation as he
from PdeSimulation import burgers as be
import numpy as np

#initial profiles
def sine(x):
    return np.sin(2*x*np.pi)
def Gaussian(x):
    return np.exp(-200*(x-0.5)**2)
def Gaussian_2(x):
    return np.exp(-200*(x-0.25)**2)
def Gaussian_2d_1(x,y):
    return -0.8*np.exp(-200*((x-1)**2+(y-0)**2))
def Gaussian_2d_2(x,y):
    return 0.8*np.exp(-200*((x-0.5)**2+(y-0.5)**2))
def Gaussian_2d_3(x,y):
    return 0.8*np.exp(-700*((x-1)**2+(y-0)**2))+0.8*np.exp(-700*((x-0)**2+(y-1)**2))
def f2(x,y):
    ic01 = np.logical_and(x >= 1/4, x <= 3/4)
    ic02 = np.logical_and(y >= 1/4, y <= 3/4)
    return np.multiply(ic01, ic02)

```
### Transport equation
Use the `solution` function    
solve 1d transport equation u_t + u_x = 0 using Lax Wendroff method periodic BCs
Input: (u_0, pl, save)     
u_0: initial condition u(x, 0) = u_0,      
pl: can be specified to "animation", "2d", or "3d" for different plot, default is animation     
save: specified to True if you want the animation(mp4 file) or figure(png file) saved, default is false.   
Output: animation or figure of the solution depending on different kind of plot    
```python
te.solution(Gaussian_2,"animation") 
```
<img src="/test/transport_equation.gif" width="320" height="240"/>    

```python
te.solution(Gaussian_2,"2d") 
```
<img src="/test/transport_equation_2d.png" width="320" height="240"/> 

```python
te.solution(Gaussian_2,"3d") 
```
<img src="/test/transport_equation_3d.png" width="320" height="240"/> 

### Wave equation
Use the `solution` function   
The 1d wave equation u_tt = u_xx in [0,1] is solved using CTCS method, with zero initial velocity and homogeneous Dirichlet BCs    
The 2d wave equation u_tt = u_xx+u_yy in [0,1]x[0,1] is solved using CTCS method, with zero initial velocity and homogeneous Neumann BCs    
Input: (u_0, dim, pl, save)    
u_0: u_0 is the initial profile   
dim: can be specified to 1 or 2 indicating you want to solve a 1d or 2d wave equation    
pl: for 1d wave equation, pl can be specified to animation, 2d or 3d for different types of plot, default is animation    
    for 2d wave equation, pl can only be specified to animation, default is animation    
save: specified to True if you want the animation(mp4 file) or figure(png file) saved, default is false    
Output:animation or figure of the solution depending on different kind of plot    

```python
we.solution(sine,1,"animation") 
```
<img src="/test/wave_1d.gif" width="320" height="240"/> 

```python
we.solution(sine,1,"2d") 
```
<img src="/test/wave_1d_2dplot.png" width="320" height="240"/> 

```python
we.solution(sine,1,"3d") 
```
<img src="/test/wave_1d_3dplot.png" width="320" height="240"/> 

```python
we.solution(Gaussian_2d_1,2,"animation") 
```
<img src="/test/wave_2d_2.gif" width="320" height="240"/> 

```python
we.solution(Gaussian_2d_2,2,"animation") 
```
<img src="/test/wave_2d.gif" width="320" height="240"/> 


### Heat Equation
Use the `solution` function    
The 1d heat equation u_t = u_xx in [0,1] with homogeneous Dirichlet BCs can be solved   
using Forward Time Central Space (explicit) or Backward Time Central Space (implicit) methods    
The 2d heat equation u_t = u_xx+u_yy in [0,1]x[0,1] with homogeneous Neumann BCs is solved using the ADI method   
Input: (u_0, dim, pl, method, save)    
u_0: u_0 is the initial profile    
dim: can be specified to 1 or 2 indicating you want to solve a 1d or 2d heat equation    
pl: for 1d heat equation, pl can be specified to animation, 2d or 3d for different types of plot, default is animation    
    for 2d heat equation, pl can only be specified to animation, default is animation    
method: for 1d heat equation, the method can be specified to implicit or explicit, default is implicit    
save: specified to True if you want the animation(mp4 file) or figure(png file) saved, default is false    
Output:animation or figure of the solution depending on different kind of plot    

```python
he.solution(Gaussian,1,method = "explicit") 
```
<img src="/test/heat_1d_FTCS.gif" width="320" height="240"/> 

```python
he.solution(Gaussian,1,method = "implicit") 
```
<img src="/test/heat_1d_BTCS.gif" width="320" height="240"/> 

```python
he.solution(Gaussian,1, pl = "2d", method = "implicit") 
```
<img src="/test/heat_1d_2dplot.png" width="320" height="240"/> 

```python
he.solution(Gaussian,1, pl = "3d", method = "implicit") 
```
<img src="/test/heat_1d_3dplot.png" width="320" height="240"/> 


```python
he.solution(f2,2) 
```
<img src="/test/heat_2d.gif" width="320" height="240"/> 

### Burgers Equation
Use the `solution` function   
solve the Riemann problem for the inviscid Burgers equation u_t + uu_x = 0 using upwind conservative method    
Input: (type, save)   
type: can be specified to shock or rarefaction   
save: specified to True if you want the animation(mp4 file) saved, default is false   
Output: animation of the solution depending on different kind of plot   


```python
be.solution("rarefaction")
```
<img src="/test/burgers_equation_rarefaction.gif" width="320" height="240"/> 

```python
be.solution("shock")
```
<img src="/test/burgers_equation_shock.gif" width="320" height="240"/> 

