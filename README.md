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
    return 0.8*np.exp(-700*((x-1)**2+(y-0)**2))
def Gaussian_2d_2(x,y):
    return 0.8*np.exp(-700*((x-0.5)**2+(y-0.5)**2))
def Gaussian_2d_3(x,y):
    return 0.8*np.exp(-700*((x-1)**2+(y-0)**2))+0.8*np.exp(-700*((x-0)**2+(y-1)**2))
def f2(x,y):
    ic01 = np.logical_and(x >= 1/4, x <= 3/4)
    ic02 = np.logical_and(y >= 1/4, y <= 3/4)
    return np.multiply(ic01, ic02)

```
### Transport equation
Use the solution function    
Input: (u_0, pl, save)
u_0: initial condition u(x, 0) = u_0, 
pl: can be specified to "animation", "2d", or "3d" for different plot, default is animation
save: specified to True if you want the animation(mp4 file) or figure(png file) saved, default is false
Output: animation or figure of the solution depending on different kind of plot
```python
te.solution(Gaussian_2,"animation") 
```
![](/test/transport_equation.mp4)

<video width="640" height="480" controls>
  <source src="/test/transport_equation.mp4" type="video/mp4">
</video>







