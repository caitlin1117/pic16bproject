"""
Finite-difference solver for transport equation:
    u_t + u_x = 0.
Initial conditions:
    u(x, 0) = u_0.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['animation.ffmpeg_path'] = '/Users/caitlin/opt/anaconda3/bin/ffmpeg'
import matplotlib.animation as animation
from matplotlib import colors
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from IPython.display import HTML
plt.style.use('seaborn-pastel')

def solution(u_0, pl="animation",save=False):
    p = plot(u_0,pl)
    if save:
        if p=="animation":
            p.save('animation.mp4')
        else:
            p.savefig('plot.png') 
            
    return p

def solve(u_0):
    N = 100
    tmax = 1.5
    h = (1 - 0) / 100
    x = np.linspace(-h, 1+h, 103)
    u0 = u_0(x)
    k = 0.009
    lamda = k/h
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*103).reshape(nsteps+1, 103)
    #print(np.shape(u))
    u[0] = u0
    temp = u0.copy()
    for i in range(1,int(nsteps)):
        for j in range(N+2):
            temp[j] = 0.5*(u[i-1,j+1]+u[i-1,j-1])-(lamda/2)*(u[i-1,j+1]-u[i-1,j-1])
        u[i] = temp
        u[i, 0] = u[i, N+1]
        u[i, N+2] = u[i,1]

    return u
def solve_wendroff(u_0):
    N = 100
    tmax = 1.5
    h = (1 - 0) / 100
    x = np.linspace(-h, 1+h, 103)
    u0 = u_0(x)
    k = 0.009
    lamda = k/h
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*103).reshape(nsteps+1, 103)
    #print(np.shape(u))
    u[0] = u0
    temp = u0.copy()
    for i in range(1,int(nsteps)):
        for j in range(N+2):
            temp[j] = u[i-1,j]-0.5*lamda*(u[i-1,j+1]-u[i-1,j-1])+0.5*(lamda**2)*(u[i-1,j+1]-2*u[i-1,j]+u[i-1,j-1])
        u[i] = temp
        u[i, 0] = u[i, N+1]
        u[i, N+2] = u[i,1]

    return u
    

def plot(u_0,pl="animation"):
    #u = solve(u_0)
    u = solve_wendroff(u_0)
    #e = exact(u_0)
    init_state = u[0,:]
    if(pl=="animation"):
        fig, ax = plt.subplots()
        x = np.linspace(-0.01, 1+0.01, 103)
        line1, = ax.plot(x, init_state)
    #line2, = ax.plot(x, init_state) 
        ax.set_xlim((0, 1))
    #ax.set_ylim((0, 0.6)) 
        ax.set_ylabel("$u$")
        ax.set_title('$u_t+u_x=0$ travelling wave simulation');
        def update(step):    
            state = u[step,:]
        #state2 = e[step,:]#u_0(x-c*step*0.009)
            line1.set_data(x, state)
        #line2.set_data(x, state2)
            return
        ani = animation.FuncAnimation(fig, update, frames = 166, interval=60)
        #plt.show()
        return ani
    elif pl=="2d":
        fig = plt.figure(figsize=(5,5))
        plt.imshow(u[165::-1, :], extent=[0,100,0,166])
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('t')
        plt.title('$u_t+u_x=0$ on a periodic domain')
        #plt.show()
        return fig
    else:
        x = np.linspace(-0.01, 1+0.01, 103)
        y = np.linspace(0,70, 71)
        x, y = np.meshgrid(x, y)
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca(projection='3d')
        print(u[0:71,:].size)
        ax.plot_surface(x,y, u[0:71,:], rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
            
        ax.set_title('$u_t+u_x=0$');
        ax.set_xlabel('x')
        ax.set_ylabel('t')
        ax.set_zlabel('u');
        #plt.show()
        return fig
    

def f(x):
    return np.exp(-200*(x-0.25)**2)

#solution(f,"3d")
 

