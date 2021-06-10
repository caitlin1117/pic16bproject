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
    """
    solve 1d transport equation u_t + u_x = 0 using Lax Wendroff method periodic BCs
    Input: (u_0, pl, save)
    u_0: initial condition u(x, 0) = u_0, 
    pl: can be specified to "animation", "2d", or "3d" for different plot,default is animation
    save: specified to True if you want the animation(mp4 file) or figure(png file) saved, default is false
    Output: animation or figure of the solution depending on different kind of plot
    """
    p = plot(u_0,pl)
    if save:
        if pl=="animation":
            p.save('animation.mp4')
        else:
            p.savefig('plot.png') 
            
    return p


def solve_wendroff(u_0):
    '''
    solve using Lax Wendroff method
    '''
    N = 100
    tmax = 1.5
    h = (1 - 0) / 100
    x = np.linspace(-h, 1+h, 103)
    u0 = u_0(x)
    k = 0.009
    lamda = k/h
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*103).reshape(nsteps+1, 103)
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
    '''
    plot solution of 1d transport equation
    '''
    u = solve_wendroff(u_0)
    init_state = u[0,:]
    if(pl=="animation"):
        fig, ax = plt.subplots()
        x = np.linspace(-0.01, 1+0.01, 103)
        line1, = ax.plot(x, init_state)
        ax.set_xlim((0, 1))
        ax.set_ylabel("$u$")
        ax.set_title('$u_t+u_x=0$ travelling wave simulation');
        def update(step): #update each time step   
            state = u[step,:]
            line1.set_data(x, state)   
            return
        ani = animation.FuncAnimation(fig, update, frames = 166, interval=60)
        plt.show()
        return ani
    elif pl=="2d":
        fig = plt.figure(figsize=(5,5))
        plt.imshow(u[165::-1, :], extent=[0,100,0,166])
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('t')
        plt.title('$u_t+u_x=0$ on a periodic domain')
        plt.show()
        return fig
    else:
        x = np.linspace(-0.01, 1+0.01, 103)
        y = np.linspace(0,70, 71)
        x, y = np.meshgrid(x, y)
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca(projection='3d')
        ax.plot_surface(x,y, u[0:71,:], rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
            
        ax.set_title('$u_t+u_x=0$');
        ax.set_xlabel('x')
        ax.set_ylabel('t')
        ax.set_zlabel('u');
        plt.show()
        return fig
    


