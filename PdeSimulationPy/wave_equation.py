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
from mpl_toolkits import mplot3d

def solution(u_0,dim, pl="animation",save=False):
    '''
    The 1d wave equation u_tt = u_xx in [0,1] is solved using CTCS method, with zero initial velocity and homogeneous Dirichlet BCs
    The 2d wave equation u_tt = u_xx+u_yy in [0,1]x[0,1] is solved using CTCS method, with zero initial velocity and homogeneous Neumann BCs
    Input: (u_0, dim, pl, save) 
    u_0: u_0 is the initial profile, 
    dim: can be specified to 1 or 2 indicating you want to solve a 1d or 2d wave equation
    pl: for 1d wave equation, pl can be specified to animation, 2d or 3d for different types of plot, default is animation
        for 2d wave equation, pl can only be specified to animation, default is animation
    save: specified to True if you want the animation(mp4 file) or figure(png file) saved, default is false
    Output:animation or figure of the solution depending on different kind of plot

    '''
    if dim==1:
        p = plot(u_0,pl)
        if save:
            if pl=="animation":
                p.save('animation.mp4')
            else:
                p.savefig("plot.png") 
        return p
    else:
        p = plot_2d(u_0)
        if save:
            p.save('animation.mp4')
        return p

    


def laplacian(u,n,i,j,im1,ip1,jm1,jp1):
    '''
    calculate Laplacian
    '''
    return u[n,ip1,j]+u[n,im1,j]+u[n,i,jp1]+u[n,i,jm1]-4*u[n,i,j]

def solve_2d(u_0):
    '''
    solve 2d wave equation
    '''
    N = 40
    tmax = 4
    h = (1 - 0) / N
    x = np.linspace(0, 1, N+1)
    y = np.linspace(0, 1, N+1)
    x, y = np.meshgrid(x, y)# sapce grid
    u0 = u_0(x,y)
    k = 0.015
    lamda = k/h
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*(N+1)*(N+1)).reshape(nsteps+1, (N+1), (N+1))# u is a 3 dim array, first dim is timestep
    u[0] = u0
    temp = u0.copy()
    # First take care of u[1] which is the solution at timestep 1
    # update all inner points
    for i in range(1,N):
            for j in range (1,N):
                temp[i,j] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i-1,i+1,j-1,j+1))/2
    # update boundary points
    i = 0
    for j in range(1,N):
        temp[0,j] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i+1,i+1,j-1,j+1))/2
    i = N
    for j in range(1,N):
        temp[N,j] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i-1,i-1,j-1,j+1))/2
    j = 0
    for i in range(1, N):
        temp[i,0] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i-1,i+1,j+1,j+1))/2
    j = N
    for i in range(1, N):
        temp[i,N] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i-1,i+1,j-1,j-1))/2
    #update all corner points
    i = 0
    j = 0
    temp[0,0] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i+1,i+1,j+1,j+1))/2
    i = N
    j = 0
    temp[N,0] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i-1,i-1,j+1,j+1))/2
    i = 0
    j = N
    temp[0,N] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i+1,i+1,j-1,j-1))/2
    i = N
    j = N
    temp[N,N] = (2*u[0,i,j]+lamda**2*laplacian(u,0,i,j,i-1,i-1,j-1,j-1))/2
    u[1] = temp

    #For the rest of the timestep
    for n in range(2,int(nsteps)):
        # update all inner points
        for i in range(1,N):
            for j in range (1,N):
                temp[i,j] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i-1,i+1,j-1,j+1)
        # update boundary points
        i = 0
        for j in range(1,N):
            temp[0,j] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i+1,i+1,j-1,j+1)
        i = N
        for j in range(1,N):
            temp[N,j] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i-1,i-1,j-1,j+1)
        j = 0
        for i in range(1, N):
            temp[i,0] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i-1,i+1,j+1,j+1)
        j = N
        for i in range(1, N):
            temp[i,N] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i-1,i+1,j-1,j-1)
        #update all corner points
        i = 0
        j = 0
        temp[0,0] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i+1,i+1,j+1,j+1)
        i = N
        j = 0
        temp[N,0] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i-1,i-1,j+1,j+1)
        i = 0
        j = N
        temp[0,N] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i+1,i+1,j-1,j-1)
        i = N
        j = N
        temp[N,N] = 2*u[n-1,i,j]-u[n-2,i,j]+lamda**2*laplacian(u,n-1,i,j,i-1,i-1,j-1,j-1)
        u[n] = temp
    return u

def plot_2d(u_0):
    '''
    plot the solution of the 2d wave equation
    '''
    u = solve_2d(u_0)
    init_state = u[0,:,:]#initial profile
    x = np.linspace(0, 1, 41)
    y = np.linspace(0, 1, 41)
    x, y = np.meshgrid(x, y)
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca(projection='3d')
    

    line1 = [ax.plot_surface(x,y,init_state, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')]#plot the initial profile using plot_surface function
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u');
    ax.set_zlim3d([-0.4,0.4])

    ax.set_title('$u_{tt} = u_{xx}+u_{yy} \, u_t(x,y,0) = 0$ homogeneous Neumann BCs');
    def update(step,u,line1): 
        state = u[step,:,:]#update at each timestep
        line1[0].remove()#remove previous plot
        ax.set_zlim3d([-0.4,0.4])
        line1[0] = ax.plot_surface(x,y,state, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')# plot surface at each timestep
        return line1,
    ani = animation.FuncAnimation(fig, update, frames = 266, interval= 50,fargs=(u, line1))
    plt.show()
    return ani




def solve_1d(u_0):
    '''
    solve 1d wave equation CTCS
    '''
    N = 100
    tmax = 1.5
    h = (1 - 0) / 100
    x = np.linspace(0, 1, 101)
    u0 = u_0(x)
    k = 0.009
    lamda = k/h
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*101).reshape(nsteps+1, 101)# create the grid in space-time, 2d array, first dim is timestep
    u[0] = u0
    u[0,0] = 0
    u[0,N] = 0
    temp = u0.copy()
    # First update solution at time step 1
    for j in range(1,N):
    	temp[j] = u[0,j]+0.5*lamda**2*(u[0,j+1]-2*u[0,j]+u[0,j-1])
    # update BC
    temp[0] = 0
    temp[N] = 0
    u[1] = temp

    #update rest of the timestep
    for i in range(2,int(nsteps)):
        for j in range(1, N):
            temp[j] = -u[i-2,j]+2*u[i-1,j]+lamda**2*(u[i-1,j+1]-2*u[i-1,j]+u[i-1,j-1])
        # update BCs
        temp[0] = 0
        temp[N] = 0
        u[i] = temp

    return u


def plot(u_0,pl="animation"):
    '''
    plot 1d wave equation
    pl can be animation,2d,or 3d
    '''
    u = solve_1d(u_0)
    
    init_state = u[0,:]
    if(pl=="animation"):
        fig, ax = plt.subplots()
        x = np.linspace(0, 1, 101)
        line1, = ax.plot(x, init_state)
        ax.set_xlim((0, 1))
        ax.set_ylim((-1,1))
        ax.set_ylabel("$u$")
        ax.set_title('$u_{tt} = u_{xx} \, u_t(x,0) = 0$ homogeneous Dirichlet BCs');
        def update(step):    
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
        plt.title('$u_{tt} = u_{xx} \, u_t(x,0) = 0$ homogeneous Dirichlet BCs')
        plt.show()
        return fig
    else:
        x = np.linspace(0, 1, 101)
        y = np.linspace(0,165, 166)
        x, y = np.meshgrid(x, y)
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca(projection='3d')
            
        ax.plot_surface(x,y, u[0:166,:], rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
            
        ax.set_title('$u_{tt} = u_{xx} \, u_t(x,0) = 0$ homogeneous Dirichlet BCs');
        ax.set_xlabel('x')
        ax.set_ylabel('t')
        ax.set_zlabel('u');
        plt.show()
        return fig
    


