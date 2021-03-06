import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
plt.rcParams['animation.ffmpeg_path'] = '/Users/caitlin/opt/anaconda3/bin/ffmpeg'
import matplotlib.animation as animation
from matplotlib import colors
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from IPython.display import HTML
plt.style.use('seaborn-pastel')
from mpl_toolkits import mplot3d
from matplotlib import cm


def solution(u_0, dim, pl="animation", method="implicit",save=False):
    '''
    The 1d heat equation u_t = u_xx in [0,1] with homogeneous Dirichlet BCs can be solved 
    using Forward Time Central Space (explicit) or Backward Time Central Space (implicit) methods 
    The 2d heat equation u_t = u_xx+u_yy in [0,1]x[0,1] with homogeneous Neumann BCs is solved using the ADI method, 
    Input: (u_0, dim, pl, method, save) 
    u_0: u_0 is the initial profile, 
    dim: can be specified to 1 or 2 indicating you want to solve a 1d or 2d heat equation
    pl: for 1d heat equation, pl can be specified to animation, 2d or 3d for different types of plot, default is animation
        for 2d heat equation, pl can only be specified to animation, default is animation
    method: for 1d heat equation, the method can be specified to implicit or explicit, default is implicit
    save: specified to True if you want the animation(mp4 file) or figure(png file) saved, default is false
    Output:animation or figure of the solution depending on different kind of plot
    '''
    if dim==1:
        p = plot(u_0,pl,method)
        if save:
            if pl=="animation":
                 p.save('animation.mp4')
            else:
                 p.savefig("plot.png")
        return p 
    else:
        p = plot_3d(u_0)
        if save:
            p.save('animation.mp4')
        return p



def solve_1d_FTCS(u_0):
    '''
    FTCS
    xi tn
    u[i,n+1] = (1-2lamda)u[i,n]+lamda(u[i+1,n]+u[i-1,n] interior points
    u[0,j] = u[M,j] = 0 boundary points
    '''
    N = 20
    tmax = 0.1
    h = (1 - 0) / N
    x = np.linspace(0, 1, N+1)
    u0 = u_0(x)
    k = 0.001
    lamda = k/h**2
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*(N+1)).reshape(nsteps+1, (N+1))
    u[0] = u0
    u[0,0] = 0
    u[0,N] = 0
    temp = u0.copy()
    # update each time step
    for i in range(1,int(nsteps)):
        for j in range(1, N):
            temp[j] = (1-2*lamda)*u[i-1,j]+lamda*(u[i-1,j+1]+u[i-1,j-1])
        # update boundary points
        temp[0] = 0
        temp[N] = 0
        u[i] = temp

    return u




def plot(u_0,pl ="animation",method ="implicit"):
    '''
    plot 1d heat equation
    '''
    if method =="implicit":
        u = solve_1d_BTCS(u_0)
        N=100
    else:
        u = solve_1d_FTCS(u_0)
        N=20

    init_state = u[0,:]
    if(pl=="animation"):
        fig, ax = plt.subplots()
        x = np.linspace(0, 1, N+1)
        line1, = ax.plot(x, init_state)
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))
        ax.set_ylabel("$u$")
        ax.set_title('$u_t = u_{xx}$ homogeneous Dirichlet BCs %s method'%(method,));
        def update(step):  #update each time step    
            state = u[step,:]
            line1.set_data(x, state)
            return
        ani = animation.FuncAnimation(fig, update, frames = 100, interval= 200)
        plt.show()
        return ani
    elif pl=="2d":
        fig = plt.figure(figsize=(9,9))
        plt.imshow(u[99::-1, :], extent=[0,1,0,1])
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('t')
        plt.title('$u_t = u_{xx}$  homogeneous Dirichlet BCs %s method'%(method,))
        plt.show()
        return fig
    else:
        x = np.linspace(0, 1, N+1)
        y = np.linspace(0,0.99, 100)
        x, y = np.meshgrid(x, y)
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca(projection='3d')
            
        ax.plot_surface(x,y, u[0:100,:], rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
            
        ax.set_title('$u_t = u_{xx}$  homogeneous Dirichlet BCs %s method'%(method,))
        ax.set_xlabel('x')
        ax.set_ylabel('t')
        ax.set_zlabel('u');
        plt.show()
        return fig


def solve_1d_BTCS(u_0):
    '''
    solve 1d heat equation using BTCS method
    solve using tri-diagonal matrices
    unconditionally stable
    '''
    N = 100
    tmax = 0.1
    h = (1 - 0) / N
    x = np.linspace(0, 1, N+1)
    u0 = u_0(x)
    k = 0.001
    lamda = k/h**2
    a = 1+2*lamda
    b = -lamda
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*(N+1)).reshape(nsteps+1, (N+1))
    u[0] = u0
    u[0,0] = 0
    u[0,N] = 0
    temp = u0.copy()
    # create tridiagonal matrix
    A = np.zeros((N-1)*(N-1)).reshape((N-1),(N-1))
    A[0,0] = a
    A[0,1] = b
    for i in range (1, N-2):
        A[i,i-1] = b
        A[i,i] = a
        A[i,i+1] = b
    A[N-2,N-3] = b
    A[N-2,N-2] = a

    for i in range(1,int(nsteps)):
        B = u[i-1,1:N]
        temp[1: N] = np.linalg.solve(A,B)#solve matrix equation Ax=B
        #update boundary points
        temp[0] = 0
        temp[N] = 0
        u[i] = temp
    return u


def solve_2d(u_0):
    '''
    Alternating-direction implicit method
    Solve heat equation using tri-diagonal matrices
    Each timestep delta t is subdivided into 2 steps
    The space derivatives are approximated implicited in the x-direction and explicited in the y-direction 
    in the first half time step. The procedure is reversed over the second one
    '''
    N = 50
    tmax = 2.5
    h = (1 - 0) / (N-1)
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    x, y = np.meshgrid(x, y)
    u0 = u_0(x,y)
    k = 0.015
    lamda = 0.01*k/(2*h**2)
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*(N)*(N)).reshape(nsteps+1, (N), (N))
    u[0] = u0
    temp = u0.copy()
    
    #Create left matrix
    maindiag = (1+2*lamda)*np.ones((1, N))
    offdiag = -lamda*np.ones((1, N-1))
    a = maindiag.shape[1]
    diagonals = [maindiag, offdiag, offdiag]
    Lx = sparse.diags(diagonals, [0, -1, 1], shape=(a, a)).toarray()
    Ix = sparse.identity(N).toarray()
    A = sparse.kron(Ix, Lx).toarray()
        
    pos1 = np.arange(0,(N)**2,N)
        
    for i in range(len(pos1)):
        A[pos1[i], pos1[i]] = 1 + lamda
            
    pos2 = np.arange(N-1, N**2, N)
        
    for j in range(len(pos2)):
        A[pos2[j], pos2[j]] = 1 + lamda

    #Create right matrix
    maindiag = (1-lamda)*np.ones((1, N))
    offdiag = lamda*np.ones((1, N-1))
    a = maindiag.shape[1]
    diagonals = [maindiag, offdiag, offdiag]
    Rx = sparse.diags(diagonals, [0, -1, 1], shape=(a, a)).toarray()
    Ix = sparse.identity(N).toarray()
    A_rhs = sparse.kron(Rx, Ix).toarray()
        
    pos3 = np.arange(N, N**2-N)
        
    for k in range(len(pos3)):
        A_rhs[pos3[k], pos3[k]] = 1 - 2*lamda

    
    for n in range(1,nsteps):
        
        b1 = np.flipud(u[n-1]).reshape(N**2, 1)
        sol = np.linalg.solve(A, np.matmul(A_rhs, b1))#solve Ax=A_rhs*b1
        temp = np.flipud(sol).reshape(N, N)
            
        b2 = np.flipud(temp).reshape(N**2, 1)
        sol = np.linalg.solve(A, np.matmul(A_rhs, b2))#solve Ax=A_rhs*b2
        temp = np.flipud(sol).reshape(N, N)
        u[n] = temp
        

    return u

def plot_3d(u_0):
    '''
    plot solution of the 2d heat equation
    '''
    u = solve_2d(u_0)
    init_state = u[0,:,:]
    x = np.linspace(0, 1, 50)
    y = np.linspace(0, 1, 50)
    x, y = np.meshgrid(x, y)
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca(projection='3d')
    

    line1 = [ax.plot_surface(x,y,init_state, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')]
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u');
    ax.set_zlim3d([0, 1])

    ax.set_title('$u_t = u_{xx}+u_{yy}$ homogeneous Neumann BCs ADI method');
    def update(step,u,line1): #update each time step  
        state = u[step,:,:]
        line1[0].remove()
        ax.set_zlim3d([0, 1])
        line1[0] = ax.plot_surface(x,y,state, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
        return line1,
    ani = animation.FuncAnimation(fig, update, frames = 166, interval= 100,fargs=(u, line1))
    plt.show()
    return ani






