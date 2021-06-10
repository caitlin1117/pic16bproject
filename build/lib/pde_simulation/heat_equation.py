'''
1d heat equation
ut=uxx x in (0,L)  u(x,0) = I(x) u(0,t) = u(L,t) = 0 
'''
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
    if dim==1:
        p = plot(u_0,pl,method)
        if save:
            if pl=="animation":
                 p.save('animation.mp4')
            else:
                 p.savefig("plot.png")
        return p 
    else:
        if pl=="heatmap":
            p = plot_2d(u_0)
        else:
            p = plot_3d(u_0)
        if save:
            p.save('animation.mp4')
        return p




def solve_1d_FTCS(u_0):
    '''
    FTCS
    xi tn
    u[i,n+1] = (1-2lamda)u[i,n]+lamda(u[i+1,n]+u[ii1,n] interior points
    u[0,j] = u[M,j] = 0 BC
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
    #print(np.shape(u))
    u[0] = u0
    u[0,0] = 0
    u[0,N] = 0
    temp = u0.copy()
    

    for i in range(1,int(nsteps)):
        for j in range(1, N):
            temp[j] = (1-2*lamda)*u[i-1,j]+lamda*(u[i-1,j+1]+u[i-1,j-1])

        temp[0] = 0
        temp[N] = 0
        u[i] = temp

    return u




def plot(u_0,pl ="animation",method ="implicit"):
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
        def update(step):    
            state = u[step,:]
            line1.set_data(x, state)
            return
        ani = animation.FuncAnimation(fig, update, frames = 100, interval= 200)
        #plt.show()
        return ani
    elif pl=="2d":
        fig = plt.figure(figsize=(9,9))
        plt.imshow(u[99::-1, :], extent=[0,1,0,1])
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('t')
        plt.title('$u_t = u_{xx}$  homogeneous Dirichlet BCs %s method'%(method,))
        #plt.show()
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
        #plt.show()
        return fig


def solve_1d_BTCS(u_0):
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
    #print(np.shape(u))
    u[0] = u0
    u[0,0] = 0
    u[0,N] = 0
    temp = u0.copy()
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
        temp[1: N] = np.linalg.solve(A,B)
        temp[0] = 0
        temp[N] = 0
        u[i] = temp
    return u


def solve_2d(u_0):
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

    ic01 = np.logical_and(x >= 1/4, x <= 3/4)
    ic02 = np.logical_and(y >= 1/4, y <= 3/4)
    ic0 = np.multiply(ic01, ic02)
     
    #u_initial = np.random.uniform(low=0, high=1, size=(1,1))

    u[0] = ic0*1
    temp = u0.copy()
    


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
        sol = np.linalg.solve(A, np.matmul(A_rhs, b1))
        temp = np.flipud(sol).reshape(N, N)
            
        b2 = np.flipud(temp).reshape(N**2, 1)
        sol = np.linalg.solve(A, np.matmul(A_rhs, b2))
        temp = np.flipud(sol).reshape(N, N)
        u[n] = temp
        

    return u

def plot_3d(u_0):
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
    def update(step,u,line1): 
        state = u[step,:,:]
        line1[0].remove()
        ax.set_zlim3d([0, 1])
        line1[0] = ax.plot_surface(x,y,state, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
        return line1,
    ani = animation.FuncAnimation(fig, update, frames = 166, interval= 100,fargs=(u, line1))
    #plt.show()
    return ani




def plot_2d(u_0):
    u = solve_2d(u_0)
    def update(step):
        plt.clf()
        plt.title("$u_t = u_{xx}+u_{yy}$ homogeneous Neumann BCs")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.pcolormesh(u[step], cmap=plt.cm.jet, vmin=0, vmax=1)
        plt.colorbar()

        return plt
    ani = animation.FuncAnimation(plt.figure(), update, frames = 166, interval= 100, repeat = False)
    #plt.show()
    return ani





def f(x):
	return 4*x-4*x**2


def Gaussian(x):
    return np.exp(-200*(x-0.5)**2)
def Gaussian_2d(x,y):
    return 0.8*np.exp(-70*((x-0.5)**2+(y-0.5)**2))

#solution(Gaussian,1,pl="3d", method="explicit")





