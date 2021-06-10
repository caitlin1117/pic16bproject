"""
Finite-difference solver for transport equation:
    u_t + uu_x = 0.
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

def solution(type,save=False):
    #u = solve_wendroff(type)
    p = plot(type)
    if save:
        p.save('animation.mp4')       
    return p


def f1(x):
    return np.piecewise(x, [x < 0, x >= 0], [1, 0])
def f2(x):
    return np.piecewise(x, [x < 0, x >= 0], [0, 1])


def solve_wendroff(type):
    if type =="shock":
        u_0 = f1
    else:
        u_0 = f2
    N = 400
    tmax = 1.5
    h = (4 - 0) / 400
    x = np.linspace(-2-h, 2+h, 403)
    u0 = u_0(x)
    k = 0.009
    lamda = k/h
    nsteps = int(tmax/k)
    u = np.ones((nsteps+1)*403).reshape(nsteps+1, 403)
    #print(np.shape(u))
    u[0] = u0
    temp = u0.copy()
    for i in range(1,int(nsteps)):
        for j in range(N+2):
            temp[j] = u[i-1,j]-lamda*(0.5*u[i-1,j]**2-0.5*u[i-1,j-1]**2)
        u[i] = temp
        u[i, 0] = u[i, N+1]
        u[i, N+2] = u[i,1]

    return u


def plot(type):
    #u = solve(u_0)
    u = solve_wendroff(type)
    #e = exact(u_0)
    init_state = u[0,99:302]
    
    fig, ax = plt.subplots()
    x = np.linspace(-1.01, 1+0.01, 203)
    line1, = ax.plot(x, init_state)
    #line2, = ax.plot(x, init_state) 
    ax.set_xlim((-1, 1))
    #ax.set_ylim((0, 0.6)) 
    ax.set_ylabel("$u$")
    ax.set_title('$u_t+uu_x=0$ %s wave'%(type,));
    def update(step):    
        state = u[step,99:302]
        #state2 = e[step,:]#u_0(x-c*step*0.009)
        line1.set_data(x, state)
        #line2.set_data(x, state2)
        return
    ani = animation.FuncAnimation(fig, update, frames = 166, interval=100)
    #plt.show()
    return ani
    
    



#solution("rarefaction")
 

