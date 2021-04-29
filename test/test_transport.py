import sys
sys.path.insert(1, '/Users/caitlin/Documents/Math/PIC16B/pic16bproject')
from pde_simulation import transport_equation as te
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['animation.ffmpeg_path'] = '/Users/caitlin/opt/anaconda3/bin/ffmpeg'
import matplotlib.animation as animation
from matplotlib import colors
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
def f(x):
    return np.exp(-200*(x-0.25)**2)
te.solution(f,"animation")
