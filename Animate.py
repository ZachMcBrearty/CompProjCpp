import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
fig.set_figheight(5)
fig.set_figwidth(5)
xs = []
ys = []
times = []
xdata = []
ydata = []
ln, = ax.plot([], [], 'o')
timetext = ax.text(0.05, 0.9,"", transform=ax.transAxes)
trace = []
scale = (1, "m")
timescale = (1, "s")
tracelength = 500

def getData(filename):
    times = [] # 1D
    kinetic = [] # 1D
    gravitational = [] # 1D
    total = [] # 1D

    objs = [] # 3D
    with open(filename, "r") as f:
        for line in f:
            time, energies, positions = line.strip().split(" # ")
            times.append(float(time))

            energies = list(map(float, energies.split(" , ")))
            kinetic.append(energies[0])
            gravitational.append(energies[1])
            total.append(energies[2])#

            pos = positions.split(" * ")
            pos = [list(map(float, q[1:-1].split(", "))) for q in pos]
            objs.append(pos)
    return times, kinetic, gravitational, total, objs

def getDataBinary(filename, p):
    '''p: expected number of particles'''
    n = 3 + 3 * p
    objs = []
    with open(filename, "rb") as f:
        data = np.fromfile(f)
        times = data[0::n]
        kinetic = data[1::n]
        gravi = data[2::n]
        total = kinetic + gravi
        data = np.delete(data, np.s_[0::n])
        data = np.delete(data, np.s_[0::n-1])
        data = np.delete(data, np.s_[0::n-2])
        objs = np.reshape(data, (-1, p, 3))
    return times, kinetic, gravi, total, objs

def init():
    global trace
    ax.set_xlim(-1.2, 1.2)
    ax.set_xlabel("x, " + scale[1])
    ax.set_ylim(-1.2, 1.2)
    ax.set_ylabel("y, " + scale[1])
    trace = [ax.plot([], [],'-', lw=1)[0] for _ in range(tracenum)]
    return ln, *trace, timetext,

def animate(i):  
    # x, y values to be plotted 
    i = int(i)
    xdata = xs[i]
    ydata = ys[i]
    ln.set_data(xdata, ydata)
    if i > tracelength:
        for j, t in enumerate(trace):
            t.set_data(xs[i-tracelength:i, j], ys[i-tracelength:i, j])
    else:
        for j, t in enumerate(trace):
            t.set_data(xs[:i, j], ys[:i, j])
    #timetext.set_text(str(times[i])+"s")
    timetext.set_text(str(round(times[i], 2))+ " " +timescale[1])
    return ln, *trace, timetext

def animateFile(filename, framesskip=1, repeat=True, scale_=(1, "m"), 
    p=2, timescale_=(1, "s")):
    global xs, ys, times, tracenum, scale, timescale
    timescale = timescale_
    scale = scale_
    if filename.endswith(".txt"):
        times, kin, grav, total, objs = getData(filename)
    elif filename.endswith(".bin"):
        times, kin, grav, total, objs = getDataBinary(filename, p)
    tracenum = len(objs[0])
    pos = np.array(objs)
    
    times = times / timescale[0]
    xs = pos[:, :, 0] / scale[0]
    ys = pos[:, :, 1] / scale[0]

    ani = FuncAnimation(fig, animate, frames=range(0, len(pos), framesskip),
                        init_func=init, blit=True, repeat=repeat)
    
    plt.show()

def graphEnergies(filename, change=False, p=2):
    fig, ax = plt.subplots()
    if filename.endswith(".txt"):
        times, kin, grav, total, objs = getData(filename)
    elif filename.endswith(".bin"):
        times, kin, grav, total, objs = getDataBinary(filename, p)
    t = times / timescale[0]
    if not change:
        ax.plot(t, kin, label="KE")
        ax.plot(t, grav, label="GPE")
        ax.plot(t, total, label="Total Energy")
        ax.set_ylabel("Energy, J")
        ax.set_xlabel("Time, " + timescale[1])
    else:
        te = np.array(total)
        initenergy = te[0]
        change = te - initenergy
        percchange = 100 * change / initenergy
        ax.plot(t, percchange, label="PercentChange in Total")
        ax.set_ylabel("Percent change in Energy, %")
        ax.set_xlabel("Time, "+timescale[1])
    plt.legend()
    plt.show()


if __name__=="__main__":
    f = "PlanetsSim.txt"
    fbin = "GravitySim.bin"
    animateFile(fbin, framesskip=1, repeat=False, scale_=10, p=2)
    #graphEnergies(f)
    #graphEnergies(f, True)