import matplotlib.pyplot as plt
import numpy as np


def plot(file,method):
    data = np.loadtxt(file)
    sorted_data = np.transpose(data)

    x = sorted_data[0]
    y = sorted_data[1]
    time = sorted_data[2]
    kin_en = sorted_data[3]
    pot_en = sorted_data[4]

    if method == "euler":
        col = "black"
        meth = "Euler"
    elif method == "verlet":
        col = "red"
        meth = "Verlet"


    plt.figure(1)
    plt.title("Orbit in %.1f year" % (time[-1]))
    plt.plot(x[0],y[0],"o",label="Earth (%s)"%meth,color=col)
    plt.plot(x,y,label="Earth orbit (%s)"%meth,color=col)
    plt.xlabel("x (AU)"); plt.ylabel("y (AU)")
    plt.legend()

    plt.figure(2)
    plt.plot(time,kin_en,"--",label="Kinetic energy (%s)"%meth,color=col)
    plt.plot(time,pot_en,"-.",label="Potential energy (%s)"%meth,color=col)
    plt.plot(time,kin_en+pot_en,label="Total energy (%s)"%meth,color=col)
    plt.xlabel("t (year)"); plt.ylabel("Energy")
    plt.legend()


euler_plot = plot("Orbit_euler.txt","euler")
verlet_plot = plot("Orbit_verlet.txt","verlet")

plt.figure(1); plt.grid()
plt.plot(0,0,"o",label="Sun",color="yellow"); plt.legend()
plt.figure(2); plt.grid()
plt.show()
