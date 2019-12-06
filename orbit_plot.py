import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot(file,method):
    data = np.loadtxt(file)
    sorted_data = np.transpose(data)

    x = sorted_data[0]
    y = sorted_data[1]
    z = sorted_data[2]
    time = sorted_data[3]
    kin_en = sorted_data[4]
    pot_en = sorted_data[5]

    if method == "euler":
        col = "black"
        meth = "Euler"
    elif method == "verlet":
        col = "red"
        meth = "Verlet"


    ax.text2D(0.4,1,"Orbit in %.1f year" % (time[-1]),transform=ax.transAxes)
    ax.scatter(x[0],y[0],z[0],"o",label="Earth (%s)"%meth,color=col)
    ax.plot(x,y,z,label="Earth path (%s)"%meth,color=col)
    ax.set_xlabel("x (AU)"); ax.set_ylabel("y (AU)"); ax.set_zlabel("z (AU)")

    ax2.plot(time,kin_en,"--",label="Kinetic energy (%s)"%meth,color=col)
    ax2.plot(time,pot_en,"-.",label="Potential energy (%s)"%meth,color=col)
    ax2.plot(time,kin_en+pot_en,label="Total energy (%s)"%meth,color=col)
    ax2.set_xlabel("t (year)"); ax2.set_ylabel("Energy")


fig = plt.figure(1)
ax = fig.gca(projection="3d")
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111); ax2.grid()
#euler_plot = plot("Orbit_euler.txt","euler")
verlet_plot = plot("Orbit_verlet.txt","verlet")

ax.scatter(0,0,0,"O",label="Sun",color="orange",s=200)
ax.legend(); ax2.legend()
plt.show()
