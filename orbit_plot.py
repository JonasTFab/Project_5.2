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
        col = "red"

    plt.figure(1); plt.grid()
    plt.title("Orbit in %.1f year" % (time[-1]))
    plt.plot(0,0,"o",label="Sun",color="yellow")
    plt.plot(x[0],y[0],"o",label="Earth",color="blue")
    plt.plot(x,y,label="Earth orbit",color="blue")
    plt.xlabel("x (AU)"); plt.ylabel("y (AU)")
    plt.legend()

    plt.figure(2); plt.grid()
    plt.plot(time,kin_en,label="Kinetic energy")
    plt.plot(time,pot_en,label="Potential energy")
    plt.plot(time,kin_en+pot_en,label="Total energy")
    plt.xlabel("t (year)"); plt.ylabel("Energy")
    plt.legend()


euler_plot = plot("Orbit_euler.txt","euler")
verlet_plot = plot("Orbit_verlet.txt","verlet")

plt.show()
