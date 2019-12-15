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


def system_plot(file):
    planet_names = open(file,"r")
    header = planet_names.readline()
    names = header.split()
    data = np.loadtxt(file,skiprows=1)
    sorted_data = np.transpose(data)
    steps = len(data)
    planets = int(len(sorted_data)/3)

    """"check if planets are in within actual orbits """
    AU = 149597870700
    Perehelion_values = (np.array([46.0,107.5,147.1,206.6,740.5,1352.6,2741.3,4444.5])*10**9)/AU
    for i in range(planets-1):
        x = sorted_data[3*i]
        y = sorted_data[3*i+1]
        z = sorted_data[3*i+2]
        r = np.sqrt(x**2+y**2+z**2)
        perehelion_model = r[np.argmin(r)]
        relative_error = abs(perehelion_model-Perehelion_values[i])/Perehelion_values[i]
        print(relative_error)

    system_fig = plt.figure(3)
    system_ax = system_fig.gca(projection="3d")
    system_ax.scatter(0,0,0,"O",label="Sun",color="orange",s=50)
    system_ax.autoscale(enable=False,axis='both')  #you will need this line to change the Z-axis
    system_ax.set_xbound(-5, 5)#system_ax.set_xbound(-0.001, 0.001)
    system_ax.set_ybound(-5, 5)#system_ax.set_ybound(-0.001, 0.001)
    system_ax.set_zbound(-2, 2)#system_ax.set_zbound(-0.0001, 0.0001)


    #for i in range(planets):
    #    system_ax.plot(sorted_data[3*i],sorted_data[3*i+1],sorted_data[3*i+2],label=names[i])
    #system_ax.set_xlabel("x (AU)"); system_ax.set_ylabel("y (AU)"); system_ax.set_zlabel("z (AU)")
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0),ncol = 4,fancybox=True,fontsize=10)





"""
fig = plt.figure(1)
ax = fig.gca(projection="3d")
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111); ax2.grid()
#euler_plot = plot("Orbit_euler.txt","euler")
verlet_plot = plot("Orbit_verlet.txt","verlet")

ax.scatter(0,0,0,"O",label="Sun",color="orange",s=200)
ax.legend(); ax2.legend()
plt.show()
"""



#plt.title('Three body problem')#,loc = 'upper center', bbox_to_anchor=(0.5, 1))
system_plot("data_orbits.txt")
#plt.show()
