import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtick

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


    ax.text2D(0.4,1,"%.1f year orbit" % (time[-1]),fontsize = 18,transform=ax.transAxes)
    #ax.scatter(x[0],y[0],z[0],"o",color=col)
    ax.plot(x,y,z,label="Orbit (%s)"%meth,color=col)
    ax.set_xlabel("x (AU)"); ax.set_ylabel("y (AU)"); ax.set_zlabel("z (AU)")

    ax2 = plt.subplot(2,1,1)
    hfont = {'fontname':'Times'}
    plt.title('Kinetic and Potential Energy', fontsize = 15,**hfont)
    ax2.plot(time,10**6*kin_en,"-",label="%s"%meth,color=col)
    ax2.set_ylabel("Kinetic Energy [$\mu$]",fontsize = 12,**hfont)
    ax2 = plt.subplot(2,1,2)
    ax2.plot(time,10**6*pot_en,"-",label="%s"%meth,color=col)

    #ax2.plot(time,kin_en+pot_en,label="Total energy (%s)"%meth,color=col)
    ax2.set_xlabel("t (year)",fontsize = 12); ax2.set_ylabel("Potential Energy [$\mu$]",fontsize = 12,**hfont)
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),ncol = 5,fancybox=True,fontsize=10)


fig = plt.figure(1)
ax = fig.gca(projection="3d")
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111); ax2.grid()
euler_plot = plot("Orbit_euler.txt","euler")
verlet_plot = plot("Orbit_verlet.txt","verlet")

ax.scatter(0,0,0,"O",label="Sun",color="orange",s=200)
ax.legend(); ax2.legend()
plt.show()
