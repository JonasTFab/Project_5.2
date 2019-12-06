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

    """ax2 = plt.subplot(2,1,1)
    hfont = {'fontname':'Times'}
    plt.title('Kinetic and Potential Energy', fontsize = 15,**hfont)
    ax2.plot(time,10**6*kin_en,"-",label="%s"%meth,color=col)
    ax2.set_ylabel("Kinetic Energy [$\mu$]",fontsize = 12,**hfont)
    ax2 = plt.subplot(2,1,2)
    ax2.plot(time,10**6*pot_en,"-",label="%s"%meth,color=col)

    #ax2.plot(time,kin_en+pot_en,label="Total energy (%s)"%meth,color=col)
    ax2.set_xlabel("t (year)",fontsize = 12); ax2.set_ylabel("Potential Energy [$\mu$]",fontsize = 12,**hfont)
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),ncol = 5,fancybox=True,fontsize=10)"""


"""fig = plt.figure(1)
ax = fig.gca(projection="3d")
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111); ax2.grid()
#euler_plot = plot("Orbit_euler.txt","euler")
verlet_plot = plot("Orbit_verlet.txt","verlet")

ax.scatter(0,0,0,"O",label="Sun",color="orange",s=200)
ax.legend(); ax2.legend()
plt.show()"""


def plot_d(filename):
    data = np.loadtxt(filename)
    sorted_data = np.transpose(data)
    colors = ["r","g","b","c","magenta"]
    cut = int(len(sorted_data[0])/5)
    for i in range(0,5):
        f_index = i*cut
        l_index = (i+1)*cut - 1
        x = sorted_data[0][f_index:l_index]
        y = sorted_data[1][f_index:l_index]
        z = sorted_data[2][f_index:l_index]
        time = sorted_data[3][f_index:l_index]
        kin_en = sorted_data[4][f_index:l_index]
        pot_en = sorted_data[5][f_index:l_index]
        col = colors[i]
        ax3.plot(x,y,z,label="Orbit",color=col)
    ax3.set_xlabel("x (AU)"); ax3.set_ylabel("y (AU)"); ax3.set_zlabel("z (AU)")
    ax3.scatter(0,0,0,"O",label="Sun",color="orange",s=200)
    ax3.legend();
    plt.show()

fig3 = plt.figure(3)
ax3 = fig3.gca(projection="3d")
diff_r_orbits = plot_d("Orbit_diff_r_exp.txt")
