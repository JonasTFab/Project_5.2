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

    fig = plt.figure(1)
    ax = fig.gca(projection="3d")
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111); ax2.grid()

    ax.text2D(0.4,1,"Orbit in %.1f year" % (time[-1]),transform=ax.transAxes)
    ax.scatter(x[0],y[0],z[0],"o",label="Earth (%s)"%meth,color=col)
    ax.plot(x[::100],y[::100],z[::100],label="Earth path (%s)"%meth,color=col)
    ax.set_xlabel("x (AU)"); ax.set_ylabel("y (AU)"); ax.set_zlabel("z (AU)")

    if (file=="mercury_perihelion.txt"):
        ax.autoscale(enable=False,axis='both')
        center = 0.3075
        diff = center/10
        ax.set_xbound(center-diff, center+diff)
        ax.set_ybound(-diff, diff)
        ax.set_zbound(-diff, diff)

    ax2.plot(time[::10],kin_en[::10],"--",label="Kinetic energy (%s)"%meth,color=col)
    ax2.plot(time[::10],pot_en[::10],"-.",label="Potential energy (%s)"%meth,color=col)
    ax2.plot(time,kin_en+pot_en,label="Total energy (%s)"%meth,color=col)
    ax2.set_xlabel("t (year)"); ax2.set_ylabel("Energy")
    ax.scatter(0,0,0,"O",label="Sun",color="orange",s=200)
    ax.legend(); ax2.legend()

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
        aphelion_model = r[np.argmax(r)]
        print(aphelion_model-perehelion_model)
        print(aphelion_model, perehelion_model)
        relative_error = abs(perehelion_model-Perehelion_values[i])/Perehelion_values[i]
        #print(relative_error)

    system_fig = plt.figure(3)
    system_ax = system_fig.gca(projection="3d")
    system_ax.scatter(0,0,0,"O",label="Sun",color="orange",s=50)
    system_ax.autoscale(enable=False,axis='both')  #you will need this line to change the Z-axis
    system_ax.set_xbound(-5, 5)#system_ax.set_xbound(-0.001, 0.001)
    system_ax.set_ybound(-5, 5)#system_ax.set_ybound(-0.001, 0.001)
    system_ax.set_zbound(-0.5, 0.5)#system_ax.set_zbound(-0.0001, 0.0001)

    for i in range(planets):
        system_ax.plot(sorted_data[3*i],sorted_data[3*i+1],sorted_data[3*i+2],label=names[i])
    system_ax.set_xlabel("x (AU)"); system_ax.set_ylabel("y (AU)"); system_ax.set_zlabel("z (AU)")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0),ncol = 4,fancybox=True,fontsize=10)

def data_analysis(filename):
    data = np.loadtxt(filename)
    sorted_data = np.transpose(data)

    x = sorted_data[0]
    y = sorted_data[1]

    #z = sorted_data[2]
    #time = sorted_data[3]
    #kin_en = sorted_data[4]
    #pot_en = sorted_data[5]
    mercury_years = int(len(x)/(365/88))
    last_x = x[-mercury_years-1000:-1] #Slice data as to only use the last year
    last_y = y[-mercury_years-1000:-1]
    #last_z = z[-mercury_years-1000:-1]

    r = np.sqrt(last_x**2 + last_y**2)
    perehelion_idx = np.argmin(r)
    last_perehelion = r[perehelion_idx]
    first_perehelion = 0.3075
    perehelion_tangens = (last_y[perehelion_idx-54]/last_x[perehelion_idx-54])
    perehelion_tangens = (last_y[perehelion_idx]/last_x[perehelion_idx])
    theta = np.arctan(perehelion_tangens)
    arcsec = np.rad2deg(theta)*3600
    print(arcsec)
    plt.plot(last_x,last_y)
    plt.plot(last_x[perehelion_idx-54],last_y[perehelion_idx-54],'o',color = 'g')
    plt.plot(last_x[perehelion_idx],last_y[perehelion_idx],'o',color = 'r')
    plt.show()

def energy(filename,subplot,jup_mass):
    data = np.loadtxt(filename)
    sorted_data = np.transpose(data)
    steps = len(data)
    planets = int(len(sorted_data)/2)

    x = np.linspace(1,steps,steps)*25/steps

    total_energy = np.zeros((2,steps))
    total_energy[0,:] = sorted_data[0]+sorted_data[2]
    total_energy[1,:] = sorted_data[1]+sorted_data[3]
    E = total_energy[0,:]+total_energy[1,:]
    E0 = E[0]

    size=25
    plt.figure(7); plt.subplot(subplot)
    plt.plot(x,E/E0,label="Jupiter+Earth")
    plt.xlabel("Time (yr)",fontsize=size); plt.ylabel(r"$E/E_0$ ($[m_\odot AU^2/yr^2]$)",fontsize=size)
    plt.title(r"$10^%i$ timesteps, Jupiter mass = %.i $m_j$" % (np.log10(steps),jup_mass),fontsize=size)
    mean_tot = sum(total_energy[0,:]+total_energy[1,:])/steps
    plt.xlim(-0.5,25.5,size); plt.ylim(0.999999,np.amax(E/E0))
    plt.tick_params(axis='both', which='both', labelsize=size/1.5)
    #plt.tick_params(axis='both', which='minor', labelsize=size/2)
    plt.grid(); plt.legend()


#plot("Orbit_euler.txt","euler")
#plot("Orbit_verlet.txt","verlet")
#plot("mercury_perihelion.txt","verlet")
#plt.show()

system_plot("data_orbits.txt")
plt.title('Three body problem')#,loc = 'upper center', bbox_to_anchor=(0.5, 1))
#plt.show()


#data_analysis("mercury_perihelion.txt")

# Works only for Earth-Jupiter-fixed sun system
# at total time = 25 years
energy("energy_1m_jupiter.txt",121,1)
energy("energy_1000m_jupiter.txt",122,1000)
plt.show()
