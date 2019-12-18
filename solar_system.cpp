#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include <sstream>
#include <string>
#include <list>
//#include "mpi.h"

//For debugging:
// compile with: g++ solar_system.cpp -o main1 -larmadillo -llapack -lblas
// execute with ./main1 "Method(euler or verlet)" "Number of integration points" "Initial velocity" "Time period in years"

//For paralellization
//mpicxx  -o main_mpi.x  main.cpp -std=c++11
//mpiexec -n 2 ./main_mpi.x

std::ofstream ofile;
std::ofstream enfile;
double const pi = 3.14159265359;
double const AU = 149597870700;               // astronomical unit in meters
double const c = 299792458.0*60*60*24*365/AU; // speed of light in vaccum AU/year
double const eps = 1e-5;                      // small error

//planetary properties
double const GM = 4*pi*pi;                    // gravitational constant times one solar mass
double const sun_rad = 0.00465047;            // in AU
double const M_sun = 2*pow(10,30);            // all masses are given in kilograms
double const M_earth = 6*pow(10,24);
double const M_jupiter = 1.9*pow(10,27);
double const M_mercury = 3.302*pow(10,23);
double const M_venus = 48.685*pow(10,23);
double const M_mars = 6.4171*pow(10,23);
double const M_saturn = 5.6834*pow(10,26);
double const M_neptune = 102.413*pow(10,24);
double const M_uranus = 8.8*pow(10,25);


class object {
private:
  double a, r, ax, ay, az, dt;
  double ax_new, ay_new, az_new, ax_prev, ay_prev, az_prev;
  double vx_half, vy_half, vz_half;
  double merc_a = 0;
  int merc_peri = 0;
public:
  double mass;

  int N;
  double tmax = 1;
  arma::Col <double> x;
  arma::Col <double> y;
  arma::Col <double> z;
  arma::Col <double> vx;
  arma::Col <double> vy;
  arma::Col <double> vz;
  arma::Col <double> kin_en;
  arma::Col <double> pot_en;

  // initialize
  object(double x0,double y0,double z0,double vx0,double vy0,double vz0,double M, int len)
  {
    N = len;
    x = arma::zeros(N);
    y = arma::zeros(N);
    z = arma::zeros(N);
    vx = arma::zeros(N);
    vy = arma::zeros(N);
    vz = arma::zeros(N);
    kin_en = arma::zeros(N);
    pot_en = arma::zeros(N);

    mass = M/M_sun;
    x(0) = x0;
    y(0) = y0;
    z(0) = z0;
    vx(0) = vx0*365;
    vy(0) = vy0*365;
    vz(0) = vz0*365;

  }

  // functions
  void mercury_perihelion_initialize()
  {
    merc_peri = 1;
    x(0) = 0.3075;
    y(0) = 0;
    z(0) = 0;
    vx(0) = 0;
    vy(0) = 12.44;
    vz(0) = 0;
  }

  void mercury_perihelion(double distance, double merc_x, double merc_y, double merc_z,
                                           double merc_vx, double merc_vy, double merc_vz)
  {
    double I,J,K,l;
    I = merc_y*merc_vz - merc_vy*merc_z;
    J = merc_x*merc_vz - merc_vx*merc_z;
    K = merc_x*merc_vy - merc_vx*merc_y;

    l = sqrt(I*I + J*J + K*K);
    //merc_a = 3*GM*l*l / (c*c*pow(distance,5));
    //merc_a = GM*7.33*pow(10,-8)/ pow(distance,3);
    merc_a = 7.33*pow(10,-8);
    //merc_a = GM*7.33*10e-8/ pow(distance,3);
    std::cout << merc_a << std::endl;
    //std::cout << 3*GM*l*l / (c*c*pow(distance,5)) << std::endl;
  }

  void velocity_verlet()
  {
    dt = tmax/N;
    r = sqrt(x(0)*x(0) + y(0)*y(0) + z(0)*z(0));
    if (merc_peri == 1){mercury_perihelion(r,x(0),y(0),z(0),vx(0),vy(0),vz(0));}      // adding general relativity to Mercury
    a = GM / (r*r*r)*(1 + merc_a);
    ax_prev = -a * x(0);
    ay_prev = -a * y(0);
    az_prev = -a * z(0);
    for (int i=1; i<N; i++){
      x(i) = x(i-1) + dt*vx(i-1) + 0.5*dt*dt*ax_prev;
      y(i) = y(i-1) + dt*vy(i-1) + 0.5*dt*dt*ay_prev;
      z(i) = z(i-1) + dt*vz(i-1) + 0.5*dt*dt*az_prev;

      r = sqrt(x(i)*x(i) + y(i)*y(i) + z(i)*z(i));
      a = GM / (r*r*r)+merc_a;
      ax_new = -a * x(i);
      ay_new = -a * y(i);
      az_new = -a * z(i);

      vx(i) = vx(i-1) + 0.5*dt*(ax_new + ax_prev);
      vy(i) = vy(i-1) + 0.5*dt*(ay_new + ay_prev);
      vz(i) = vz(i-1) + 0.5*dt*(az_new + az_prev);
      //if (merc_peri == 1){mercury_perihelion(r,x(i),y(i),z(i),vx(i),vy(i),vz(i));}
      //std::cout << merc_a << std::endl;
      ax_prev = ax_new;
      ay_prev = ay_new;
      az_prev = az_new;
    }
  }

  void euler()
  {
    dt = tmax/N;
    for (int i=1; i<N; i++){
      r = sqrt(x(i-1)*x(i-1) + y(i-1)*y(i-1) + z(i-1)*z(i-1));
      a = GM / (r*r*r);
      ax = -a*x(i-1);
      ay = -a*y(i-1);
      az = -a*z(i-1);
      vx(i) = vx(i-1) + ax*dt;
      vy(i) = vy(i-1) + ay*dt;
      vz(i) = vz(i-1) + az*dt;
      x(i) = x(i-1) + vx(i-1)*dt;
      y(i) = y(i-1) + vy(i-1)*dt;
      z(i) = z(i-1) + vz(i-1)*dt;
    }
  }

  void kinetic_energy()
  {
    for (int i; i<N; i++){
      kin_en(i) = 0.5*mass*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i));
    }
  }

  void potential_energy()
  {
    for (int i; i<N; i++){
      pot_en(i) = -GM*mass/sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i));
    }
  }


  void write_to_file(std::string filename)
  {
    dt = tmax/N;
    ofile.open(filename);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    int start = N/100*99;
    for (int i=start; i<N; i++){
      ofile << std::setw(15) << x(i);
      ofile << std::setw(15) << y(i) << "\n";
      ofile << std::setw(15) << z(i);
      ofile << std::setw(15) << dt*i;
      ofile << std::setw(15) << kin_en(i);
      ofile << std::setw(15) << pot_en(i) << "\n";
    }
    ofile.close();

  }

};



class solar_system {
private:
  int sun_fix,num_planets;
  double tot_mass,G,tot_momentum,dt,rad,GM_fixed,fixed_force;
  double ax_prev,ay_prev,az_prev,sum_forcex,sum_forcey,sum_forcez;
  arma::Col <double> pos,vel,mass,p_en,k_en;
  arma::Col <double> com = arma::vec(3);
  std::list<std::string> planet_name;
public:
  int N;
  double T;
  solar_system()
  {
    N = 1000;
    T = 10;
    G = GM;                                           // divided by one sun mass
    num_planets=tot_mass=tot_momentum=sun_fix=0;
    sum_forcex=sum_forcey=sum_forcez=GM_fixed=0;
  }


  void add_planet(object planet, std::string name)
  {
    num_planets += 1;
    planet_name.push_back(name);
    tot_mass += planet.mass;
    mass.resize(num_planets);
    mass(num_planets-1) = planet.mass;

    pos.resize(3*num_planets);
    pos(3*(num_planets-1)) = planet.x(0);
    pos(3*(num_planets-1)+1) = planet.y(0);
    pos(3*(num_planets-1)+2) = planet.z(0);

    vel.resize(3*num_planets);
    vel(3*(num_planets-1)) = planet.vx(0);
    vel(3*(num_planets-1)+1) = planet.vy(0);
    vel(3*(num_planets-1)+2) = planet.vz(0);
  }

  // Sun is stationary at origin
  void sun_fixed()
  {
    GM_fixed = GM;                // usually zero
    sun_fix = 1;                  // usually zero
  }

  // includes the Sun as a planet
  // making sure the Sun is placed as the last index
  void sun_included()
  {
    num_planets += 1;
    planet_name.push_back("Sun");
    tot_mass += 1;
    mass.resize(num_planets);
    mass(num_planets-1) = 1;                                  // sun mass relative to the sun
    pos.resize(3*num_planets);
    pos(3*(num_planets-1)) = 0;
    pos(3*(num_planets-1)+1) = 0;
    pos(3*(num_planets-1)+2) = 0;
    vel.resize(3*num_planets);

    for (int i=0; i<num_planets-1; i++){                      // making sure the total momentum of the system is zero
      vel(3*(num_planets-1)) -= mass(i)*vel(3*i);             // divided by one sun mass
      vel(3*(num_planets-1)+1) -= mass(i)*vel(3*i+1);         // divided by one sun mass
      vel(3*(num_planets-1)+2) -= mass(i)*vel(3*i+2);         // divided by one sun mass
    }
  }

  // finds center of mass in all dimensions
  void center_of_mass()
  {
    com(0)=com(1)=com(2) = 0;
    for (int i=0; i<num_planets; i++){
      com(0) += mass(i)*pos(3*i);
      com(1) += mass(i)*pos(3*i+1);
      com(2) += mass(i)*pos(3*i+2);
    }
    com /= (tot_mass+sun_fix);
    std::cout << "Center of mass is (x,y,z) \n" << com << std::endl;
  }

  // finds total momentum in solar system. Used in unit test
  void total_momentum()
  {
    tot_momentum=0;
    for (int i=0; i<num_planets; i++){
      tot_momentum += mass(i)*(vel(3*i)+vel(3*i+1)+vel(3*i+2));
    }
  }

  // calculate distance between two planets based on the coordinate
  double distance(int planet1, int planet2)
  {
    double X,Y,Z,distance;
    std::string p1,p2;

    X = pos(3*planet1)-pos(3*planet2);
    Y = pos(3*planet1+1)-pos(3*planet2+1);
    Z = pos(3*planet1+2)-pos(3*planet2+2);
    distance = sqrt(X*X+Y*Y+Z*Z);
    if (distance < eps){                // Unit test if distance between planets are too small
      std::cout << "Distance between planet " << planet1 << " and " << planet2 << " is under "
                << eps << " AU! The program is stopped!" << std::endl;
      std::exit(0);
    }
    else{
      return distance;
    }
  }

  // Newton's law of motion due to gravity
  double force(int planet1, int planet2)
  {
    return G*mass(planet1)*mass(planet2) / pow(distance(planet1,planet2),3);
  }

  void solve()
  {
    std::string fileout = "data_orbits.txt";
    std::string en_file = "energy.txt";
    ofile.open(fileout);
    enfile.open(en_file);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    enfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    dt = T/N;
    p_en = arma::vec(num_planets);
    k_en = arma::vec(num_planets);

    ofile << std::setw(15);
    for (auto v : planet_name)
          ofile << v << std::setw(45);
    ofile << "\n";


    for (int k=0; k<num_planets; k++){
      ofile << std::setw(15) << pos(3*k);
      ofile << std::setw(15) << pos(3*k+1);
      ofile << std::setw(15) << pos(3*k+2);
    }
    ofile << "\n";

    for (int iter=0; iter<N; iter++){
      for (int i=0; i<num_planets; i++){
        p_en(i) = -sun_fix*GM*mass(i)/sqrt(pos(3*i)*pos(3*i)+pos(3*i+1)*pos(3*i+1)+pos(3*i+2)*pos(3*i+2));
        k_en(i) = 0.5*mass(i)*(vel(3*i)*vel(3*i) + vel(3*i+1)*vel(3*i+1) + vel(3*i+2)*vel(3*i+2));

        for (int j=0; j<num_planets; j++){
          if (i!=j){
            p_en(i) -= G*mass(i)*mass(j)/distance(i,j);

            sum_forcex += force(i,j) * (pos(3*i) - pos(3*j));
            sum_forcey += force(i,j) * (pos(3*i+1) - pos(3*j+1));
            sum_forcez += force(i,j) * (pos(3*i+2) - pos(3*j+2));
          }
        }

          rad = sqrt(pow(pos(3*i),2) + pow(pos(3*i+1),2) + pow(pos(3*i+2),2));
          if (fabs(rad) > eps){fixed_force = GM_fixed/pow(rad,3);}

          ax_prev = -sum_forcex/mass(i) - fixed_force*pos(3*i);
          ay_prev = -sum_forcey/mass(i) - fixed_force*pos(3*i+1);
          az_prev = -sum_forcez/mass(i) - fixed_force*pos(3*i+2);

          pos(3*i) = pos(3*i) + dt*vel(3*i) + 0.5*dt*dt*ax_prev;
          pos(3*i+1) = pos(3*i+1) + dt*vel(3*i+1) + 0.5*dt*dt*ay_prev;
          pos(3*i+2) = pos(3*i+2) + dt*vel(3*i+2) + 0.5*dt*dt*az_prev;

          sum_forcex=sum_forcey=sum_forcez=0;
          for (int j=0; j<num_planets; j++){
            if (i!=j){
              sum_forcex += force(i,j) * (pos(3*i) - pos(3*j));
              sum_forcey += force(i,j) * (pos(3*i+1) - pos(3*j+1));
              sum_forcez += force(i,j) * (pos(3*i+2) - pos(3*j+2));
            }
          }

          rad = sqrt(pow(pos(3*i),2) + pow(pos(3*i+1),2) + pow(pos(3*i+2),2));
          if (fabs(rad) > eps){fixed_force = GM_fixed/pow(rad,3);}

          vel(3*i) = vel(3*i) + 0.5*dt*(-sum_forcex/mass(i) - fixed_force*pos(3*i) + ax_prev);
          vel(3*i+1) = vel(3*i+1) + 0.5*dt*(-sum_forcey/mass(i) - fixed_force*pos(3*i+1) + ay_prev);
          vel(3*i+2) = vel(3*i+2) + 0.5*dt*(-sum_forcez/mass(i) - fixed_force*pos(3*i+2) + az_prev);

          sum_forcex=sum_forcey=sum_forcez=0;

      }

      if (sun_fix==0 && iter%100==0){
        total_momentum();
        if (fabs(tot_momentum)>eps){
          std::cout << "Total momentum is too large! Program stopped! \n" <<
                       "Total momentum: " << tot_momentum << "       Tolerance(+/-): " << eps << std::endl;
          std::exit(0);
        }
      }

      for (int k=0; k<num_planets; k++){
        ofile << std::setw(15) << pos(3*k);
        ofile << std::setw(15) << pos(3*k+1);
        ofile << std::setw(15) << pos(3*k+2);
        enfile << std::setw(15) << p_en(k);
        enfile << std::setw(15) << k_en(k);
      }
      ofile << "\n";
      enfile << "\n";
    }
    ofile.close();


  }

};



/*
void planet(arma::Col <double> &x, arma::Col <double> &y, arma::Col <double> &z,
            arma::Col <double> &vx, arma::Col <double> &vy,
            arma::Col <double> &vz, std::string method, double mass){
  double a,r,ax,ay,az;

  int n = x.n_elem;

  arma::Col <double> t = arma::vec(n);
  double tmin = 0;
  double tmax = 5;
  double dt = (tmax-tmin)/n;
  for (int i=0; i<n; i++){
    t(i)=i*dt;
  }

  double red_mass = mass/M_sun;
  //set up vectors containing potential and kinetic energy
  arma::Col <double> kin_en = arma::vec(n); kin_en(0) = 0.5*red_mass*(vx(0)*vx(0)+vy(0)*vy(0)+vz(0)*vz(0));
  arma::Col <double> pot_en = arma::vec(n); pot_en(0) = -GM*red_mass/sqrt(x(0)*x(0)+y(0)*y(0)+z(0)*z(0));

  if (method == "euler"){
    //Solving ode using forward Euler
    for (int i=1; i<n; i++){
      r = sqrt(x(i-1)*x(i-1) + y(i-1)*y(i-1) + z(i-1)*z(i-1));
      a = GM / (r*r*r);
      ax = -a*x(i-1);
      ay = -a*y(i-1);
      az = -a*z(i-1);
      vx(i) = vx(i-1) + ax*dt;
      vy(i) = vy(i-1) + ay*dt;
      vz(i) = vz(i-1) + az*dt;
      x(i) = x(i-1) + vx(i-1)*dt;
      y(i) = y(i-1) + vy(i-1)*dt;
      z(i) = z(i-1) + vz(i-1)*dt;

      kin_en(i) = 0.5*red_mass*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i));
      pot_en(i) = -GM*red_mass/sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i));
    }
    //write results to file for analysis
    std::string fileout = "Orbit_euler.txt";
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    for (int i=0; i<n; i++){
      ofile << std::setw(15) << x(i);
      ofile << std::setw(15) << y(i);
      ofile << std::setw(15) << z(i);
      ofile << std::setw(15) << t(i);
      ofile << std::setw(15) << kin_en(i);
      ofile << std::setw(15) << pot_en(i) << "\n";
    }
  }

  else if (method == "verlet"){
    //Solving ode using velocity Verlet
    double ax_new,ay_new,az_new,ax_prev,ay_prev,az_prev;
    double vx_half,vy_half,vz_half;
    double r0 = sqrt(x(0)*x(0) + y(0)*y(0) + z(0)*z(0));
    double a0 = GM / (r0*r0*r0);
    ax_prev = -a0*x(0);
    ay_prev = -a0*y(0);
    az_prev = -a0*z(0);
    for (int i=1; i<n; i++){
      x(i) = x(i-1) + dt*vx(i-1) + 0.5*dt*dt*ax_prev;
      y(i) = y(i-1) + dt*vy(i-1) + 0.5*dt*dt*ay_prev;
      z(i) = z(i-1) + dt*vz(i-1) + 0.5*dt*dt*az_prev;

      r = sqrt(x(i)*x(i) + y(i)*y(i) + z(i)*z(i));
      a = GM / (r*r*r);
      ax_new = -a * x(i);
      ay_new = -a * y(i);
      az_new = -a * z(i);

      vx(i) = vx(i-1) + 0.5*dt*(ax_new + ax_prev);
      vy(i) = vy(i-1) + 0.5*dt*(ay_new + ay_prev);
      vz(i) = vz(i-1) + 0.5*dt*(az_new + az_prev);

      ax_prev = ax_new;
      ay_prev = ay_new;
      az_prev = az_new;


      kin_en(i) = 0.5*red_mass*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i));
      pot_en(i) = -GM*red_mass/sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i));
    }
    //write results to file for analysis
    std::string fileout = "Orbit_verlet.txt";
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    for (int i=0; i<n; i++){
      ofile << std::setw(15) << x(i);
      ofile << std::setw(15) << y(i);
      ofile << std::setw(15) << z(i);
      ofile << std::setw(15) << t(i);
      ofile << std::setw(15) << kin_en(i);
      ofile << std::setw(15) << pot_en(i) << "\n";
    }
  }

  else{
    std::cout << "You must choose either (euler) or (verlet)\n";
    exit(1);
    }
  ofile.close();
} // end of function planet()
*/



int main(int argc, char* argv[]){

  ofile << "\n";
  std::string method = argv[1];
  int len = atoi(argv[2]);        //number of integration points
  double int_vel = atof(argv[3]); //initial velocity
  double time = atof(argv[4]);
  double esc_vel = sqrt(2*GM);
  double jup_scale;

/*
  std::cout << esc_vel << "\n";
  arma::Col <double> x = arma::vec(len); x(0)=1;
  arma::Col <double> y = arma::vec(len); y(0)=0;
  arma::Col <double> z = arma::vec(len); z(0)=0;
  arma::Col <double> vx = arma::vec(len); vx(0)=0;
  arma::Col <double> vy = arma::vec(len); vy(0)=2*pi;
  arma::Col <double> vz = arma::vec(len); vz(0)=0;

  planet(x,y,z,vx,vy,vz,method,M_earth);
*/


  // object( x0, y0, z0, vx0, vy0, vz0, M, points )
  // distance is given in AU, velocity in AU/day and mass in kg
  // Positions and velocity data gathered from https://ssd.jpl.nasa.gov/horizons.cgi
  // at A.D. 2019-Dec-09 00:00:00.0000 TDB
  object mercury(-3.985847784965280E-01,-8.678484206044013E-02,2.854818397850725E-02,6.808061022618487E-04,-2.615697349455290E-02, -2.200251411818960E-03,M_mercury,len);
  object venus(6.140254422268607E-01,-3.889781916531376E-01,-4.077096546312043E-02,1.070185289465133E-02,1.700477808956028E-02,-3.842439888550384E-04,M_venus,len);
  object mars(-1.485032517264654E+00, -6.306157101254950E-01, 2.322202328310920E-02,5.992165013982506E-03,-1.168365481307998E-02,-3.918498445436787E-04,M_mars,len);
  object saturn(3.685089251558515E+00,-9.335591564910553E+00,1.562158057974095E-02,4.889009775915366E-03,2.032733431539527E-03,-2.295408335647753E-04,M_saturn,len);
  jup_scale = 1;     // increase the mass of Jupiter
  object jupiter(3.551315858851771E-01,-5.223858708443553E+00,1.375193093344411E-02,7.445397359016055E-03,8.688615308896841E-04,-1.701937692576648E-04,M_jupiter*jup_scale,len);
  object earth(2.328416719695888E-01, 9.570420225654582E-01,-4.193306777199945E-05,-1.699305780122259E-02,3.997104358502586E-03,-4.831893976607005E-07,M_earth,len);
  object uranus(1.627777749498813E+01,1.130905239963674E+01,-1.688216806579894E-01,-2.265866949228651E-03,3.047569009304266E-03,4.052178469796985E-05,M_uranus,len);
  object neptune(2.922766815589142E+01,-6.438194386201971E+00,-5.410875794296358E-01,6.618180582706258E-04,3.085812272712285E-03,-7.886168713184974E-05,M_neptune,len);


  // The earth-sun system where the sun is fixed at origin
  // using the object class
  /*
  earth.velocity_verlet();
  earth.kinetic_energy();
  earth.potential_energy();
  earth.write_to_file("Orbit_verlet.txt");
  */

  // Setting up the Solar System. Used in following tasks
  solar_system system;

  // The Earth-Jupiter-Sun system with Sun fixed at origin
  system.N = len;
  system.T = time;


  // Adding all planets in the Solar system that's
  // been initialized through the class "object".
  // Choose the Sun to be stationary or in motion (sun_included() or sun_fixed())
  /*
  system.N = len;
  system.T = time;
  system.add_planet(mercury,"Mercury");
  system.add_planet(venus,"Venus");
  system.add_planet(earth, "Earth");
  system.add_planet(mars,"Mars");
  system.add_planet(jupiter, "Jupiter");
  system.add_planet(saturn,"Saturn");
  system.add_planet(uranus,"Uranus");
  system.add_planet(neptune,"Neptune");
  system.sun_included();
  //system.sun_fixed();
  system.solve();
  */

  // General theory of relativity correction on Mercury
  /*
  mercury.tmax = time;
  mercury.mercury_perihelion_initialize();
  mercury.velocity_verlet();
  mercury.kinetic_energy();
  mercury.potential_energy();
  mercury.write_to_file("mercury_perihelion.txt");
  */

  return 0;
} // end of function main()
