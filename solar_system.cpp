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
// execute with ./main1

//For paralellization
//mpicxx  -o main_mpi.x  main.cpp -std=c++11
//mpiexec -n 2 ./main_mpi.x

std::ofstream ofile;
double const pi = 3.14159265359;
double const GM = 4*pi*pi;                  // gravitational constant times one solar mass
double const sun_rad = 0.00465047;          // in AU
double const M_sun = 2*pow(10,30);          // in kilograms
double const M_earth = 6*pow(10,24);        // in kilograms
double const M_jupiter = 1.9*pow(10,27);    // in kilograms



class object {
private:
  double a, r, ax, ay, az, dt;
  double ax_new, ay_new, az_new, ax_prev, ay_prev, az_prev;
  double vx_half, vy_half, vz_half;
public:
  //double x0,y0,z0,vx0,vy0,vz0,mass;
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
    vx(0) = vx0;
    vy(0) = vy0;
    vz(0) = vz0;
    dt = tmax/N;

  }

  // functions
  void velocity_verlet()
  {
    r = sqrt(x(0)*x(0) + y(0)*y(0) + z(0)*z(0));
    a = GM / (r*r*r);
    ax_prev = -a * x(0);
    ay_prev = -a * y(0);
    az_prev = -a * z(0);
    for (int i=1; i<N; i++){
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
    }
  }

  void euler()
  {
    for (int i=1; i<N; i++){
      r = sqrt(x(i-1)*x(i-1) + y(i-1)*y(i-1) + z(i-1)*z(i-1));
      a = GM / (r*r);
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
    ofile.open(filename);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);

    for (int i=0; i<N; i++){
      ofile << std::setw(15) << x(i);
      ofile << std::setw(15) << y(i);
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
  int num_planets,N,nothing;
  double tot_mass,G,tot_momentum,com,r_com,dt;
  double ax_prev,ay_prev,az_prev;
  arma::Col <double> pos,vel,mass;
  std::string = word;
  std::list<std::string> planet_name;
public:
  double T;
  solar_system(int len)
  {
    N = len;
    T = 10;
    G = GM;                                           // divided by one sun mass
    num_planets=tot_mass=tot_momentum=0;
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

  void sun_fixed()
  {

  }


  void sun_included()
  {
    num_planets += 1;
    tot_mass += 1;
    mass.resize(num_planets);
    mass(num_planets-1) = 1;                          // sun mass relative to the sun
    pos.resize(3*num_planets);
    pos(3*(num_planets-1)) = 0;
    pos(3*(num_planets-1)+1) = 0;
    pos(3*(num_planets-1)+2) = 0;
    vel.resize(3*num_planets);

    for (int i=0; i<num_planets-1; i++){
      vel(3*(num_planets-1)) -= mass(i)*vel(3*i);             // divided by one sun mass
      vel(3*(num_planets-1)+1) -= mass(i)*vel(3*i+1);         // divided by one sun mass
      vel(3*(num_planets-1)+2) -= mass(i)*vel(3*i+2);         // divided by one sun mass
    }

    for (int i=0; i<num_planets; i++){
      tot_momentum += mass(i)*vel(3*i);
      tot_momentum += mass(i)*vel(3*i+1);
      tot_momentum += mass(i)*vel(3*i+2);
    }
    //std::cout << tot_momentum << std::endl << "\n";
  }

  double center_of_mass()
  {
    for (int i=0; i<num_planets; i++){
      r_com = sqrt(pos(i)*pos(i)+pos(i+1)*pos(i+1)+pos(i+2)*pos(i+2));
      com += r_com*mass(i);
    }
    com /= tot_mass;
    return com;
  }

  double distance(int planet1, int planet2)
  {
    double x,y,z;
    x = pos(3*planet1)-pos(3*planet2);
    y = pos(3*planet1+1)-pos(3*planet2+1);
    z = pos(3*planet1+2)-pos(3*planet2+2);
    return sqrt(x*x+y*y+z*z);
  }

  double mass_prod(int planet1, int planet2)
  {
    return mass(planet1)*mass(planet2);
  }

  double force(int planet1, int planet2)
  {
    return G*mass_prod(planet1,planet2) / pow(distance(planet1,planet2),3);
  }

  void solve()
  {
    std::string fileout = "data_orbits.txt";
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    dt = T/N;

    //for (int k=0; k<num_planets; k++){
    //  ofile << std::setw(15) << planet_name(k);
    //  ofile << std::setw(15) << planet_name(k);
    //  ofile << std::setw(15) << planet_name(k);
    //}
    ofile << "\n";
    for (int k=0; k<num_planets; k++){
      ofile << std::setw(15) << pos(3*k);
      ofile << std::setw(15) << pos(3*k+1);
      ofile << std::setw(15) << pos(3*k+2);
    }
    ofile << "\n";

    for (int iter=0; iter<N; iter++){

      for (int i=0; i<num_planets; i++){
        for (int j=0; j<num_planets; j++){
          if (i==j) {nothing=0;}
          else{
            ax_prev = -force(i,j)*pos(3*i)/mass(i);
            ay_prev = -force(i,j)*pos(3*i+1)/mass(i);
            az_prev = -force(i,j)*pos(3*i+2)/mass(i);

            pos(3*i) = pos(3*i) + dt*vel(3*i) + 0.5*dt*dt*ax_prev;
            pos(3*i+1) = pos(3*i+1) + dt*vel(3*i+1) + 0.5*dt*dt*ay_prev;
            pos(3*i+2) = pos(3*i+2) + dt*vel(3*i+2) + 0.5*dt*dt*az_prev;

            vel(3*i) = vel(3*i) + 0.5*dt*(-force(i,j)*pos(3*i)/mass(i) + ax_prev);
            vel(3*i+1) = vel(3*i+1) + 0.5*dt*(-force(i,j)*pos(3*i+1)/mass(i) + ay_prev);
            vel(3*i+2) = vel(3*i+2) + 0.5*dt*(-force(i,j)*pos(3*i+2)/mass(i) + az_prev);
          }
        }
      }

      for (int k=0; k<num_planets; k++){
        ofile << std::setw(15) << pos(3*k);
        ofile << std::setw(15) << pos(3*k+1);
        ofile << std::setw(15) << pos(3*k+2);
      }
      ofile << "\n";
    }


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
  double tmax = 20;
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
      a = GM / (r*r);
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

  std::string method = argv[1];
  int len = atoi(argv[2]);//number of integration points
  double int_vel = atof(argv[3]);//initial velocity
  double time = atof(argv[4]);
  double esc_vel = sqrt(2*GM);

  //std::cout << esc_vel << "\n";
  //arma::Col <double> x = arma::vec(len); x(0)=1;
  //arma::Col <double> y = arma::vec(len); y(0)=0;
  //arma::Col <double> z = arma::vec(len); z(0)=0;
  //arma::Col <double> vx = arma::vec(len); vx(0)=0;
  //arma::Col <double> vy = arma::vec(len); vy(0)=int_vel;
  //arma::Col <double> vz = arma::vec(len); vz(0)=0;

  //planet(x,y,z,vx,vy,vz,method,M_earth);


  // object( x0, y0, z0, vx0, vy0, vz0, M )
  // distance is given in AU and mass is given in kg
  object earth(1,0,0,0,int_vel,0,M_earth,len);
  object jupiter(0,4,0,2,0,0,M_jupiter,len);
  solar_system system(len);
  system.T = time;
  system.add_planet(earth, "Earth");
  system.add_planet(jupiter, "Jupiter");
  system.sun_included();
  system.solve();

  //earth.velocity_verlet();
  //earth.kinetic_energy();
  //earth.write_to_file("earth.txt");



  return 0;
} // end of function main()
