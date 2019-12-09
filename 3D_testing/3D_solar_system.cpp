#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include <sstream>
#include <string>
#include <ctime>
//#include "mpi.h"
//#include "orbit.h"

//For debugging:
// compile with: g++ solar_system.cpp -o main1 -larmadillo -llapack -lblas
// execute with ./main1

//For paralellization
//mpicxx  -o main_mpi.x  main.cpp -std=c++11
//mpiexec -n 2 ./main_mpi.x


std::ofstream ofile;
double const pi = 3.14159265359;
double const M_sun = 2*pow(10,30);        // in kilograms
double const M_earth = 6*pow(10,24);      // in kilograms
double const sun_rad = 0.00465047;        // in AU
double const GM = 4*pi*pi;

void planet(arma::Col <double> &x, arma::Col <double> &y, arma::Col <double> &z,
            arma::Col <double> &vx, arma::Col <double> &vy,
            arma::Col <double> &vz, std::string method, double mass, double beta,double tmax){
  double a,r,ax,ay,az;

  std::clock_t start;
  double duration;

  int n = x.n_elem;

  arma::Col <double> t = arma::vec(n);
  double tmin = 0;
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
    start = std::clock();
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


    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Euler time with "<< n << " iterations:" << duration <<'\n';
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
    double a0 = GM / pow(r0,beta+1);
    ax_prev = -a0*x(0);
    ay_prev = -a0*y(0);
    az_prev = -a0*z(0);
    start = std::clock();
    for (int i=1; i<n; i++){
      x(i) = x(i-1) + dt*vx(i-1) + 0.5*dt*dt*ax_prev;
      y(i) = y(i-1) + dt*vy(i-1) + 0.5*dt*dt*ay_prev;
      z(i) = z(i-1) + dt*vz(i-1) + 0.5*dt*dt*az_prev;

      r = sqrt(x(i)*x(i) + y(i)*y(i) + z(i)*z(i));
      a = GM / pow(r,beta+1);
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
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Verlet time with "<< n << " iterations:" << duration <<'\n';
    //write results to file for analysis
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
} // end of function planet()




int main(int argc, char* argv[]){


  std::string method = argv[1];
  int len = atoi(argv[2]);//number of integration points
  double int_vel = atof(argv[3]);//initial velocity
  double beta = atof(argv[4]);// exponent of the radius in the gravitational force
  double orbital_time = atof(argv[5]);//simulation time in years
  double esc_vel = sqrt(2*GM);
  std::cout << esc_vel << "\n";
  arma::Col <double> x = arma::vec(len); x(0)=1.0;
  arma::Col <double> y = arma::vec(len); y(0)= 0.0;
  arma::Col <double> z = arma::vec(len); z(0)= 0.0;
  arma::Col <double> vx = arma::vec(len); vx(0)=0.0;
  arma::Col <double> vy = arma::vec(len); vy(0)= int_vel;
  arma::Col <double> vz = arma::vec(len); vz(0)=0.0;

  std::string fileout = "Orbit_diff_r_exp.txt";
  ofile.open(fileout);
  for (double b = 2; b <= beta+0.1; b += 0.2){
      std::cout << b << std::endl;
      planet(x,y,z,vx,vy,vz,method,M_earth,b,orbital_time);
  }
  ofile.close();


  return 0;
} // end of function main()
