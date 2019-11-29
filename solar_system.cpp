#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include <sstream>
#include <string>
//#include "mpi.h"

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


void planet(arma::Col <double> &x, arma::Col <double> &y, arma::Col <double> &vx,
            arma::Col <double> &vy, std::string method){
  double a,r,ax,ay;
  double GM = 4*pi*pi;

  int n = x.n_elem;

  arma::Col <double> t = arma::vec(n);
  double tmin = 0;
  double tmax = 5;
  double dt = (tmax-tmin)/n;
  for (int i=0; i<n; i++){
    t(i)=i*dt;
  }

  double red_earth = M_earth/M_sun;
  //set up vectors containing potential and kinetic energy
  arma::Col <double> kin_en = arma::vec(n); kin_en(0) = 0.5*red_earth*(vx(0)*vx(0)+vy(0)*vy(0));
  arma::Col <double> pot_en = arma::vec(n); pot_en(0) = -GM*red_earth/sqrt(x(0)*x(0)+y(0)*y(0));

  if (method == "euler"){
    //Solving ode using forward Euler
    for (int i=1; i<n; i++){
      r = sqrt(x(i-1)*x(i-1) + y(i-1)*y(i-1));
      a = GM / (r*r);
      ax = -a*x(i-1);
      ay = -a*y(i-1);
      vx(i) = vx(i-1) + ax*dt;
      vy(i) = vy(i-1) + ay*dt;
      x(i) = x(i-1) + vx(i-1)*dt;
      y(i) = y(i-1) + vy(i-1)*dt;

      kin_en(i) = 0.5*red_earth*(vx(i)*vx(i)+vy(i)*vy(i));
      pot_en(i) = -GM*red_earth/sqrt(x(i)*x(i)+y(i)*y(i));
    }
    //write results to file for analysis
    std::string fileout = "Orbit_euler.txt";
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    for (int i=0; i<n; i++){
      ofile << std::setw(15) << x(i);
      ofile << std::setw(15) << y(i);
      ofile << std::setw(15) << t(i);
      ofile << std::setw(15) << kin_en(i);
      ofile << std::setw(15) << pot_en(i) << "\n";
    }
  }

  else if (method == "verlet"){
    //Solving ode using velocity Verlet
    double ax_new, ay_new;
    double ax_prev,ay_prev;
    double vx_half,vy_half;
    double r0 = sqrt(x(0)*x(0) + y(0)*y(0));
    double a0 = GM / (r0*r0);
    ax_prev = -a0*x(0);
    ay_prev = -a0*y(0);
    for (int i=0; i<n-1; i++){
      vx_half = vx(i) + 0.5*ax_prev*dt;
      vy_half = vy(i) + 0.5*ay_prev*dt;

      x(i+1) = x(i) + vx_half*dt;
      y(i+1) = y(i) + vy_half*dt;

      r = sqrt(x(i+1)*x(i+1) + y(i+1)*y(i+1));
      ax_new = -GM*x(i)/pow(r,3);
      ay_new = -GM*y(i)/pow(r,3);

      vx(i+1) = vx_half + 0.5*ax_new*dt;
      vy(i+1) = vy_half + 0.5*ay_new*dt;

      ax_prev = ax_new;
      ay_prev = ay_new;

      /*
      ax_new = -a*x(i-1);
      ay_new = -a*y(i-1);

      x(i) = x(i-1) + dt*vx(i-1) + dt*dt*ax_prev/2;
      y(i) = y(i-1) + dt*vy(i-1) + dt*dt*ay_prev/2;

      vx(i) = vx(i-1)+dt/2*(ax_new + ax_prev);
      vy(i) = vy(i-1)+dt/2*(ay_new + ay_prev);

      ax_prev = ax_new;
      ay_prev = ay_new;
      */

      kin_en(i) = 0.5*red_earth*(vx(i)*vx(i)+vy(i)*vy(i));
      pot_en(i) = -GM*red_earth/sqrt(x(i)*x(i)+y(i)*y(i));
    }
    //write results to file for analysis
    std::string fileout = "Orbit_verlet.txt";
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    for (int i=0; i<n; i++){
      ofile << std::setw(15) << x(i);
      ofile << std::setw(15) << y(i);
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




int main(int argc, char* argv[]){
  int len = 100000;
  std::string method = argv[1];
  arma::Col <double> x = arma::vec(len); x(0)=1;
  arma::Col <double> y = arma::vec(len); y(0)=0;
  arma::Col <double> vx = arma::vec(len); vx(0)=0;
  arma::Col <double> vy = arma::vec(len); vy(0)=5;//sqrt(4*pi*pi);


  planet(x,y,vx,vy,method);




  return 0;
} // end of function main()
