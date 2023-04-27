#include<iostream>
#include<cmath>
#include<chrono>
#include<iomanip>

using namespace std;
using namespace std::chrono;
auto start =high_resolution_clock::now();

double f(double);

void Euler_Explicit();
void Euler_Implicit();
void Leap_Frog();
void Crank_Nicholson();
void Second_Adams_Bashworth();
void Second_Runge_Kutta();
void Fourth_Runge_Kutta();

double t;
double dt;
double theta_0;
double Error;
double theta_an;

int main(){
    t=0;
    theta_0 = 100;
    theta_an = 100*exp(-2.5);       
    //cout<<" Enter the measure of time step :\t"<<setprecision(15);
   // cin>>dt;
    Error=0;
    dt= 1.95626e-06;
  
  auto start =high_resolution_clock::now();
      Euler_Implicit();
    cout<<" dt "<<dt << " Error "<<Error<<endl;
    auto stop = high_resolution_clock::now();
    auto duration =duration_cast<microseconds>(stop -start);
    cout<<'\n'<<"The time taken for convergence is : "<<duration.count() <<"\n";
  
  
  
  
    //Euler_Explicit();
    //Euler_Implicit();
    //Leap_Frog();
    // Crank_Nicholson();
    //  Second_Adams_Bashworth();
    // Second_Runge_Kutta();
    // Fourth_Runge_Kutta();

}

void Euler_Explicit(){
    double theta_np1 =0; double theta_n= theta_0 ;
    while( t< 5.0 ){
        theta_np1 = theta_n - theta_n*(dt)*0.5;
        theta_n = theta_np1;
        t = t + dt;
    }
    Error = abs( theta_np1 - 100*exp(-t/2.0) ) ;
    //cout<<" EULER EXPLICIT "<<"\t : " <<theta_np1<<" ERROR : "<<Error<<endl;

}


void Euler_Implicit(){
    t=0; 
    double theta_np1=0; double theta_n = theta_0;
    while( t < 5.0 ){
        theta_np1 = 2*theta_n/(2 + dt);
        theta_n = theta_np1;
        t = t+dt;   
    }
    Error = abs( theta_np1 - 100*exp(-t/2.0) ) ;
    cout<<" EULER IMPLICIT "<<"\t : " <<theta_np1<<" ERROR : "<<Error<<endl;

}

void Leap_Frog(){
    t=0; 
    double theta_np1=0; double theta_np_half ; double theta_n = theta_0;
    theta_np_half = 100*exp(-dt*0.25);
    while( t < 5.0 ){
          theta_np1 = theta_n + dt*f(theta_np_half);
          theta_np_half = theta_np_half + dt*f(theta_np1);
          theta_n = theta_np1;
          t= t+dt;
    }
    Error = abs( theta_np1 - 100*exp(-t/2.0) ) ;
    cout<<" LEAP FROG "<<"\t\t : " <<theta_np1<<" ERROR : "<<Error<<endl;

}

void Crank_Nicholson(){
    t=0; 
    double theta_np1=0; double theta_n = theta_0;
    while( t < 5.0 ){
        theta_np1 = 2*theta_n/(2 + dt);
        theta_np1 = theta_n + dt*0.5*( f(theta_n) + f(theta_np1) );        
        theta_n = theta_np1;
        t = t + dt;
    }
    Error = abs( theta_np1 - 100*exp(-t/2.0) ) ;
    cout<<" CRANK NICHOLSON "<<"\t : " <<theta_np1<<" ERROR : "<<Error<<endl;
}

void Second_Adams_Bashworth(){
    t=0;
    double theta_n_minus_1 = theta_0;  double theta_n = 100*exp(-dt*0.5); double theta_np1 = 0;
    t = t+dt;
    while( t < 5.0 ){
        theta_np1 = theta_n + dt*( 1.5*f(theta_n) - 0.5*f(theta_n_minus_1) ) ;
        theta_n_minus_1 = theta_n;
        theta_n = theta_np1;
        t = t +dt;
    } 
    Error = abs( theta_np1 - 100*exp(-t/2.0) ) ;
    cout<<" 2nd Or Adams Bashworth "<<" : " <<theta_np1<<" ERROR : "<<Error<<endl;
}

void Second_Runge_Kutta(){
    t=0;
    double theta_n = theta_0; double theta_np1 = 0; double theta_np_half_star=0;
    while( t < 5.0 ){
        theta_np_half_star = theta_n + dt/2*f(theta_n);
        theta_np1 = theta_n + dt*f(theta_np_half_star);
        theta_n = theta_np1;
        t = t +dt;
    }
    Error = abs( theta_np1 - 100*exp(-t/2.0) ) ;
    cout<<" 2nd Order Runge Kutta "<<"\t : " <<theta_np1<<" ERROR : "<<Error<<endl;

}

void Fourth_Runge_Kutta(){
    t=0;
    double theta_n = theta_0; double theta_np1 = 0; double theta_np_half_star = 0; double theta_np_half_star_star = 0; double theta_np1_star=0;
    while( t < 5.0 ){
        theta_np_half_star = theta_n + dt*0.5*f(theta_n);
        theta_np_half_star_star = theta_n + dt*0.5*f(theta_np_half_star);
        theta_np1_star = theta_n + dt*f(theta_np_half_star_star);
        theta_np1 = theta_n + (dt/6.0)*( f(theta_n) + 2*f(theta_np_half_star) + 2*f(theta_np_half_star_star) + f(theta_np1_star) );
        // theta_np_half_star = theta_n*( 1 - dt*0.25 ); //changed 0.25 to 0.5
        // theta_np_half_star_star = theta_n - dt*0.25*theta_np_half_star;
        // theta_np1_star = theta_n - dt*0.5*theta_np_half_star_star;
        // theta_np1 = theta_n - (theta_n + 2*theta_np_half_star + 2*theta_np_half_star_star + theta_np1_star )*dt/12 ; 
        theta_n = theta_np1;
        t = t +dt;
        //cout<<"TIME "<<t<<endl;
    }
    
    Error = abs( theta_np1 - 100*exp(-t/2.0) ) ;
    cout<<" 4th order Runge Kutta "<<"\t : " <<theta_np1<<" ERROR : "<<Error<<endl;
    cout<<"TIME  "<<t<<endl;
}


double f( double theta ){
    return -theta/2;
}