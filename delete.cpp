#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <utility>
#include <algorithm>


using namespace std;

/*  PARAMETERS  */
int nx =128 , ny=nx;
double lx = 1 , ly=lx ;
double dx = lx/nx , dy = ly/ny;

/* PHYSICAL PARAMETERS */
double Re = 1000 ;
double dt = 1.25E-4  ;

int main(){
    cout<<"CFL CONITION "<<dx<<endl;
    cout<<"NEUMANN CONDITION "<<0.25*pow( dx,2 )*Re<<endl;
    cout<<" NON LINEAR CFD "<<pow(2,2.0/3.0)*pow(dx, 4.0/3.0 )<<endl;
    double t3 = pow(2,2.0/3.0)*pow(dx, 4.0/3.0 );
    cout<<min( {dx , 0.25*pow( dx,2 )*Re , t3}  )<<endl;
}