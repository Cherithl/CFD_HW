#include<iostream>
#include<cmath>
using namespace std;

void GS();
void TDMA();

int n=161;
double dx=1.0/(n-1);
double Pe_L = 50;
double Pe = Pe_L*dx;

int main(){
//GS();
TDMA();
}


void GS(){
    double T[n] = {}; T[0] =0 ; T[n-1] = 1;
    double res[n] = {}; int iter=0;
    double res_avg =1;double tol=-8;


while( log10(res_avg) > tol ){
    res_avg=0;
    for( int i=1 ; i < n-1 ; ++i){
        T[i] =   1.0/(2+Pe)*( (1+Pe)*T[i-1] + T[i+1]  );
    } 
//RESIDUAL
    for( int i=1 ; i<n-1 ; ++i){
        res[i]= -(1+Pe)*T[i-1]  + (2+Pe)*T[i] - T[i+1]  ;
        res_avg += pow(res[i],2)/(n-2); 
    }
    res_avg += pow(res_avg,0.5); 
    ++iter;
}

    for( int i=0 ; i<n ;++i){
        cout<<T[i]<<endl;
    }

    cout<<iter<<endl;

}



void TDMA(){  

double Q[n-2] = {};
double W[n-2] = {};
double D[n-2] = {};
double E[n-2] = {};
double sol[n-2] = {};

E[0] = -1; D[0] = 2+Pe;    W[n-3] = -1-Pe; D[n-3] = -2;
Q[n-3] =  pow( 2+Pe ,-1 ) ; 

for(int i=1 ; i<n-2 ; ++i){
  W[i]= -1-Pe  ; D[i]=2+Pe ; E[i] = -1;
}
//FORWARD SUBSTITUTION to reduce to UPPER TRIANGULAR MATRIX
    for(int i=1;i<n-2;++i){
      D[i]=  D[i]  -  (W[i]/D[i-1])*E[i-1];
      Q[i]=Q[i]    -  (W[i]/D[i-1])*Q[i-1];     
    }

//BACKWRARD SUBSTITUTION for SOLUTION
    sol[n-3]= Q[n-3]/D[n-3];
    for(int i=n-4;i>-1;--i){
      sol[i]= (Q[i]- E[i]*sol[i+1])/D[i];   
    } 
cout<<0<<endl;
for(int i=0 ; i<n-2 ; ++i){
    cout<<sol[i]<<endl;
}
cout<<1<<endl;


}