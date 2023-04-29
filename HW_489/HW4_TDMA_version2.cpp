#include<iostream>
#include<cmath>
using namespace std;
//GAUSSIAN ELIMINATION :
void TDMA(double [],double [],double [],double [],double [], int);
void printmatrix(double [], int );
double k=25;double T_s=365.15;

int main(){
//INPUTS	
      int n;
      cout<<"Enter number of nodes :"<<endl;
      cin>>n;

//PROCESSED PARAMETERS
        double del_x=0.1/(n-1);     double q=(-(0.3E+06)/k)*pow(del_x,2);     n=n-1; //Enforcing BC at right end and we need to compute only for n-1 nodes
        double W[n]={}; double D[n]={}; double  E[n]={}; double   Q[n]={}; double sol[n]={}; //D=DIAGONAL COEFF , E=EAST COEFF , W=WEST COEFF , Q=LOAD VECTOR , sol=SOLUTION VEC 

//WRITING THE COEFFICIENT MATRIX
        D[0]=-1;     E[0]=1;     Q[0]=0;
        W[n-1]=1;   D[n-1]=-2;  Q[n-1]=q-T_s;
        for(int i=1;i<n-1;++i){
          E[i]=W[i]=1; D[i]=-2; 
            Q[i]=q; 
        }

TDMA(W,D,E,Q,sol,n);
printmatrix(sol,n);
}


void TDMA(double W[],double D[], double E[] , double Q[] , double sol[], int n){  
//FORWARD SUBSTITUTION to reduce to UPPER TRIANGULAR MATRIX
    for(int i=1;i<n;++i){
      D[i]=  D[i]  -  (W[i]/D[i-1])*E[i-1];
      Q[i]=Q[i]    -  (W[i]/D[i-1])*Q[i-1];     
    }

//BACKWRARD SUBSTITUTION for SOLUTION
    sol[n-1]= Q[n-1]/D[n-1];
    for(int i=n-2;i>-1;--i){
      sol[i]= (Q[i]- E[i]*sol[i+1])/D[i];   
    } 

}

void printmatrix(double dummy[], int n){
  cout<<'\n'<<"TEMPERATURE AT THE NODES ARE (degrees Centigrade):";
  for(int i=0;i<n;++i){
    cout<<dummy[i]-273.15<<',';
  } 
  cout<<T_s-273.15<<endl<<endl;
}
