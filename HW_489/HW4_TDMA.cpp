#include<iostream>
#include<cmath>
using namespace std;

int main(){
//INPUTS	
          int n;
          cout<<"Enter number of nodes :"<<endl;
          cin>>n;
//PROCESSED PARAMETERS
          double del_x=0.1/(n-1);   double k=25;  double q=(-(0.3E+06)/k)*pow(del_x,2);    double T_s=365.15;
          n=n-1;
          double A[n][3]={};  double Q[n]={};  double sol[n]={}; 
//WRITING THE COEFFICIENT MATRIX
 
         A[0][1]=-1;     A[0][2]=1;     Q[0]=0;
         A[n-1][0]=1;   A[n-1][1]=-2;   Q[n-1]=q-T_s;

      for(int i=1;i<n-1;++i){
         A[i][0]=1 ; A[i][1]=-2 ; A[i][2]=1;
         Q[i]=q; 
      }


//GAUSSIAN ELIMINATION :
//FORWARD SUBSTITUTION to reduce to UPPER TRIANGULAR MATRIX
    for(int i=1;i<n;++i){
      A[i][1]=  A[i][1]  -  (A[i][0]/A[i-1][1])*A[i-1][2];
      Q[i]=Q[i]  -   (A[i][0]/A[i-1][1])*Q[i-1];
    }


//BACKWRARD SUBSTITUTION for SOLUTION
    sol[n-1]= Q[n-1]/A[n-1][1];
    for(int i=n-2;i>-1;--i){
      sol[i]= (Q[i]- A[i][2]*sol[i+1])/A[i][1];   
    }



cout<<"TEMPERATURE AT THE NODES ARE (degrees Centigrade):";
    for(int i=0;i<n;i++){
      cout<<sol[i]-273.15<<',';
    }
    cout<<T_s-273.15;

}


