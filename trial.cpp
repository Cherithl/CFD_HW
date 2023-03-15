#include<iostream>
#include<cmath>

using namespace std;

int main(){

    const double lx=1;const double ly=1;                      //dimensions of the 2D Domain
  const int nx=21; const int ny=nx;           //No of nodes in x and y direction
  const int n=nx*ny;                         //No of nodes
double D[nx-2]={};double W[nx-2]={}; double E[nx-2]={};
            D[0]=4;     E[0]=-1;     
            W[nx-3]=-1;   D[nx-3]=4; 
            for(int i=1;i<nx-3;++i){
            E[i]=W[i]=-1; D[i]=4; 
            }
    double Q[nx-2]={};
for(int j=0;j<nx-2;++j){
        if(j==0)
        Q[j]=S[i][j+1]*dx*dy+T[i-1][j+1]+T[i+1][j+1]+T[i][j];
        else if(j==nx-3)
        Q[j]=S[i][j+1]*dx*dy+T[i-1][j+1]+T[i+1][j+1]+T[i][j+2];
        else{
            Q[j]=S[i][j+1]*dx*dy+T[i-1][j+1]+T[i+1][j+1];//VERIFY AGAIN
            //cout<<T[i-1][j+1]<<'\t';
            }

    }
 
        
    for(int j=1;j<nx-2;++j){
      D[j]=  D[j]  -  (W[j]/D[j-1])*E[j-1];
      Q[j]=Q[j]    -  (W[j]/D[j-1])*Q[j-1];
    }

    for(int i=0;i<nx-2;++i){
        cout<<E[i]<<endl
    }
    for(int i=0;i<nx-2;++i){
        cout<<D[i]<<endl
    }
    for(int i=0;i<nx-2;++i){
        cout<<W[i]<<endl
    }
    for(int i=0;i<nx-2;++i){
        cout<<Q[i]<<endl
    }


//BACKWRARD SUBSTITUTION for SOLUTION
    //cout<<(T[i][nx-2]=Q[nx-3]/D[nx-3])<<'\t';
    for(int j=nx-3;j>-1;--j){
        
      T[i][j]= (Q[j]- E[j]*T[i][j+1])/D[j];
    }
}