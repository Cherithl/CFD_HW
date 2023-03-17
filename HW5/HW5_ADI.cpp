#include<iostream>
#include<cmath>
#include <fstream>
#include <chrono>
using namespace std::chrono;
using namespace std;

auto start =high_resolution_clock::now();
void printmatrix(int);
void ADI();
void TDMAX(int);
void TDMAY(int);
void init();
//INPUTS
const double lx=1;const double ly=1;         //dimensions of the 2D Domain
const int nx=81; const int ny=nx;           //No of nodes in x and y direction
const int n=nx*ny;                         //No of nodes
double w;                         //RELAXATION FACTOR

// Calculated Parameters
const double dx=lx/(nx-1); const double dy=ly/(ny-1);

//DECLARATION OF NODAL TEMPERATURES
double X[ny][nx],Y[ny][nx];
double T[ny][nx]={}; 
double S[ny][nx]={};
double res[ny][nx]={};

int main(){
    w=1.32;
    init();
}

void init(){   

//INITIALISING
   for(int i=0;i<ny-1;++i){
    for(int j=0;j<nx-1;++j){
        T[i][j]=0;res[i][j]=0;
    }
   }

//MESHING IN 2D
        for(int i=0;i<ny;i++){
            for(int j=0;j<nx;j++){
                X[i][j]=j*dx; Y[i][j]=i*dy;
            }
        }
//FORCING BOUNDARY VALUES 
    //LEFT BOUNDARY
        for(int row=0;row<ny;++row){
            T[row][0]=500*exp(-50*(pow(Y[row][0],2)+1));
        }
    //BOTTOM BOUNDARY
        for(int col=0;col<nx;++col){
            T[0][col]=500*exp(-50*(pow(1-X[0][col],2)))+100*X[0][col];
        }
    //RIGHT BOUNDARY
        for(int row=0;row<ny;++row){
            T[row][nx-1]=500*exp(-50*pow(Y[row][nx-1],2))+100*(1-Y[row][nx-1]);
        }
    //TOP BOUNDARY
        for(int col=0;col<nx;++col){
            T[ny-1][col]=500*exp(-50*(pow(1-X[ny-1][col],2)+1));
        }
//DECLARING SOURCE TERMS
        for(int i=1;i<ny-1;++i){
            for(int j=1;j<nx-1;++j){
                S[i][j]=50000*(100*( pow(1-X[i][j],2)+ pow(Y[i][j],2))-2)*exp(-50*( pow(1-X[i][j],2)+pow(Y[i][j],2)));  
            }
        }
//CALLING ADI()
        ADI();
    auto stop = high_resolution_clock::now();
    auto duration =duration_cast<milliseconds>(stop -start);
    //cout<<'\n'<<"The time taken for convergence is : "<<duration.count() <<"\n";
}


//*******************************************************************************************************************************************************************************************//
void ADI(){
        //cout<<"Start of ADI"<<endl;
        double tol=1e-8;
        double res_max=1;
        int iter=0;
        ofstream fout("HW5_ADI_residuals_81.txt");
        fout<<"ITERATION:"<<"\t"<<"log(MAX RESIDUAL) :"<<"\t\t"<<"AVG RESIDUAL"<<"\t\t"<<"AVG ERROR"<<'\n';
    double res_avg=1;
    while(log10(res_avg)>log10(tol)){
                res_max=0;                     //INITIALISING MAX RESIDUAL FOR THE RESPECTIVE ITERATION
                double err_avg=0; res_avg=0;// INITIALISING AVERAGE ERROR AND AVERAGE RESIDUAL 
    //CALLING TDMAX AND TDMAY
               double dummy[ny][nx];
               for(int i=1;i<ny-1;++i){
                for(int j=1;j<ny-1;++j){
                    dummy[i][j]=T[i][j];
                }
               }

                for(int i=1;i<ny-1;++i){
                    TDMAX(i);
                    TDMAY(i);
                }

    //RESIDUAL AND ERROR CALCULATION    
                for(int i=1;i<ny-1;++i){
                    for(int j=1;j<nx-1;++j){
                    res[i][j]= abs( (T[i][j+1] + T[i][j-1] + T[i+1][j] + T[i-1][j] - 4*T[i][j])/(dx*dy)- S[i][j] );
                    res_avg+=pow( res[i][j] , 2 )/((nx-2)*(ny-2));
                    err_avg+=pow( dummy[i][j]-T[i][j] , 2 )/n;
                    if(res[i][j]>res_max){ res_max = res[i][j]; }
                    }
                }
                err_avg=pow( err_avg , 0.5 );
                res_avg=pow( res_avg , 0.5 );
                ++iter;
            fout<<iter<<"\t"<<res_max<<"\t"<<res_avg<<"\t"<<err_avg<<'\n';
                
    }
    cout<<w<<'\t'<<iter<<endl;
            fout.close();
            //cout<<"ADI Done:"<<endl;
            printmatrix(iter);
}
//*******************************************************************************************************************************************************************************************//
void TDMAX(int row ){
//DECLARING DIAGONAL FOR TDMA
    double D[nx-2]={};double E[nx-2]={};double W[nx-2]={};
    W[0]=0;D[0]=4;E[0]=-1*w;   W[nx-3]=-1*w;D[nx-3]=4;E[nx-3]=0;
    
    for(int i=1;i<nx-3;++i){
        W[i]=-1*w;D[i]=4;E[i]=-1*w;
    }
//DECLARING SOURCE FOR TDMA
    double Q[nx-2]={};
    for(int i=0;i<nx-2;++i){
        if(i==0)
        Q[i]=w*(-S[row][i+1]*dx*dy+ T[row+1][i+1] + T[row-1][i+1] + T[row][i])  +4*(1-w)*T[row][i+1] ;
        else if(i==nx-3)
        Q[i]=w*( -S[row][i+1]*dx*dy+ T[row+1][i+1] + T[row-1][i+1] + T[row][i+2]) +4*(1-w)*T[row][i+1] ;
        else
        Q[i]=w*(-S[row][i+1]*dx*dy+ T[row+1][i+1] + T[row-1][i+1])+4*(1-w)*T[row][i+1];
    }

//FORWARD SUBSTITUION
    for(int i=1;i<nx-2;++i){
        D[i]= D[i] - (W[i]/D[i-1])*E[i-1] ;
        Q[i]= Q[i] - (W[i]/D[i-1])*Q[i-1] ;
    }

//BACKWARD SUBSTITUTION
    for(int i=nx-3;i>-1;--i){
        T[row][i+1]= ( Q[i] - E[i]*T[row][i+2] )/D[i];
       
    }
}
//*********************************************************************************************************************************************************************************************//
void TDMAY(int col ){
//DECLARING DIAGONAL FOR TDMA
    double D[ny-2]={};double E[ny-2]={};double W[ny-2]={};
    W[0]=0;D[0]=4;E[0]=-1*w;   W[ny-3]=-1*w;D[ny-3]=4;E[ny-3]=0;
     
    for(int i=1;i<ny-3;++i){
        W[i]=-1*w;D[i]=4;E[i]=-1*w;
    }
//DECLARING SOURCE FOR TDMA  
    double Q[ny-2]={};
    for(int i=0;i<ny-2;++i){
        if(i==0)
        Q[i]= w*( -S[i+1][col]*dx*dy+ T[i][col] + T[i+1][col-1] + T[i+1][col+1] ) + 4*(1-w)*T[i+1][col]  ;
        else if(i==ny-3)
        Q[i]=w*( -S[i+1][col]*dx*dy+ T[i+2][col] + T[i+1][col-1] + T[i+1][col+1] ) + 4*(1-w)*T[i+1][col]  ;
        else
        Q[i]=w*( -S[i+1][col]*dx*dy+ T[i+1][col-1] + T[i+1][col+1] ) + 4*(1-w)*T[i+1][col]  ;
    }

//FORWARD SUBSTITUION
    for(int i=1;i<ny-2;++i){
        D[i]= D[i] - (W[i]/D[i-1])*E[i-1] ;
        Q[i]= Q[i] - (W[i]/D[i-1])*Q[i-1] ;
    }

//BACKWARD SUBSTITUTION
    for(int i=ny-3;i>-1;--i){
        T[i+1][col]= ( Q[i] - E[i]*T[i+2][col] )/D[i];

       
    }
}
//*******************************************************************************************************************************************************************************************//
void printmatrix(int iter){

          //cout<<iter<<'\n';
        //cout<<"No.Iterations taken to converge :"<<iter<<'\n';
        //cout<<"The NODAL VLAUES are printed to : HW5_ADI_sol_81.txt"<<'\n';
        //cout<<"The Residual VLAUES are printed to : HW5_ADI_residuals_81.txt"<<'\n';
        ofstream fout("HW5_ADI_sol_81.txt");
        fout<<"X coord:"<<"\t\t"<<"Y coord:"<<"\t\t"<<"Temperature:"<<'\n';
        fout<<"-----------------------------------------------------------"<<'\n';
        for(int i=0;i<ny;++i){
            for(int j=0;j<nx;++j){
                fout<<X[i][j]<<"\t"<<Y[i][j]<<"\t"<<T[i][j]<<'\n';
            }
        }

        fout.close();
}
//********************************************************************** END OF ALTERNATING DIRECTION IMPLICIT SCHEME************************************************************************//




  