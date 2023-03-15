#include<iostream>
#include<cmath>
#include <fstream>
#include <chrono>
#include<iomanip>
using namespace std::chrono;
using namespace std;

auto start =high_resolution_clock::now();
void printmatrix(int);
void init();
void GS();
//INPUTS
  const double lx=1;const double ly=1;                      //dimensions of the 2D Domain
  const int nx=81; const int ny=nx;           //No of nodes in x and y direction
  const int n=nx*ny;                         //No of nodes
  double w;                          //RELAXATION FACTOR

// Calculated Parameters
const double dx=lx/(nx-1); const double dy=ly/(ny-1);

//DECLARATION OF NODAL TEMPERATURES
double X[nx][ny],Y[nx][ny];
double T[nx][ny]={}; 
double S[nx][ny]={};
double res[nx][ny]={};

int main(){
    w=1;
    init();
}

void init(){   
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
//BOUNDARY VALUES DECLARATION
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
//SOURCE TERMS 
        for(int i=1;i<ny-1;++i){
            for(int j=1;j<nx-1;++j){
                S[i][j]=50000*(100*( pow(1-X[i][j],2)+ pow(Y[i][j],2))-2)*exp(-50*( pow(1-X[i][j],2)+pow(Y[i][j],2)));
            }
        }
        GS();
            auto stop = high_resolution_clock::now();
            auto duration =duration_cast<milliseconds>(stop -start);
            cout<<'\n'<<"The time taken for convergence is : "<<duration.count() <<"\n";
}


//************************************************************************************************************************************************************************************//
void GS(){
        //cout<<"Start of Gauss Seidel"<<endl;
        double tol=1e-8;
        double res_max=1; double res_avg=1;
        int iter=0;

        ofstream fout("HW5_GS_residuals_81.txt");
        fout<<"ITERATION:"<<"\t"<<"log(MAX RESIDUAL) :"<<"\t\t"<<"AVG RESIDUAL"<<"\t\t"<<"AVG ERROR"<<'\n'; 
    while((res_avg)>tol){
                res_max=0; res_avg=0;                    //INITIALISING MAX RESIDUAL FOR THE RESPECTIVE ITERATION
                double err=0; // INITIALISING AVERAGE ERROR AND AVERAGE RESIDUAL 
                for(int i=1;i<ny-1;++i){
                    for(int j=1;j<nx-1;++j){
                        double dummy=0.25*( 4*(1-w)*T[i][j] + w*( T[i][j+1] + T[i][j-1] + T[i+1][j] + T[i-1][j] - dx*dy*(S[i][j]) ));
                        err+=pow(T[i][j]-dummy,2); T[i][j]=dummy; //ERROR CALCULATION    
                    }
                }
                err=pow((err/n),0.5);
    //RESIDUAL AND ERROR CALCULATION    
                for(int i=1;i<ny-1;++i){
                    for(int j=1;j<nx-1;++j){
                    res[i][j]= abs(T[i][j+1]+T[i][j-1]+T[i+1][j]+T[i-1][j]-4*T[i][j]-dx*dy*S[i][j]);
                    res_avg+=pow(res[i][j],2)/((nx-2)*(ny-2));
                    if(res[i][j]>res_max){ res_max=res[i][j];}
                    }
                }
            res_avg=pow((res_avg),0.5);
            ++iter;
    fout<<iter<<"\t"<<res_max<<"\t"<<res_avg<<"\t"<<err<<'\n'<<setprecision(15);
            
    }
    
            fout.close();
            cout<<"Gauss Seidel Done:"<<endl;
            printmatrix(iter);
}

//*******************************************************************************************************************************************************************************************//
void printmatrix(int iter){

            cout<<"No.Iterations taken to converge :"<<iter<<'\n';
            cout<<"The NODAL VLAUES are printed to : HW5_GS_sol_81.txt"<<'\n';
            ofstream fout("HW5_GS_sol_81.txt");
            fout<<"X coord:"<<"\t\t"<<"Y coord:"<<"\t\t"<<"Temperature:"<<'\n';
            fout<<"-----------------------------------------------------------"<<'\n';
            for(int i=0;i<ny;++i){
                for(int j=0;j<nx;++j){
                    fout<<X[i][j]<<"\t"<<Y[i][j]<<"\t"<<T[i][j]<<'\n';
                }
            }

            fout.close();

}
//********************************************************************** END OF GAUSS SEIDEL ********************SCHEME************************************************************************// 
