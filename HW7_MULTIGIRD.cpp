#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include <chrono>
#include<iomanip>
using namespace std;
using namespace std::chrono;
typedef vector<vector<double>> matrix;

auto start =high_resolution_clock::now();

//INPUTS
const double lx=1; const double ly=1;  //dimensions of Domain
const int nx=81; const int ny=nx;     //Number of nodes along X and Y
const int n=nx*ny;                   //Total number of Nodes
int w=1;

//DERIVED INPUTS
const double dx=lx/(nx-1); const double dy=ly/(ny-1);

//DECLARING ADDITIONAL MATRICES
double X_h[ny][nx],Y_h[ny][nx];
matrix T(ny,vector<double>(nx,0.0));
double S[ny][nx]={};


//FUNCITONS DECLARATIONS..............................................................................................................................................................
matrix Smooth();
//matrix Coarse_Smooth( matrix );
//matrix Restrict(matrix);
//matrix Prolongate(matrix);
void Multigrid();
//void Correct(matrix);
void printmatrix();

//GLOBAL VARIABLES
double res_avg=1;
double total_iter=0;
double total_iter_fine=0;
double total_iter_coarse=0;



int main(){

//MESHING IN 2D
    for(int i=0;i<ny;++i){
        for(int j=0;j<nx;++j){
            X_h[i][j]=j*dx; Y_h[i][j]=i*dy;
        }
    }

//BOUNDARY VALUES DECLARATION
    //LEFT BOUNDARY
            for(int row=0;row<ny;++row){
                T[row][0]=500*exp(-50*(pow(Y_h[row][0],2)+1));
            }
    //BOTTOM BOUNDARY
            for(int col=0;col<nx;++col){
                T[0][col]=500*exp(-50*(pow(1-X_h[0][col],2)))+100*X_h[0][col];
            }
    //RIGHT BOUNDARY
            for(int row=0;row<ny;++row){
                T[row][nx-1]=500*exp(-50*pow(Y_h[row][nx-1],2))+100*(1-Y_h[row][nx-1]);
            }
    //TOP BOUNDARY
            for(int col=0;col<nx;++col){
                T[ny-1][col]=500*exp(-50*(pow(1-X_h[ny-1][col],2)+1));
            }
    //SOURCE TERMS
            for(int i=1;i<ny-1;++i){
                for(int j=1;j<nx-1;++j){
                    S[i][j]=50000*(100*( pow(1-X_h[i][j],2)+ pow(Y_h[i][j],2))-2)*exp(-50*( pow(1-X_h[i][j],2)+pow(Y_h[i][j],2)));
                }
            }

    cout<<"before Multigrid\n";
//MULTIGRID FROM HERE
Multigrid();
            auto stop = high_resolution_clock::now();
            auto duration =duration_cast<milliseconds>(stop -start);
            cout<<'\n'<<"The time taken for convergence is : "<<duration.count() <<"\n";
            cout<<"TOTAL FINE ITERATIONS :"<<total_iter_fine<<endl;
            cout<<"TOAL COARSE ITERATIONS :"<<total_iter_coarse<<endl;

}

//**************************************************** END OF INT MAIN()*******************************************************************************************************//

void Multigrid(){
    matrix res_h( ny, vector<double>( nx , 0.0) );
    matrix res_H ( (ny+1)/2 , vector<double>( (nx+1)/2 , 0.0 ) ) ;
    matrix E_h( ny, vector<double>( nx , 0.0) );
    matrix E_H( (ny+1)/2 , vector<double>( (nx+1)/2 , 0.0 ) ) ;
double tol=-8;
cout<<"Before Smooth\n";
        res_h = Smooth();
        cout<<"Start of Multigrid"<<endl;

while( log10(res_avg) > tol ){ //        log10(res_avg) > tol
       // res_H =  Restrict(res_h);        
        //E_H   =  Coarse_Smooth(res_H); 
       //E_h   =  Prolongate(E_H);
       // Correct(E_h);
        res_h =  Smooth();
}
cout<<res_avg<<endl;
printmatrix();
cout<<"End of Multigrid"<<endl;
}






//GAUSS SEIDEL SMOOOTHENING-------------------------------------------------------------------------------------------------------------------------------------------------------------------//
matrix Smooth(){
        matrix res(ny , vector<double>( nx , 0.0 ));
        double res_new=1; double res_old=1;
        int iter=0; double tol=-8;

    while( log10(res_new) > tol ){  //        res_new/res_old < 0.5 || iter <2 
        res_old = res_new; res_new=0;
                for(int i=1;i<ny-1;++i){
                    for(int j=1;j<nx-1;++j){
                        T[i][j]=0.25*( 4*(1-w)*T[i][j] + w*( T[i][j+1] + T[i][j-1] + T[i+1][j] + T[i-1][j] - dx*dy*(S[i][j]) ));
                    }
                }                           
    //RESIDUAL  CALCULATION
                for(int i=1;i<ny-1;++i){
                    for(int j=1;j<nx-1;++j){    
                    res[i][j] = -(T[i][j+1]+T[i][j-1]+T[i+1][j]+T[i-1][j]-4*T[i][j]-dx*dy*S[i][j]); //removed the negative sign
                    res_new += pow( res[i][j] , 2 )/((nx-2)*(ny-2));
                    }
                }
                res_new = pow( res_new , 0.5 ) ;
            ++iter;
            ++total_iter;
            ++total_iter_fine;
            //cout<<total_iter<<endl;
    }
    
    res_avg=res_new;
    //cout<<res_avg<<endl;
    return res;

}


//COARSE LEVEL SMOOTHENING -------------------------------------------------------------------------------------------------------------------------------------------------------------------//


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//RESTRICTION FUNCTION




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//PROLONGATION FUNCTION



//CORRECTION FUNCTION------------------------------------------------------------------------------------------------------------------------------------------------------------//


//PRINTING OUT SOLUTION-----------------------------------------------------------------------------------------------------------------------------------------------------------//
void printmatrix(){

            //cout<<"No.Iterations taken to converge :"<<iter<<'\n';
            cout<<"The NODAL VLAUES are printed to : HW6_sol.txt"<<'\n';
            ofstream fout("HW6_sol.txt");
            fout<<"X coord:"<<"\t\t"<<"Y coord:"<<"\t\t"<<"Temperature:"<<'\n';
            fout<<"-----------------------------------------------------------"<<'\n';
            for(int i=0;i<ny;++i){
                for(int j=0;j<nx;++j){
                    fout<<X_h[i][j]<<"\t"<<Y_h[i][j]<<"\t"<<T[i][j]<<'\n';
                }
            }

            fout.close();

}

//**************************************************  END OF MULTIGRID WITH GAUSS SEIDEL *******************************************************************************************//