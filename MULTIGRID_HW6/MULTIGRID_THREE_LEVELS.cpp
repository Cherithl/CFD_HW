#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include <chrono>
#include<iomanip>
#include <utility>
using namespace std;
using namespace std::chrono;
typedef vector<vector<double>> matrix;

auto start =high_resolution_clock::now();

//INPUTS
const double lx=1; const double ly=1;  //dimensions of Domain
const int nx=81; const int ny=nx;     //Number of nodes along X and Y
const int n=nx*ny;                   //Total number of Nodes
const int w=1;

//DERIVED INPUTS
const double dx=lx/(nx-1); const double dy=ly/(ny-1);

//DECLARING ADDITIONAL MATRICES
double X_h[ny][nx],Y_h[ny][nx];
double S[ny][nx]={};

matrix T(ny,vector<double>(nx,0.0));

//FUNCITONS DECLARATIONS..............................................................................................................................................................
matrix Smooth();
std::tuple<matrix,matrix> Coarse_Smooth( matrix , matrix , int);
matrix Restrict(matrix);
matrix Prolongate(matrix);
matrix Correct(matrix,matrix);
matrix init(matrix);

void Multigrid();
void printmatrix();

//GLOBAL VARIABLES
double res_avg=1;
double total_iter=0;
double wus =0;
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
            cout<<"TOAL WORK UNITS TAKEN :"<<wus<<endl;

}

//**************************************************** END OF INT MAIN()*******************************************************************************************************//

void Multigrid(){
    ofstream fout("MULTIGRID_THREE_LEVELS_residuals.txt");
    fout<<"Residual:"<<"\t\t"<<"TOTAL FINE ITERATIONS"<<"\t\t"<<"TOTAL COARSE ITERATIONS"<<'\n';
    fout<<"-----------------------------------------------------------"<<'\n';

    matrix res_h( ny, vector<double>( nx , 0.0) );
    matrix res_H_1 ( (ny+1)/2 , vector<double>( (nx+1)/2 , 0.0 ) ) ;
    matrix res_coarse_1 ( (ny+1)/2 , vector<double>( (nx+1)/2 , 0.0 ) ) ;
    matrix res_H_2 ( ((ny+1)/2+1)/2, vector<double>( ((ny+1)/2+1)/2 , 0.0 )   );
    matrix res_coarse_2( ((ny+1)/2+1)/2, vector<double>( ((ny+1)/2+1)/2 , 0.0 ) );
    matrix res_H_3 ( (((ny+1)/2+1)/2+1)/2, vector<double>( (((ny+1)/2+1)/2+1)/2 , 0.0 )  );
    matrix res_coarse_3( (((ny+1)/2+1)/2+1)/2, vector<double>( (((ny+1)/2+1)/2+1)/2 , 0.0 ) );

    matrix E_h( ny, vector<double>( nx , 0.0) );
    matrix E_h_1( (ny+1)/2 , vector<double>( (nx+1)/2 , 0.0 ) );
    matrix E_h_2( ((ny+1)/2+1)/2, vector<double>( ((nx+1)/2+1)/2 , 0.0 ) );
    matrix E_H_1( (ny+1)/2 , vector<double>( (nx+1)/2 , 0.0 ) ) ;
    matrix E_H_2( ((ny+1)/2+1)/2, vector<double>( ((nx+1)/2+1)/2 , 0.0 ) );
    matrix E_H_3( (((ny+1)/2+1)/2+1)/2, vector<double>( (((ny+1)/2+1)/2+1)/2 , 0.0 ) );
    
    double tol=-8;

cout<<"Multigrid Start"<<endl;
cout<<"PreSmoothening..\n";

res_h = Smooth();
 while( log10(res_avg) > tol ){ //   total_iter<3

//INITIALISING MATRICES FOR THE CYCLE
     E_H_1 = init(E_H_1); E_H_2 = init(E_H_2) ; E_H_3 = init(E_H_3);
//START OF THE CYCLE 
        res_H_1 =  Restrict(res_h);
    std::tie(res_coarse_1,E_H_1)   =  Coarse_Smooth( E_H_1, res_H_1 , 1); 
        res_H_2 = Restrict( res_coarse_1 );
    std::tie(res_coarse_2,E_H_2)   =  Coarse_Smooth( E_H_2, res_H_2, 2 ); 
        res_H_3 = Restrict( res_coarse_2 );
    std::tie(res_coarse_3,E_H_3)   =  Coarse_Smooth( E_H_3, res_H_3, 3 );
        E_h_2 =  Prolongate( E_H_3 );
        E_H_2 =  Correct( E_H_2 , E_h_2 );
    std::tie(res_coarse_2,E_H_2)   =  Coarse_Smooth( E_H_2, res_H_2, 2 );   
        E_h_1 =  Prolongate(E_H_2);
        E_H_1 =  Correct( E_H_1 , E_h_1 );
    std::tie(res_coarse_1,E_H_1)  =  Coarse_Smooth( E_H_1, res_H_1 ,1 );       
         E_h   =  Prolongate(E_H_1);
         T     =  Correct( T , E_h );

    res_h =  Smooth();
 fout<<res_avg<<"\t"<<total_iter_fine<<"\t"<<total_iter_coarse<<'\n';
 }
    fout.close();
    cout<<res_avg<<endl;
    printmatrix();
    cout<<"End of Multigrid"<<endl;
 }






//GAUSS SEIDEL SMOOOTHENING-------------------------------------------------------------------------------------------------------------------------------------------------------------------//
matrix Smooth(){
        matrix res(ny , vector<double>( nx , 0.0 ));
        double res_new=1; double res_old=1;
        int iter=0; double tol=-8;
    while(  res_new/res_old < 0.5 || iter < 2  ){   //    log10(res_new) > tol
        res_old = res_new; res_new=0;
                for(int i=1;i<ny-1;++i){
                    for(int j=1;j<nx-1;++j){
                        T[i][j]=0.25*( 4*(1-w)*T[i][j] + w*( T[i][j+1] + T[i][j-1] + T[i+1][j] + T[i-1][j] - dx*dy*(S[i][j]) ));
                    }
                }                           
    //RESIDUAL  CALCULATION
                for(int i=1;i<ny-1;++i){
                    for(int j=1;j<nx-1;++j){    
                    res[i][j] = -( (T[i][j+1]+T[i][j-1]+T[i+1][j]+T[i-1][j]-4*T[i][j])/(dx*dy)    -S[i][j]); //removed the negative sign
                    res_new += pow( res[i][j] , 2 )/((nx-2)*(ny-2));
                    }
                }
                res_new = pow( res_new , 0.5 ) ;
            ++iter;
            ++total_iter;
            ++total_iter_fine;
            ++wus;
    }
    res_avg=res_new;
    return res;

}


//COARSE LEVEL SMOOTHENING -------------------------------------------------------------------------------------------------------------------------------------------------------------------//
std::tuple<matrix,matrix> Coarse_Smooth( matrix E_H , matrix res_H , int level ){
 
    int Rows_H = res_H.size();
    int Cols_H = res_H[0].size();
    double dx_H = dx*pow(2,level) ; 
    double dy_H = dy*pow(2,level);

//matrix E_H( Rows_H , vector<double>(Cols_H ,0.0));
matrix res_coarse( Rows_H , vector<double>(Cols_H , 0.0));   
        double res_coarse_new=1; 
        double res_coarse_old=1;
        int iter_coarse=0;
    while( res_coarse_new/res_coarse_old <0.5 || iter_coarse<1 ){  
                res_coarse_old = res_coarse_new ;                    //INITIALISING MAX RESIDUAL FOR THE RESPECTIVE ITERATION
                res_coarse_new = 0; 
                for( int i=1 ; i<Rows_H-1 ; ++i ){
                    for( int j=1 ; j< Cols_H-1 ; ++j ){
                        E_H[i][j] =0.25*( 4*(1-w)*E_H[i][j] + w*( E_H[i][j+1] + E_H[i][j-1] + E_H[i+1][j] + E_H[i-1][j] - dx_H*dy_H*(res_H[i][j]) ));
                    }
                }               

                for(int i=1 ; i<Rows_H-1 ; ++i){
                    for( int j=1 ; j<Cols_H-1 ; ++j){
                        res_coarse[i][j] = -( (E_H[i][j+1]+E_H[i][j-1]+E_H[i+1][j]+E_H[i-1][j]-4*E_H[i][j] )/(dx_H*dy_H)  -res_H[i][j])   ;                         
                        res_coarse_new += pow( res_coarse[i][j] , 2 )/( (Rows_H-2)*(Cols_H -2) );
                    }
                }
                res_coarse_new = pow(res_coarse_new , 0.5);
                //cout<<res_coarse_new<<endl;
                ++iter_coarse;
                ++total_iter;
                ++total_iter_coarse;
                wus += pow(0.5,level+1);
                //cout<<total_iter<<endl;
    }

return std::make_tuple( res_coarse , E_H  ) ;

}



//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//RESTRICTION FUNCTION

matrix Restrict(matrix res_h){
    int Rows_h = res_h.size();     
    int Cols_h = res_h[0].size(); 

    int Rows_H = (Rows_h+1)/2;
    int Cols_H = (Cols_h+1)/2; 
    
double res_max=0; //DELETE THIS MAN
    matrix res_H( Rows_H , vector<double>( Cols_H , 0.0 ) );
    
        for(int i=1;i<Rows_H-1;++i){
            for(int j=1;j<Cols_H-1;++j){
                res_H[i][j]= (4*res_h[2*i][2*j]  + res_h[2*i-1][2*j] + res_h[2*i+1][2*j] + res_h[2*i][2*j +1] + res_h[2*i][2*j -1] )/8;// (4*P + S + N + E + W)*(1/8)    
            }
        }  

return res_H;

}


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//PROLONGATION FUNCTION

matrix Prolongate(matrix E_H){
    int Rows_H = E_H.size();
    int Cols_H = E_H[0].size();

    int Rows_h=Rows_H*2-1; int Cols_h = Cols_H*2 -1 ;
matrix E_h( Rows_h, vector<double>( Cols_h , 0.0)) ; 
double E_h_max=0;
  for(int i=0 ; i<Rows_h ; ++i){
    for( int j=0 ; j<Cols_h ; ++j){
        if( i%2==0 && j%2==0 ){
            E_h[i][j] = E_H[i/2][j/2];
        }
        else if( i%2==0 && j%2!=0){
            E_h[i][j] = ( E_H[i/2][(j+1)/2] + E_H[i/2][(j-1)/2] )/2 ; 
        }
        else if( i%2!=0 && j%2==0){
            E_h[i][j] = ( E_H[(i-1)/2][j/2] + E_H[(i+1)/2][j/2] )/2 ;
        }
        else{
            E_h[i][j] = ( E_H[(i-1)/2][(j-1)/2] + E_H[(i-1)/2][(j+1)/2] + E_H[(i+1)/2][(j-1)/2] + E_H[(i+1)/2][(j+1)/2] )/4;
        }
        if( abs(E_h[i][j]) > abs(E_h_max) ){ E_h_max=E_h[i][j]; }
    }
  }     
//cout<<"MAX CORRECTION : "<<E_h_max<<endl;

return E_h;

}

//CORRECTION FUNCTION------------------------------------------------------------------------------------------------------------------------------------------------------------//
matrix Correct( matrix U , matrix E_h){

    int Rows = E_h.size();
    int Cols = E_h[0].size();

    double E_h_max=0;
    for(int i=1;i<Rows-1;++i){
        for(int j=1;j<Cols-1;++j){
            U[i][j] += E_h[i][j];
        }
    }
    return U;
 }

//----------------------------------------------------------------------------------------INITIALISING MATRICES FUNCTIONS---------------------------------------------------//

matrix init(matrix mat){
    int Rows= mat.size();
    int Cols= mat[0].size();
    for( int i=0; i<Rows ; ++i){
        for( int j=0 ; j<Cols ; ++j ){
            mat[i][j]=0;
        }
    }
    return mat;
}


//PRINTING OUT SOLUTION-----------------------------------------------------------------------------------------------------------------------------------------------------------//
void printmatrix(){

            //cout<<"No.Iterations taken to converge :"<<iter<<'\n';
            cout<<"The NODAL VLAUES are printed to : MULTIGRID_THREE_LEVELS_sol.txt"<<'\n';
            ofstream fout("MULTIGRID_THREE_LEVELS_sol.txt");
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