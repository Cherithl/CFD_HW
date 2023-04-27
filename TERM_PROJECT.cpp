#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<chrono>
#include<iomanip>
#include<utility>
using namespace std;
using namespace std::chrono;
typedef vector<vector<double>> matrix;


//FUNCTIONS IN ALGORITHM
void FV_Momentum();
void Pressure_Poisson();
void Velocity_Correction();

//FUNCTIONS FOR OPERATIONS
std::tuple<matrix,matrix> Update_Boundaries( matrix , matrix );
void printmatrix( matrix );

//DOMAIN PARAMETERS
const double lx=1 , ly =1;
int nx =6 , ny =6;
int n = (nx+2)*(ny+2);
const double dx = lx/nx , dy = ly/ny;

//PHYSICAL PARAMETERS
double Re=100 ;


int TIME_ITER=0;

//TIME STEP
double dt = 0.02;


//VELOCITY AND PRESSURE MATRICES
matrix u_nm1(ny+2,vector<double>(nx+2,0.0));
matrix v_nm1(ny+2,vector<double>(nx+2,0.0));
matrix u_n(ny+2,vector<double>(nx+2,0.0));
matrix v_n(ny+2,vector<double>(nx+2,0.0));
matrix u_tilda(ny+2,vector<double>(nx+2,0.0));
matrix v_tilda(ny+2,vector<double>(nx+2,0.0));

matrix P_corr(ny+2,vector<double>(nx+2,0.0));




int main(){
//ENFORCING BOUNDARY CONDITIONS    
    for( int i=1 ; i<ny+1 ; ++i ){
        u_nm1[i][1] = 0;    v_nm1[i][0]   = 2*0 -  v_nm1[i][1];   //LEFT BOUNDARY 
        u_nm1[i][nx+1]= 0;  v_nm1[i][nx+1]= 2*0 - v_nm1[i][nx];  //RIGHT BOUNDARY
    }
    for(int j=1 ; j<nx+1 ; ++j){
        u_nm1[ny+1][j] = 2*1.0 - u_nm1[ny][j] ;  v_nm1[ny+1][j] = 0; //TOP WALL
        u_nm1[0][j]    = 2*0 - u_nm1[1][j];      v_nm1[0][j]    = 0; //BOTTOM WALL
    }
    u_n = u_nm1; v_n = v_nm1; //for the time being so that we get EXPLCIT TIME INTEGRATION in the first time integration step

double err_u =1 ; double err_v=1 ;
double vel_tol = -8;

    while( log10(err_u) > vel_tol && log10(err_v) > vel_tol && TIME_ITER<2){
        cout<<" I AM HERE MF "<<endl<<endl<<endl;
        err_u=0; err_v=0;
        FV_Momentum();
        Pressure_Poisson();
        Velocity_Correction();
        for( int i=1 ; i<ny+1 ; ++i){
            for( int j=3 ; j<nx+1 ; ++j){
                err_u = pow( u_n[i][j] - u_nm1[i][j] , 2 )/(nx*ny);
            }
        }
        for( int i=2 ; i<ny+1 ; ++i){
            for( int j=1 ; j<nx+1 ; ++j){
                err_v = pow( v_n[i][j] - v_nm1[i][j] , 2 )/(nx*ny);
            }
        }

    err_u = pow( err_u , 0.5 ); err_v = pow( err_v , 0.5 );
    TIME_ITER++;
    }

}


std::tuple<matrix,matrix> Update_Boundaries( matrix u , matrix v ){
       for( int i=1 ; i<ny+1 ; ++i ){
        u[i][1] = 0;    v[i][0]   = 2*0 -  v[i][1];   //LEFT BOUNDARY 
        u[i][nx+1]= 0;  v[i][nx+1]= 2*0 - v[i][nx];  //RIGHT BOUNDARY
    }
    for(int j=1 ; j<nx+1 ; ++j){
        u[ny+1][j] = 2*1.0 - u[ny][j] ;  v[ny+1][j] = 0; //TOP WALL
        u[0][j]    = 2*0 - u[1][j];      v[0][j]    = 0; //BOTTOM WALL
    }   
    return std::make_tuple( u , v  ) ; 
}

void FV_Momentum(){
//X-MOMENTUM TIME INTEGRATION
printmatrix(u_n);cout<<endl<<endl;
        for( int i=1 ; i<ny+1 ; ++i){
            for( int j=2 ; j<nx+1 ; ++j){
                double ue_nm1 = 0.5*( u_nm1[i][j] + u_nm1[i][j+1] ),     ue_n = 0.5*( u_n[i][j] + u_n[i][j+1] );
                double uw_nm1 = 0.5*( u_nm1[i][j-1] + u_nm1[i][j] ),     uw_n = 0.5*( u_n[i][j-1] + u_n[i][j] );
                double un_nm1 = 0.5*( u_nm1[i][j] + u_nm1[i+1][j] ),     un_n = 0.5*( u_n[i][j] + u_n[i+1][j] );
                double us_nm1 = 0.5*( u_nm1[i-1][j] + u_nm1[i][j] ),     us_n = 0.5*( u_n[i-1][j] + u_n[i][j] );
                double vn_nm1 = 0.5*( v_nm1[i+1][j-1] + v_nm1[i+1][j]),  vn_n = 0.5*( v_n[i+1][j-1] + v_n[i+1][j]);
                double vs_nm1 = 0.5*( v_nm1[i][j-1] + v_nm1[i][j]),      vs_n = 0.5*( v_n[i][j-1] + v_n[i][j]);

                double C_nm1 = (1.0/(dx*dy))*( -( ue_nm1*ue_nm1 - uw_nm1*uw_nm1 )*dy -(un_nm1*vn_nm1 - us_nm1*vs_nm1)*dx ) ; 
                double C_n   = (1.0/(dx*dy))*( -( ue_n*ue_n - uw_n*uw_n )*dy -(un_n*vn_n - us_n*vs_n)*dx ) ;
                double D_nm1 = ( 1.0/Re) *( 1.0/(dx*dy) )*(  (u_nm1[i][j+1] -2*u_nm1[i][j] + u_nm1[i][j-1] )*(dy/dx)  + ( u_nm1[i+1][j] -2*u_nm1[i][j] +u_nm1[i-1][j] )*(dx/dy) ); 
                double D_n   = ( 1.0/Re )*( 1.0/(dx*dy) )*(  (u_n[i][j+1] -2*u_n[i][j] + u_n[i][j-1] )*(dy/dx)  + ( u_n[i+1][j] -2*u_n[i][j] +u_n[i-1][j] )*(dx/dy) );

            double  f_nm1 = C_nm1 + D_nm1 ; double f_n = C_n + D_n ;
            u_tilda[i][j] = u_n[i][j] + dt*( (1.5)*f_n + (0.5)*f_nm1 ) ; // ADAMS BASHFORTH 2nd order Integration

            }
        }
printmatrix(u_tilda);cout<<endl<<endl;
//Y-MOMENTUM TIME INTEGRATION
        for( int i=2 ; i<ny+1 ; ++i){
            for( int j=1 ; j<nx+1 ; ++j){
                double ve_nm1 = 0.5*( v_nm1[i][j] + v_nm1[i][j+1] ),     ve_n = 0.5*( v_n[i][j] + v_n[i][j+1] );
                double vw_nm1 = 0.5*( v_nm1[i][j-1] + v_nm1[i][j] ),     vw_n = 0.5*( v_n[i][j-1] + v_n[i][j] );
                double vn_nm1 = 0.5*( v_nm1[i][j] + v_nm1[i+1][j] ),     vn_n = 0.5*( v_n[i][j] + v_n[i+1][j] );
                double vs_nm1 = 0.5*( v_nm1[i-1][j] + v_nm1[i][j] ),     vs_n = 0.5*( v_n[i-1][j] + v_n[i][j] );
                double ue_nm1 = 0.5*( u_nm1[i][j+1] + u_nm1[i-1][j+1] ), ue_n = 0.5*( u_n[i][j+1] + u_n[i-1][j+1] );
                double uw_nm1 = 0.5*( u_nm1[i][j] + u_nm1[i-1][j] ),     uw_n = 0.5*( u_n[i][j] + u_n[i-1][j] );

                double C_nm1 = (1.0/(dx*dy))*( -( ve_nm1*ue_nm1 - vw_nm1*uw_nm1 )*dy -(vn_nm1*vn_nm1 - vs_nm1*vs_nm1)*dx ) ;
                double C_n   = (1.0/(dx*dy))*( -( ve_n*ue_n - vw_n*uw_n )*dy -(vn_n*vn_n - vs_n*vs_n)*dx ) ;
                double D_nm1 = (1.0/(dx*dy))*(  (v_nm1[i][j+1] -2*v_nm1[i][j] + v_nm1[i][j-1] )*(dy/dx)  + ( v_nm1[i+1][j] -2*v_nm1[i][j] +v_nm1[i-1][j] )*(dx/dy) ); //CHECK THIS
                double D_n   = ( 1.0/Re )*( 1.0/(dx*dy) )*(  (v_n[i][j+1] -2*v_n[i][j] + v_n[i][j-1] )*(dy/dx)  + ( v_n[i+1][j] -2*v_n[i][j] +v_n[i-1][j] )*(dx/dy) );

            double f_nm1 = C_nm1 + D_nm1 ; double f_n = C_n + D_n ;
            v_tilda[i][j] = v_n[i][j] + dt*( (1.5)*f_n + (0.5)*f_nm1 ); // ADAMS BASHFORTH 2nd order Integration

            }
        }
    u_nm1 = u_n; v_nm1 = v_n; //SWAPPING (n-1) st step
    
}


void Pressure_Poisson(){
//BOUNDARY CONDITIONS 
matrix S(ny+2,vector<double>(nx+2,0.0));
matrix res( ny+2,vector<double>(nx+2,0.0));
    for( int i=1 ; i<ny+1 ; ++i){
        for( int j=1 ; j<nx+1 ; ++j){
            S[i][j] =   (1.0/dt)*( ( u_tilda[i][j+1] - u_tilda[i][j] )/dx + ( v_tilda[i+1][j] - v_tilda[i][j] )/dy );
        }
    }
double AE = 1.0/pow(dx,2), AW =AE;
double AN = 1.0/pow(dy,2), AS =AN;
double AP = AE+AW+AN+AS;
double tol = -5; double res_avg =1;
        while ( log10(res_avg) > tol ){
            for( int i=1 ; i<ny+1 ; ++i ){
                P_corr[i][0] = P_corr[i][1];
                P_corr[i][nx+1] = P_corr[i][nx];
            }    

            for( int j=1 ; j<nx+1 ; ++j){
                P_corr[0][j] = P_corr[1][j];
                P_corr[nx+1][j] = P_corr[nx][j];
            }
        //GAUSS SEIDEL 
                for( int i=1 ; i<ny+1 ; ++i){
                    for( int j=1; j<nx+1 ; ++j){
                        P_corr[i][j] = (1.0/AP)*( AW*P_corr[i][j-1] + AE*P_corr[i][j+1] + AN*P_corr[i+1][j] + AS*P_corr[i-1][j] - S[i][j]  ) ;
                    }
                }
            res_avg=0;      
                for( int i=1; i<ny+1 ; ++i){
                    for( int j=1 ; j<nx+1 ; ++j){
                        res[i][j] = AP*P_corr[i][j] - ( AW*P_corr[i][j-1] + AE*P_corr[i][j+1] + AN*P_corr[i+1][j] + AS*P_corr[i-1][j] - S[i][j]  );
                        res_avg += pow( res[i][j] ,2)/(nx*ny);
                    }
                }
            res_avg = pow(res_avg,0.5);
        //cout<<res_avg<<endl;
        }

}

//CHECK n-1 n swapping and also error calculations
void Velocity_Correction(){
    for( int i=1 ; i< ny+1 ; ++i){
        for( int j=2 ; j<nx+1 ; ++j){
            u_n[i][j] = u_tilda[i][j] - dt*( P_corr[i][j] - P_corr[i][j-1] )/dx;
        }
    }
    for( int i=2 ; i<ny+1 ; ++i){
        for( int j=1 ; j<nx+1 ; ++j){
            v_n[i][j] = v_tilda[i][j] - dt*( P_corr[i][j] - P_corr[i-1][j] )/dy;
        }
    }

}


void printmatrix( matrix mat){
    int Rows = mat.size();
    int Cols = mat[0].size();
    for( int i=0 ; i<Rows ; ++i){
        for( int j=0 ; j<Cols ; ++j){
            cout<< mat[i][j] <<'\t';
            if( j == Cols-1 ) cout<<endl;
        }
    }
}