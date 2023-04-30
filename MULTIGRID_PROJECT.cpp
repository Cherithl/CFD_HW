#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <utility>

using namespace std;
typedef vector<vector<double>> matrix;

/* FUNCTIONS FOR ALGORITHM */ /*------------------------------------------------------------*/
void FV_Momentum();
void Pressure_Poisson();
matrix Smooth();
std::tuple<matrix,matrix> Coarse_Smooth( matrix , matrix , int);
matrix Restrict(matrix);
matrix Prolongate(matrix);
matrix Correct(matrix,matrix);
matrix init(matrix);


void Velocity_Correction();
void Velocity_L2Norm();
void Output_Solution();


/* FUNCTION FOR OPERATIONS */
std::tuple<matrix,matrix> Update_Boundaries( matrix , matrix );
void printmatrix(matrix);

/*------------------------------------------------------------------------------------------*/
/*  PARAMETERS  */
int nx = 128 , ny=nx;
double lx = 1 , ly=lx ;
double dx = lx/nx , dy = ly/ny;

/* PHYSICAL PARAMETERS */
double Re = 100;
double dt = 5.0E-4;

/* PARAMETERS FOR CONVERGENCE AND REFERENCE */
int TIME_ITER =0 ; double res_avg = 1.0 ; double w=1.93;
double err_u = 1 , err_v = err_u;

/*------------------------------------------------------------------------------------------*/
/* MATRICES FOR THE PROBLEM */
matrix    u_nm1( ny+2 , vector<double>( nx+2 , 0.0) );
matrix    v_nm1( ny+2 , vector<double>( nx+2 , 0.0) );
matrix      u_n( ny+2 , vector<double>( nx+2 , 0.0 ));
matrix      v_n( ny+2 , vector<double>( nx+2 , 0.0 ));
matrix  u_tilda( ny+2 , vector<double>( nx+2 , 0.0 ));
matrix  v_tilda( ny+2 , vector<double>( nx+2 , 0.0 ));

/* MATRICES FOR THE PRESSURE POISSON PROBLEM */
matrix  P_corr( ny+2 , vector<double>( nx+2 , 0.0 ));
matrix      Ae( ny+2 , vector<double>( nx+2 , 1.0/pow(dx,2) ));
matrix      Aw( ny+2 , vector<double>( nx+2 , 1.0/pow(dx,2) ));
matrix      An( ny+2 , vector<double>( nx+2 , 1.0/pow(dy,2) ));
matrix      As( ny+2 , vector<double>( nx+2 , 1.0/pow(dy,2) ));

int main(){
    for( int i=1 ; i<ny+1 ; ++i){
        Aw[i][1]  = 0.0; 
        Ae[i][nx] = 0.0;
    }
    for( int j=1 ; j<nx+1 ; ++j){
        An[ny][j] = 0.0;
        As[1][j] = 0.0;
    }
    std::tie( u_n , v_n ) = Update_Boundaries( u_n , v_n );
    u_nm1 = u_n ; v_nm1 = v_n;
    double vel_tol = 1E-8 ; 
    while( err_u > vel_tol || TIME_ITER < 4 ){
        ++TIME_ITER;
        FV_Momentum();
        Pressure_Poisson();
        Velocity_Correction();
        Velocity_L2Norm();
    cout<<" TIME ITER :"<<TIME_ITER<<endl<<'\t';    
    cout<<" err_u : "<<err_u<<endl<<'\t';
    cout<<" err_v : "<<err_u<<endl;
    }
    printmatrix( v_n ); 
    cout<<endl<<endl;
    printmatrix( u_n );
    Output_Solution();

}


std::tuple<matrix,matrix> Update_Boundaries( matrix u , matrix v ){
    for( int i=1 ; i<ny+1 ; ++i){
        u[i][1]    = 0.0; v[i][0]    = 2*0.0 -  v[i][1]; //0
        u[i][nx+1] = 0.0; v[i][nx+1] = 2*0.0 - v[i][nx]; //0
    }
    for( int j=1 ; j<nx+1 ; ++j){
        u[0][j]    = 2*0.0 -  u[1][j]; v[0][j]    = 0.0; //0
        u[ny+1][j] = 2*1.0 - u[ny][j]; v[ny+1][j] = 0.0;
    }

    return std::tuple(u,v);
}

void FV_Momentum(){
    /* X-MOMENTUM Time Integration */
    for( int i=1 ; i<ny+1 ; ++i ){
        for( int j=2 ; j<nx+1 ; ++j ){
            double ue_n = 0.5*( u_n[i][j] + u_n[i][j+1] );                   double ue_nm1 = 0.5*( u_nm1[i][j] + u_nm1[i][j+1] );
            double uw_n = 0.5*( u_n[i][j-1] + u_n[i][j] );                   double uw_nm1 = 0.5*( u_nm1[i][j-1] + u_nm1[i][j] );
            double un_n = 0.5*( u_n[i][j] + u_n[i+1][j] );                   double un_nm1 = 0.5*( u_nm1[i][j] + u_nm1[i+1][j] );
            double us_n = 0.5*( u_n[i-1][j] + u_n[i][j] );                   double us_nm1 = 0.5*( u_nm1[i-1][j] + u_nm1[i][j] );
            double vn_n = 0.5*( v_n[i+1][j-1] + v_n[i][j] );                 double vn_nm1 = 0.5*( v_nm1[i+1][j-1] + v_nm1[i][j] );
            double vs_n = 0.5*( v_n[i][j-1] + v_n[i][j] );                   double vs_nm1 = 0.5*( v_nm1[i][j-1] + v_nm1[i][j] );

        double Conv_n = -( ue_n*ue_n - uw_n*uw_n )/dx -( un_n*vn_n - us_n*vs_n )/dy;  
        double Conv_nm1 = -( ue_nm1*ue_nm1 - uw_nm1*uw_nm1 )/dx -( un_nm1*vn_nm1 - us_nm1*vs_nm1 )/dy ;
        double Diff_n = (1.0/Re)*( ( u_n[i][j+1] -2*u_n[i][j] + u_n[i][j-1] )/pow( dx ,2 ) + ( u_n[i+1][j] -2*u_n[i][j] + u_n[i-1][j] )/pow( dy ,2 ) ) ;   
        double Diff_nm1 = (1.0/Re)*( ( u_nm1[i][j+1] -2*u_nm1[i][j] + u_nm1[i][j-1] )/pow( dx ,2 ) + ( u_nm1[i+1][j] -2*u_nm1[i][j] + u_nm1[i-1][j] )/pow( dy ,2 ) ) ;
        u_tilda[i][j] = u_n[i][j] + dt*( 1.5*( Conv_n + Diff_n ) + 0.5*( Conv_nm1 + Diff_nm1 ) );
        }
    }
    
    /* Y-MOMENTUM Time Integration */
    for( int i=2 ; i<ny+1 ; ++i ){
        for( int j=1 ; j<nx+1 ; ++j){
            double ve_n = 0.5*( v_n[i][j] + v_n[i][j+1] );                     double ve_nm1 = 0.5*( v_nm1[i][j] + v_nm1[i][j+1] );
            double vw_n = 0.5*( v_n[i][j-1] + v_n[i][j] );                     double vw_nm1 = 0.5*( v_nm1[i][j-1] + v_nm1[i][j] );
            double vn_n = 0.5*( v_n[i][j] + v_n[i+1][j] );                     double vn_nm1 = 0.5*( v_nm1[i][j] + v_nm1[i+1][j] );
            double vs_n = 0.5*( v_n[i-1][j] + v_n[i][j] );                     double vs_nm1 = 0.5*( v_nm1[i-1][j] + v_nm1[i][j] );
            double ue_n = 0.5*( u_n[i-1][j+1] + u_n[i][j+1] );                 double ue_nm1 = 0.5*( u_nm1[i-1][j+1] + u_nm1[i][j+1] );
            double uw_n = 0.5*( u_n[i-1][j] + u_n[i][j] );                     double uw_nm1 = 0.5*( u_nm1[i-1][j] + u_nm1[i][j] );

        double Conv_n = -( ue_n*ve_n - uw_n*vw_n )/dx - ( vn_n*vn_n - vs_n*vs_n )/dy ;
        double Conv_nm1 = -( ue_nm1*ve_nm1 - uw_nm1*vw_nm1 )/dx - ( vn_nm1*vn_nm1 - vs_nm1*vs_nm1 )/dy ;
        double Diff_n = (1.0/Re)*( (v_n[i][j+1] - 2*v_n[i][j] + v_n[i][j-1])/pow( dx ,2 ) + ( v_n[i+1][j] -2*v_n[i][j] + v_n[i-1][j])/pow( dy , 2 ) );
        double Diff_nm1 = (1.0/Re)*( (v_nm1[i][j+1] - 2*v_nm1[i][j] + v_nm1[i][j-1])/pow( dx ,2 ) + ( v_nm1[i+1][j] -2*v_nm1[i][j] + v_nm1[i-1][j])/pow( dy , 2 ) ) ;
        v_tilda[i][j] = v_n[i][j] + dt*( 1.5*( Conv_n + Diff_n ) + 0.5*( Conv_nm1 + Diff_nm1 ) );
        }
    }

    u_nm1 = u_n ;
    v_nm1 = v_n ;
    std::tie( u_tilda , v_tilda ) = Update_Boundaries( u_tilda , v_tilda );
      
}


/*========================================================PRESSURE EQUATION SOLVER=================================================================================================== */
void Pressure_Poisson(){
    /* PRESSURE Initialization */
    for( int i=1 ; i<ny+1 ; ++i ){
        for( int j=1 ; j<nx+1 ; ++j ){
            P_corr[i][j] = 0.0;
        }
    }

    matrix res_h( ny, vector<double>( nx , 0.0) );

res_h = Smooth();

    
}



/* -------------------===============================FINE SMOOTHING=============================-----------------------*/
matrix Smooth(){
    matrix  S( ny+2 , vector<double>( nx+2 , 0.0 ) );
    matrix res( ny+2, vector<double>( nx+2 , 0.0 ) );
    /* SOURCE TERMS */
    for( int i=1 ; i<ny+1 ; ++i ){
        for( int j=1 ; j<nx +1 ; ++j ){
            S[i][j] = (1.0/dt)*( ( u_tilda[i][j+1] - u_tilda[i][j] )/dx + ( v_tilda[i+1][j] - v_tilda[i][j] )/dy );
        }
    }    
    /* Solving Discreet POISSON Equation */
double tol = 1E-5; int ITER=0; res_avg =1.0 ;
    while( res_avg>tol ){
        ++ITER;
        //cout<<" PRESSURE ITER : "<<ITER<<endl;
        for( int i=1 ; i<ny+1 ; ++i ){
            for( int j=1 ; j<nx+1 ; ++j){
                double AE = Ae[i][j] , AW = Aw[i][j] , AN = An[i][j] , AS = As[i][j] ;
                double AP = AE + AW + AN + AS ;
                P_corr[i][j] = (1-w)*P_corr[i][j] + w*(1.0/AP)*( AE*P_corr[i][j+1] + AW*P_corr[i][j-1] + AN*P_corr[i+1][j] + AS*P_corr[i-1][j] - S[i][j] ) ;
            }
        }
        res_avg = 0.0;
        for( int i=1 ; i<ny+1 ; ++i ){
            for( int j=1 ; j<nx+1 ; ++j ){
                double AE = Ae[i][j] , AW = Aw[i][j] , AN = An[i][j] , AS = As[i][j] ;
                double AP = AE + AW + AN + AS ;
                res[i][j]    =  AP*P_corr[i][j] -( AE*P_corr[i][j+1] + AW*P_corr[i][j-1] + AN*P_corr[i+1][j] + AS*P_corr[i-1][j] - S[i][j] );     
                res_avg +=  pow( res[i][j] , 2 )/(nx*ny);
            }
        }
        res_avg = pow( res_avg , 0.5 );

    }
    cout<<endl<<endl;
    cout<<" ITER "<<ITER<<'\t'; cout<<" res_avg "<<res_avg<<endl;
    return res;
}


/* ---------------------------=========================RESTRICTION==============================-------------------*/
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

/* --------------------------------========================== COARSE LEVEL SMOOTHENING =========================----------------------------------------- */
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
    }

return std::make_tuple( res_coarse , E_H  ) ;

}

/* -------------------------============================PROLONGATION==============================--------------*/
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

/* --------------==============================CORRECTION===================================---------------------*/
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
/*----------------------------------------INITIALIZATION OF ERROR VECTORS------------------------------------------*/
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


/*==================================================VELOCITY CORRECTIONS AND ERRORS DONT GO DOWN BRO=========================================================================================*/
void Velocity_Correction(){
    for( int i=1 ; i<ny+1 ; ++i ){
        for( int j=2 ; j<nx+1 ; ++j ){
            u_n[i][j] = u_tilda[i][j] - dt*( P_corr[i][j] - P_corr[i][j-1]  )/dx ; 
        }
    }
    for( int i=2 ; i<ny+1 ; ++i ){
        for( int j=1 ; j<nx+1 ; ++j ){
            v_n[i][j] = v_tilda[i][j] - dt*( P_corr[i][j] - P_corr[i-1][j] )/dy ; 
        }
    }
    std::tie( u_n , v_n ) = Update_Boundaries( u_n , v_n );
}

void Velocity_L2Norm(){
    err_u =0 ; err_v =0;
    for( int i=1 ; i<ny+1 ; ++i ){
        for( int j=2 ; j<nx+1 ; ++j ){
            err_u += pow( (u_n[i][j] - u_nm1[i][j] ) , 2 )/((nx-1)*(ny-1));
        }
    }
    for( int i=2 ; i<ny+1 ; ++i ){
        for( int j=1 ; j<nx+1 ; ++j){
            err_v += pow( (v_n[i][j] - v_nm1[i][j] ) , 2 )/((nx-1)*(ny-1));      
        }
    }
    err_u = pow( err_u , 0.5 ) ; err_v = pow( err_v , 0.5 ); 
}



/* ============================================================ OUTPUT ===================================================================================================== */
void printmatrix( matrix mat ){
    int Rows = mat.size();
    int Cols = mat[0].size();
    for( int i=0 ; i<Rows ; ++i){
        for( int j=0 ; j<Cols ; ++j){
            cout<<mat[i][j]<<'\t';
            if( j==Cols-1 ) cout<<endl;
        }
    } 
}

void Output_Solution(){
    ofstream fout("u_128.txt");
    for( int i=0 ; i<ny+2 ; ++i  ){
        for( int j=0 ; j<nx+2 ; ++j){
            fout<<u_n[i][j]<<'\t';
            if( j==ny+1 ) fout<<endl;
        }
    }
    fout.close();

    fout.open("v_128.txt");
    for( int i=0 ; i<ny+2 ; ++i  ){
        for( int j=0 ; j<nx+2 ; ++j){
            fout<<v_n[i][j]<<'\t';
            if( j==ny+1 ) fout<<endl;
        }
    }
    fout.close();


}