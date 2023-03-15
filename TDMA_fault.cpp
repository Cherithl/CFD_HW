void TDMA( int row){
//FORWARD SUBSTITUTION to reduce to UPPER TRIANGULAR MATRIX
    //cout<<"TDMA"<<endl<<endl;
double D[nx-2]={};double W[nx-2]={}; double E[nx-2]={};
D[0]=4;     E[0]=-1;          W[nx-3]=-1;   D[nx-3]=4; 
    for(int i=1;i<nx-3;++i){
            E[i]=W[i]=-1; D[i]=4; 
    }
double Q[nx-2]={};
    for(int j=0;j<nx-2;++j){
        if(j==0)

        Q[j]=S[row][j+1]*dx*dy+T[row-1][row+1]+T[row+1][row+1]+T[row][j];
        else if(j==nx-3)
        Q[j]=S[row][j+1]*dx*dy+T[row-1][row+1]+T[row+1][j+1]+T[row][j+2];
        else{
            Q[j]=S[row][j+1]*dx*dy+T[row-1][row+1]+T[row+1][j+1];//VERIFY AGAIN
            //cout<<T[i-1][j+1]<<'\t';
            }
    if(row==1){
        cout<<"Hello"<<endl;
        for(int j=0;j<nx-2;++j){
            cout<<Q[j]<<'\t';
        } 
    cout<<endl;
    }
    }
    
    // cout<<"0"<<"\t"<<D[0]<<endl;    
    for(int j=1;j<nx-2;++j){
      D[j]=  D[j]  -  (W[j]/D[j-1])*E[j-1];
    //   cout<<j<<'\t'<<D[j]<<endl;
      Q[j]=Q[j]    -  (W[j]/D[j-1])*Q[j-1];
    }
    if(row==1){
        cout<<"fuck off"<<endl;
                for(int j=0;j<nx-2;++j){
                    cout<<Q[j]<<'\t';
                } 
                cout<<endl;
            }


//BACKWRARD SUBSTITUTION for SOLUTION
    //T[i][nx-2]=Q[nx-3]/D[nx-3];
    for(int j=nx-2;j>0;--j){
      T[row][j]= (Q[j]- E[j]*T[row][j+1])/D[j];
    }
}
