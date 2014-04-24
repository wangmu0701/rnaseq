#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Nexon 2
#define Ntran 3

double ll;
double l[Ntran][2];
double p[Ntran][3];
double* eu[Nexon][Nexon][Nexon];
double* y[Nexon];
unsigned int elen[Nexon];
unsigned int estart[Nexon];
unsigned int estop[Nexon];
unsigned int trans[Nexon][2];

double pois_fun(double l, double y){
    return (exp((-l)*pow(l,y)));
}

double f_ind(unsigned int e){
    unsigned int i,j;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    unsigned int t0=trans[e][0];
    unsigned int t1=trans[e][1];
    double ret=0;
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=0;k<elen[e];k++){
                ret+=eu[e][i][j][k]*(-l[t0][i]-l[t1][i]+y[e][estart[e]+k]*log(l[t0][i])+l[t1][i]);
            }
        }
    }
    return (-ret);
}

double f_func(){
  double ret=0;
  unsigned int e;
  for(e=0;e<Nexon;e++){
      ret+=f_ind(e);
  }
  return ret;
}

double p_ind(unsigned int e){
    unsigned int i,j;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    unsigned int t0=trans[e][0];
    unsigned int t1=trans[e][0];
    double ret=0;
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=0;k<elen[e];k++){
                ret+=eu[e][i][j][k]*log(p[t0][i]+p[t1][j]);
            }
        }
    }
    return (-ret);
}

double p_func(){
    double ret=0;
    unsigned int e;
    for(e=0;e<Nexon;e++){
        ret+=p_ind(e);
    }
    return ret;
}

double likelihood_ind(unsigned int e){
    unsigned int i,j;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    unsigned int t0=trans[e][0];
    unsigned int t1=trans[e][0];
    double ret=0;
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=0;k<elen[e];k++){
                ret+=eu[e][i][j][k]*log(p[t0][i]+p[t1][j]);
            }
        }
    }
}
double likelihood_fund(){
    double ret=0;
    unsigned int e;
    for(e=0;e<Nexon;e++){
        ret+=likelihood_ind(e);
    }
    return ret;
}

void update_eu(unsigned int e){
    unsigned int t0=trans[e][0];
    unsigned int t1=trans[e][1];
    double *sumtmp=(double*)malloc(sizeof(double)*elen[e]);
    unsigned int i,j,k;

//Initial to 0
    for(i=0;i<elen[e];i++){
        sumtmp[i]=0;
    }
//Compute sum temporary
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=0;k<elen[e];k++){
                sumtmp[k]+=pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];
            }
        }
    }
//Store averages
    double tmp=0;
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=0;k<elen[e];k++){
                tmp=p_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];
                eu[e][i][j][k]=tmp/sumtmp[k];
            }
        }
    }
}

int main(){

}

