#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Nexon   2
#define Ntran   3
#define Nlen    4000    
double ll;
double l[Ntran][2];
double p[Ntran][3];
double eu[Nexon][Nexon][Nexon][Nlen];
double y[Nexon][Nlen];
unsigned int elen[Nexon];
unsigned int estart[Nexon];
unsigned int estop[Nexon];
unsigned int trans[Nexon][2];

void readData(){
    FILE *fp=fopen("y.dat","r");
    char buf[20000];
    unsigned int i,j,k;
    for(i=0;i<Nexon;i++){
        fgets(buf,20000,fp);
        k=0;
        for(j=0;j<Nlen;j++){
            y[i][j]=atoi(&buf[k]);
            while((buf[k]!=' ') && (buf[k]!='\0')){
                k++;
            }
            if (buf[k]!='\0'){
                k++;
            }
        }
    }
    trans[0][0]=0;trans[0][1]=2;
    trans[1][0]=1;trans[1][1]=2;
    elen[0]=2000;estart[0]=0;estop[0]=2000;
    elen[1]=2000;estart[1]=2000;estop[1]=4000;

//Set a good initial value for start
    ll=0.05;
    l[0][0]=7.0;l[1][0]=8.0;l[2][0]=9.0;
    p[0][0]=0.4;p[0][1]=1.0-p[0][0];
    p[1][0]=0.3;p[1][1]=1.0-p[1][0];
    p[2][0]=0.2;p[2][1]=1.0-p[2][0];
    fclose(fp);
/*
    for(i=0;i<Nexon;i++){
        for(j=0;j<Nlen;j++){
            printf("%d ",(int)y[i][j]);
        }
        printf("\n");
    }
*/
}

double pois_fun(double l, double y){
    return (exp(-l)*pow(l,y));
}

double f_ind(unsigned int e){
    unsigned int i,j,k;
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
    unsigned int i,e;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        ret+=f_ind(e);
    }
  return ret;
}

double p_ind(unsigned int e){
    unsigned int i,j,k;
    unsigned int t0=trans[e][0];
    unsigned int t1=trans[e][1];
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
    unsigned int i,e;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        ret+=p_ind(e);
    }
    return ret;
}

double likelihood_ind(unsigned int e){
    unsigned int i,j,k;
    unsigned int t0=trans[e][0];
    unsigned int t1=trans[e][1];
    double ret=0;
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=0;k<elen[e];k++){
                ret+=eu[e][i][j][k]*( log(p[t0][i]*p[t1][j])-l[t0][i]-l[t1][j]+y[e][estart[e]+k]*log(l[t0][i]+l[t1][j])  );
            }
        }
    }
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=0;k<10;k++){
                printf("li[%d][%d][%d][%d]=%15.5f\n",e,i,j,k, eu[e][i][j][k]*( log(p[t0][i]*p[t1][j])-l[t0][i]-l[t1][j]+y[e][estart[e]+k]*log(l[t0][i]+l[t1][j])));
            }
        }
    }
    printf("li[%d]=%15.5f\n",e,ret);
    return ret;
}
double likelihood_func(){
    double ret=0;
    unsigned int i,e;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        ret+=likelihood_ind(e);
    }
    return ret;
}

void update_eu_ind(unsigned int e){
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
                tmp=pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];
                eu[e][i][j][k]=tmp/sumtmp[k];
            }
        }
    }

    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=0;k<10;k++){
                printf("eu[%d][%d][%d][%d]=%15.5f\n",e,i,j,k,eu[e][i][j][k]);
//                printf("y=%15.5f, l=%15.5f, posi=%15.5f\n",l[t0][i]+l[t1][j],y[e][estart[e]+k],pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j]);
            }
        }
    }
//exit(-1);
    free(sumtmp);
}
void update_eu(){
    unsigned int i,e;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        update_eu_ind(e);
    }
}

int main(){
    readData();
    update_eu();
    double li=likelihood_func();
    printf("likelihood=%15.5f\n",li);
    return 0;
}

