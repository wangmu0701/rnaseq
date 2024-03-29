#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <adolc/adolc.h>
#include <adolc/sparse/sparsedrivers.h>
#include <adolc/hessian/edge_main.h>
#include <sys/time.h>

#define Nexon   2
#define Ntran   3
#define Nlen    4000    

#define tag     1

#define pois_fun(l,y)   (exp(-(l))*pow((l),(y)))


adouble ll;
adouble l[Ntran][2];
adouble p[Ntran][2];
adouble eu[Nexon][Nexon][Nexon][Nlen];
adouble sumtmp[Nlen];

double ll0;
double l0[Ntran][2];
double p0[Ntran][2];
double x[7];

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
    ll0=0.05;
    l0[0][0]=7.0;l0[1][0]=8.0;l0[2][0]=9.0;
    p0[0][0]=0.4;
    p0[1][0]=0.3;
    p0[2][0]=0.2;
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

//adouble pois_fun(adouble l, adouble y){
//    return (exp(-l)*pow(l,y));
//}



adouble ff_func(){
    adouble ret=0;
    unsigned int i,e,j,k;
    unsigned int t0;
    unsigned int t1;
    adouble tmp=0;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        t0=trans[e][0];
        t1=trans[e][1];
        tmp=0;
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                for(k=0;k<elen[e];k++){
                tmp+=eu[e][i][j][k]*(-l[t0][i]-l[t1][j]+y[e][estart[e]+k]*log(l[t0][i]+l[t1][j]));
                }
            }
        }
        ret-=tmp;
    printf("ff[%d]=%15.10e\n",e,ret.getValue());
    }
  return ret;
}


adouble pp_func(){
    adouble ret=0;
    unsigned int i,j,k,e;
    unsigned int t0;
    unsigned int t1;
    adouble tmp=0;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        t0=trans[e][0];
        t1=trans[e][1];
        tmp=0;
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                for(k=0;k<elen[e];k++){
                    tmp+=eu[e][i][j][k]*log(p[t0][i]*p[t1][j]);
                }
            }
        }
        ret-=tmp;
    }
    return ret;
}


adouble likelihood_func(){
    adouble ret=0;
    unsigned int i,j,k,e;
    unsigned int t0;
    unsigned int t1;
    adouble tmp=0;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        t0=trans[e][0];
        t1=trans[e][1];
        tmp=0;
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                for(k=0;k<elen[e];k++){
                    tmp+=eu[e][i][j][k]*( log(p[t0][i]*p[t1][j])-l[t0][i]-l[t1][j]+y[e][estart[e]+k]*log(l[t0][i]+l[t1][j])  );
                }
            }
        }
        ret+=tmp;
    }
    return ret;
}

void update_eu(){
    unsigned int i,j,k,e;
    unsigned int t0;
    unsigned int t1;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        t0=trans[e][0];
        t1=trans[e][1];
        
        //Initial to 0
        for(i=0;i<elen[e];i++){
            sumtmp[i]=0;
        }
        //Compute sum temporary
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                for(k=0;k<elen[e];k++){
                    sumtmp[k]+=pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];
                    //                printf("sum[%d]=%15.5f\n",k,sumtmp[k].getValue());
                }
            }
        }
        //    exit(-1);
        //Store averages
        adouble tmp=0;
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                for(k=0;k<elen[e];k++){
                    tmp=pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];
                    eu[e][i][j][k]=tmp/sumtmp[k];
                }
            }
        }
        
    }
}

int main(){
    struct timeval tv1, tv2;
    double time_elapsed = 0;
    readData();

    trace_on(tag);
    int i,j;
    int n=7;

    ll<<=ll0;
    x[0]=ll0;
    for(i=0;i<3;i++){
        l[i][0]<<=l0[i][0];
        x[1+i]=l0[i][0];
    }
    for(i=0;i<3;i++){
        p[i][0]<<=p0[i][0];
        x[4+i]=p0[i][0];
    }
    for(i=0;i<3;i++){
        p[i][1]=1.0-p[i][0];
    }
    printf("ll=%15.5f\n",ll.getValue());
    printf("l[0][0]=%15.5f\n",l[0][0].getValue());
    printf("l[1][0]=%15.5f\n",l[1][0].getValue());
    printf("l[2][0]=%15.5f\n",l[2][0].getValue());
    printf("p[0][0]=%15.5f\n",p[0][0].getValue());
    printf("p[1][0]=%15.5f\n",p[1][0].getValue());
    printf("p[2][0]=%15.5f\n",p[2][0].getValue());

    gettimeofday(&tv1, NULL);
    update_eu();
    adouble li=likelihood_func();
    gettimeofday(&tv2, NULL);
    time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                   (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
    double li0;
    li >>= li0;

//    adouble ff=ff_func();
//    double ff0;
//    ff >>= ff0;

//    adouble pp=pp_func();
//    double pp0;
//    pp >>= pp0;

    trace_off();

    printf("likelihood=%.10e, time elapsed = %.10f",li0, time_elapsed);
//    printf("ff_func=%15.10e\n",ff0);
//    printf("pp_func=%15.10e\n",pp0);

    double **H;
    H=myalloc2(n,n);
    gettimeofday(&tv1,NULL);
    hessian(tag,n,x,H);
    gettimeofday(&tv2,NULL);
    printf("Computing the full hessian cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("H[%d][%d]=%10.5f",i,j,H[i][j]);
        }
        printf("\n");
    }
    myfree2(H);
    unsigned int *rind=NULL;
    unsigned int *cind=NULL;
    double *values=NULL;
    int nnz;
    int options[2]={0,1};
    gettimeofday(&tv1,NULL);
    edge_hess(tag,1,n,x,&nnz,&rind,&cind,&values,options);
    gettimeofday(&tv2,NULL);
    printf("Sparse Hessian: edge pushing cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
    for(i=0;i<nnz;i++){
        printf("EH[%i][%i]=%10.5f\n",rind[i],cind[i],values[i]);
    }

    free(rind);rind=NULL;
    free(cind);cind=NULL;
    free(values);values=NULL;
    return 0;
}

