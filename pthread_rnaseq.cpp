#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h> 
#include <pthread.h>
#include <errno.h>

#define Nexon   2
#define Ntran   3
#define Nlen    4000  
#define MAX_THREADS 512
#define EPSILON 5E-5
#define CHUNKSIZE 100  
using namespace std;

struct readThreadParam{
  int tid;
  int para;
};
double l01,l11,l21,ll;  //parameters to optimize

double l[Ntran][2];
double p[Ntran][3];

double oldl[Ntran][2];
double oldp[Ntran][3];
double eu[Nexon][Nexon][Nexon][Nlen];
double y[Nexon][Nlen]; //read count data
unsigned int elen[Nexon]; //length of exon
unsigned int estart[Nexon]; 
unsigned int estop[Nexon];
unsigned int trans[Nexon][2]; // indicates from which transcripts, the exon come
int numThreads;
int nItr=0;
double likelihood_value;
double propResolution =0.001;
double lambdaResolution =0.1;
double minfuncval=999999.0, minlambda,minprop;
pthread_mutex_t likelihood_mutex;
pthread_mutex_t parameter_mutex;

//parameters are p[0][0], p[1][0], p[2][0], l[0][1], l[1][1], l[2][1], ll = l[0][1] =l[1][1]=l[2][1], total 7 parameteres
double newparam[7], oldparam[7];//parameter list

double optimvalue;

double elapsedtime(timeval start, timeval end){
  return (end.tv_sec - start.tv_sec) * 1000.0+ (end.tv_usec - start.tv_usec) / 1000.0;
}
struct threadParameter{
  int tid;
  int paramterid;
};

//likelihood_func() the likelihood function, line 160 in solve.2.exon.R
//pp_func() the proportion function
//ff_func() the intensity funcion
//update_eu()

void printData(){
  int sum=0;
  int i,j;
    for(i=0;i<Nexon-1;i++){
         sum=0;
         for(j=0;j<Nlen;j++){
            sum+=y[i][j];
            printf("%d ",(int)y[i][j]);
        }
        printf("%d\n",sum);
    }
    //sum1 =7636, sum2=6388
}

//readdata from y.dat
//reading the read counts
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
    ll=0.05;// let l[0][1] =l[1][1]=l[2][1] =0.05 , assumption of the model
    l[0][0]=7.0;l[1][0]=8.0;l[2][0]=9.0;
    p[0][0]=0.4;p[0][1]=1.0-p[0][0];
    p[1][0]=0.3;p[1][1]=1.0-p[1][0];
    p[2][0]=0.2;p[2][1]=1.0-p[2][0];
    l[0][1] =l[1][1]=l[2][1] =ll;
    fclose(fp);
   // printData();
    

}

double pois_fun(double lambda, double k){
    return (exp(-lambda)*pow(lambda,k));
}

void* lambda_func(void* params ){  
    unsigned int i,j,k,e, t0, t1;
    struct threadParameter *readParams = ((struct threadParameter*) params);
    int tID = readParams->tid;
    int paraID = readParams->paramterid;
    //assume range of lambda parameter 0....10
    double funcValue=0;
    double numPoints = (30-0)/lambdaResolution;
    double param;
    double locall[Ntran][3];

    int len;
    int pointsPerThread = len = (int) ( numPoints/numThreads);
    if(tID == numThreads-1 && pointsPerThread  * numThreads < (int)numPoints)
      len+= ((int)numPoints % numThreads );
    //printf("tid:%d paraid:%d, Points:%d\n", tID, paraID,pointsPerThread );
    double startPoint = (tID*pointsPerThread)*lambdaResolution  ;//avoid 0 probability
    double endPoint = startPoint + len *lambdaResolution - lambdaResolution;
    //printf("tid:%d paraid:%d, Points:%d, startP:%.3f endP:%.3f\n", tID, paraID,pointsPerThread, startPoint, endPoint );
    for(i=0;i<Ntran;i++){//initialization
      locall[i][0]= l[i][0];
      locall[i][1]= l[i][1];
    }
     
    for( param = startPoint+1e-2; param<endPoint; param+= lambdaResolution){
      funcValue=0;
      if(paraID <=5 && paraID >=3)//l[0][0], l[1][0], l[2][0]
        locall[paraID-3][0] = param;//this changes
      else if(paraID ==6){// l[0][1]= l[1][1]=l[2][1] (model assumption)
        for(i=0;i<Ntran;i++)
	  locall[i][1]=param;
      }  
      else{
        printf("invalid parameter\n");
      }
      for(e=0; e<Nexon;e++){
        t0=trans[e][0];
        t1=trans[e][1];
        for(i=0;i<2;i++){
          for(j=0;j<2;j++){
             for(k=0;k<elen[e];k++){
	        //printf("locall[t0][i]:%.6f locall[t1][j]:%.6f, sum=%.3f\n", locall[t0][i],locall[t1][j], locall[t0][i]+locall[t1][j]);
		//printf("sum :%.3f\n", locall[t0][i]+locall[t1][j]);
                funcValue+=eu[e][i][j][k]*(-locall[t0][i]-locall[t1][j]+y[e][estart[e]+k]*log(locall[t0][i]+locall[t1][j]));
		//printf("eu:%.3f\n", eu[e][i][j][k]);
             }
          }
        }
      } 
      //printf("-funcValue:%.3f\n", -funcValue);
      if(-funcValue < minlambda){
        pthread_mutex_lock(&parameter_mutex); 
	//oldparam[paraID]=newparam[paraID];
        newparam[paraID] = param;
        minlambda = -funcValue;
        if(paraID>=3 && paraID <=5)
  	  l[paraID-3][0] = param;
        else{//paraID 7
	  for(i=0;i<Ntran;i++){
            l[i][1] = param;//all of them are same(model assumption)
	  }
	}  
	//printf("id:%d min param :%0.3f val:%.3f\n",paraID,param, minfuncval);
        
	pthread_mutex_unlock(&parameter_mutex);
      }       
    }//for*/
    return NULL;
}

void* proportion_func(void* params ){  
    unsigned int i,j,k,e, t0, t1;
    struct threadParameter *readParams = ((struct threadParameter*) params);
    int tID = readParams->tid;
    int paraID = readParams->paramterid;
    //range of proportion is parameter 0....1
    double funcValue=0;
    double numPoints = (1-0)/propResolution;
    double param;
    double localp[Ntran][3];

    int len;
    int pointsPerThread = len = (int) ( numPoints/numThreads);
    if(tID == numThreads-1 && pointsPerThread  * numThreads < (int)numPoints)
      len+= ((int)numPoints % numThreads );
    //printf("tid:%d paraid:%d, Points:%d\n", tID, paraID,pointsPerThread );
    double startPoint = (tID*pointsPerThread)*propResolution  ;//avoid 0 probability
    double endPoint = startPoint + len *propResolution - propResolution;
    //printf("tid:%d paraid:%d, Points:%d, startP:%.3f endP:%.3f\n", tID, paraID,pointsPerThread, startPoint, endPoint );
    for(i=0;i<Ntran;i++){
      localp[i][0]= p[i][0];
      localp[i][1]= p[i][1];
    }
     
    for( param = startPoint+1e-2; param<endPoint; param+= propResolution){
      funcValue=0;
      localp[paraID][0] = param;//this changes
      localp[paraID][1]=1-param;

      for(e=0; e<Nexon;e++){
        t0=trans[e][0];
        t1=trans[e][1];
        for(i=0;i<2;i++){
          for(j=0;j<2;j++){
             for(k=0;k<elen[e];k++){
                  funcValue+=eu[e][i][j][k]*log(localp[t0][i]*localp[t1][j]);
             }
          }
        }
      } 
      if(-funcValue < minprop){
        pthread_mutex_lock(&parameter_mutex); 
        //oldparam[paraID]=newparam[paraID];
	newparam[paraID] = param;
        minprop = -funcValue;
        p[paraID][0] = param;
        p[paraID][1] = 1 - param; //update for next iteration
	//printf("id:%d min param :%0.3f val:%.3f\n",paraID,param, minfuncval);
        pthread_mutex_unlock(&parameter_mutex);
      }       
    }//for
    return NULL;
}

double normParamDiff(){
  double sum=0;
  int i;
  for(i=0;i<7;i++)
    sum += (oldparam[i]- newparam[i])* (oldparam[i]- newparam[i]);
  return sqrt(sum);  
}

void saveParam(){
  int i;
  for(i=0;i<7; i++)
    oldparam[i] = newparam[i];
}
void printParam(){
  int i;
  for(i=0;i<7;i++)
    printf("new:%.3f, ", newparam[i]);
  printf("\n");  
}

double likelihood_ind(unsigned int e, int tid){
    unsigned int i,j,k;
    unsigned int t0=trans[e][0];
    unsigned int t1=trans[e][1];
    double ret=0;
    int nodePerThread = elen[e]/numThreads, len;
    int start = tid * nodePerThread ;//start 
    len = nodePerThread;// length of chunk to work on
    if( tid == numThreads-1)//last thread
      len += (elen[e]%numThreads);

    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=start;k<start+len;k++){
                ret+=eu[e][i][j][k]*( log(p[t0][i]*p[t1][j])-l[t0][i]-l[t1][j]+y[e][estart[e]+k]*log(l[t0][i]+l[t1][j])  );//equation 7,8
            }
        }
    }
    return ret;
}

void* likelihood_func(void *x){
    double ret=0;
    unsigned int i,e;
     int tid =*((int *) x);
    for(e=0;e<Nexon;e++){
        ret+=likelihood_ind(e,tid );
    }
    pthread_mutex_lock(&likelihood_mutex);
    likelihood_value+=ret;
    pthread_mutex_unlock(&likelihood_mutex);
}

void update_eu_ind(unsigned int e, int tid){ //likelihood calculation
    unsigned int t0=trans[e][0];
    unsigned int t1=trans[e][1];
    int nodePerThread = elen[e]/numThreads, len;
    int start = tid * nodePerThread ;//start 
    len = nodePerThread;// length of chunk to work on
    if( tid == numThreads-1)//last thread
      len += (elen[e]%numThreads);
    double *sumtmp=(double*)malloc(sizeof(double)*elen[e]);
    unsigned int i,j,k;
    //printf("tid:%d , start:%d end:%d\n", tid, start, start+len);
    //Initial to 0
    for(i=0;i<elen[e];i++){
        sumtmp[i]=0;
    }
//Compute sum temporary
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
          for(k=start;k<start+len;k++){
                sumtmp[k]+=pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];//P1iP3j Poi(lamba1i +lambda3j) for each base pair position
            }
        }
    }
//Store averages
    double tmp=0;
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            for(k=start;k<start+len;k++){
                tmp=pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];
                eu[e][i][j][k]=tmp/sumtmp[k];//normalization
            }
        }
    }
    free(sumtmp);
}

//update.Eu
void* update_eu(void *x){
    unsigned int i,e,j;
    int tid =*((int *) x);
    for(e=0;e<Nexon;e++)
      update_eu_ind(e,tid );
    
}

void callE(){
  int *taskIDPtr, i,j;
  pthread_t p_threads[MAX_THREADS];

  for( i=0;i<numThreads;i++){
      taskIDPtr = (int *) malloc(sizeof(int));
      *taskIDPtr = i; 
      pthread_create(&p_threads[i], NULL, update_eu, (void*)taskIDPtr);
  }
  for( i=0;i<numThreads;i++){
      pthread_join(p_threads[i], NULL);//wait for threads to terminate
  }
}

void callM(){
  int *taskIDPtr, i,j;
  pthread_t p_threads[MAX_THREADS];
  struct threadParameter* threadParam;
  
  for(i=0;i<3;i++){//proportion parameters
    for( j=0;j<numThreads;j++){
      threadParam = (struct threadParameter*  )malloc (sizeof(struct threadParameter )); 
      threadParam->paramterid = i;
      threadParam->tid =j;
      pthread_create(&p_threads[j], NULL, proportion_func, threadParam);
    }
    for( j=0;j<numThreads;j++){
      pthread_join(p_threads[j], NULL);//wait for threads to terminate
    }
    //reset for nex parameter 
  }//proportion params
  //printf("p[0][0]:%.3f, p[1][0]:%.3f, p[2][0]:%.3f\n", p[0][0], p[1][0], p[2][0]);
  
  //minfuncval=999999.0;
  for(i=3;i<7;i++){//lambda parameters
    for( j=0;j<numThreads;j++){
      threadParam = (struct threadParameter*  )malloc (sizeof(struct threadParameter )); 
      threadParam->paramterid = i;
      threadParam->tid =j;
      pthread_create(&p_threads[j], NULL, lambda_func, threadParam);
    }
    for( j=0;j<numThreads;j++){
      pthread_join(p_threads[j], NULL);//wait for threads to terminate
    }
    //minfuncval=999999.0;//reset for nex parameter 
    
 }

}

void funcEM(){
  int *taskIDPtr, i,j;
  pthread_t p_threads[MAX_THREADS];
  double li;
  double diff=9999.0;
  minlambda=999999.0;
  minprop=minlambda;
  while( diff > EPSILON){  
    saveParam();//save parameters
    callE();//E step (parallelized)
    callM(); //M step (parallelized)   
    nItr++;
    diff= normParamDiff();
    //printf("Iteration:%d %.3f\n", nItr, diff);
    //printParam();
      
  }
  //printf("likelihood=%15.10e\n",li);
}

//upto line 170 of solve.two.exon.R
//The any.exon() function, which solve the problem is not implemented
//That's what we need to do. 
int main(){
    timeval t1, t2;
    double totalTime,li;
    int i=0;
    pthread_mutex_init(&(parameter_mutex), NULL);
    readData();
    printf("Enter number of threads: ");
    scanf("%d", &numThreads);
    printf("Starting iteration\n");
    gettimeofday(&t1, NULL);
    funcEM();
    gettimeofday(&t2, NULL);
    printf("Ending iteration, TimeTaken: %0.5f\n", elapsedtime(t1,t2)/1000.0);
    printParam();
    printf("Iterations: %d minimum: %.3f\n",nItr, minlambda+ minprop);
    return 0;
}
