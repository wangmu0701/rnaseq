#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <cpmseq.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  N.allocate("N");
  n.allocate("n");
  y.allocate(1,N,1,n,"y");
  tran.allocate(1,N,1,8,"tran");
  start.allocate(1,N,"start");
  end.allocate(1,N,"end");
  Ntran.allocate(1,N,"Ntran");
  NT.allocate("NT");
  Unitran.allocate(1,NT,"Unitran");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  lambda0.allocate(0.001,1,"lambda0");
  lambda1.allocate(0.001,50,"lambda1");
  lambda2.allocate(0.001,50,"lambda2");
  lambda3.allocate(0.001,50,"lambda3");
  lambda4.allocate(0.001,50,"lambda4");
  lambda5.allocate(0.001,50,"lambda5");
  lambda6.allocate(0.001,50,"lambda6");
  lambda7.allocate(0.001,50,"lambda7");
  lambda8.allocate(0.001,50,"lambda8");
  lambda9.allocate(0.001,50,"lambda9");
  lambda10.allocate(0.001,50,"lambda10");
  lambda11.allocate(0.001,50,"lambda11");
  lambda12.allocate(0.001,50,"lambda12");
  lambda13.allocate(0.001,50,"lambda13");
  lambda14.allocate(0.001,50,"lambda14");
  lambda15.allocate(0.001,50,"lambda15");
  pi11.allocate(0.001,0.999,"pi11");
  pi21.allocate(0.001,0.999,"pi21");
  pi31.allocate(0.001,0.999,"pi31");
  pi41.allocate(0.001,0.999,"pi41");
  pi51.allocate(0.001,0.999,"pi51");
  pi61.allocate(0.001,0.999,"pi61");
  pi71.allocate(0.001,0.999,"pi71");
  pi81.allocate(0.001,0.999,"pi81");
  pi91.allocate(0.001,0.999,"pi91");
  pi101.allocate(0.001,0.999,"pi101");
  pi111.allocate(0.001,0.999,"pi111");
  pi121.allocate(0.001,0.999,"pi121");
  pi131.allocate(0.001,0.999,"pi131");
  pi141.allocate(0.001,0.999,"pi141");
  pi151.allocate(0.001,0.999,"pi151");
  lamb.allocate(1,15,1,2,"lamb");
  #ifndef NO_AD_INITIALIZE
    lamb.initialize();
  #endif
  pi.allocate(1,15,1,2,"pi");
  #ifndef NO_AD_INITIALIZE
    pi.initialize();
  #endif
  ftmp.allocate("ftmp");
  #ifndef NO_AD_INITIALIZE
  ftmp.initialize();
  #endif
  f.allocate("f");
}

void model_parameters::userfunction(void)
{
  lamb(1,1)=lambda0;
  lamb(2,1)=lambda0;
  lamb(3,1)=lambda0;
  lamb(4,1)=lambda0;
  lamb(5,1)=lambda0;
  lamb(6,1)=lambda0;
  lamb(7,1)=lambda0;
  lamb(8,1)=lambda0;
  lamb(9,1)=lambda0;
  lamb(10,1)=lambda0;
  lamb(11,1)=lambda0;
  lamb(12,1)=lambda0;
  lamb(13,1)=lambda0;
  lamb(14,1)=lambda0;
  lamb(15,1)=lambda0;
  lamb(1,2)=lambda1;
  lamb(2,2)=lambda2;
  lamb(3,2)=lambda3;
  lamb(4,2)=lambda4;
  lamb(5,2)=lambda5;
  lamb(6,2)=lambda6;
  lamb(7,2)=lambda7;
  lamb(8,2)=lambda8;
  lamb(9,2)=lambda9;
  lamb(10,2)=lambda10;
  lamb(11,2)=lambda11;
  lamb(12,2)=lambda12;
  lamb(13,2)=lambda13;
  lamb(14,2)=lambda14;
  lamb(15,2)=lambda15;
  pi(1,1)=pi11;
  pi(1,2)=1-pi11;
  pi(2,1)=pi21;
  pi(2,2)=1-pi21;
  pi(3,1)=pi31;
  pi(3,2)=1-pi31;
  pi(4,1)=pi41;
  pi(4,2)=1-pi41;
  pi(5,1)=pi51;
  pi(5,2)=1-pi51;
  pi(6,1)=pi61;
  pi(6,2)=1-pi61;
  pi(7,1)=pi71;
  pi(7,2)=1-pi71;
  pi(8,1)=pi81;
  pi(8,2)=1-pi81;
  pi(9,1)=pi91;
  pi(9,2)=1-pi91;
  pi(10,1)=pi101;
  pi(10,2)=1-pi101;
  pi(11,1)=pi111;
  pi(11,2)=1-pi111;
  pi(12,1)=pi121;
  pi(12,2)=1-pi121;
  pi(13,1)=pi131;
  pi(13,2)=1-pi131;
  pi(14,1)=pi141;
  pi(14,2)=1-pi141;
  pi(15,1)=pi151;
  pi(15,2)=1-pi151;
  f=0;
  for(int exon=1;exon<=N;exon++)
  {
    int a1=tran(exon,1);
    int a2=tran(exon,2);
    int a3=tran(exon,3);
    int a4=tran(exon,4);
    int a5=tran(exon,5);
    int a6=tran(exon,6);
    int a7=tran(exon,7);
    int a8=tran(exon,8);
    for(int m=start(exon);m<=end(exon);m++)
    {
      ftmp=0;
      if(Ntran(exon)==8)
      {
        for(int i1=1;i1<=2;i1++)
        {
        for(int i2=1;i2<=2;i2++)
        {
        for(int i3=1;i3<=2;i3++)
        {
        for(int i4=1;i4<=2;i4++)
        {
        for(int i5=1;i5<=2;i5++)
        {
        for(int i6=1;i6<=2;i6++)
        {
        for(int i7=1;i7<=2;i7++)
        {
        for(int i8=1;i8<=2;i8++)
        {
          ftmp=ftmp+pi(a1,i1)*pi(a2,i2)*pi(a3,i3)*pi(a4,i4)*pi(a5,i5)*pi(a6,i6)*pi(a7,i7)*pi(a8,i8)*exp(-lamb(a1,i1)-lamb(a2,i2)-lamb(a3,i3)-lamb(a4,i4)-lamb(a5,i5)-lamb(a6,i6)-lamb(a7,i7)-lamb(a8,i8))*pow(lamb(a1,i1)+lamb(a2,i2)+lamb(a3,i3)+lamb(a4,i4)+lamb(a5,i5)+lamb(a6,i6)+lamb(a7,i7)+lamb(a8,i8),y(exon,m));
        }
        }
        }
        }
        }
        }
        }
        }
        f=f-log(ftmp);
      }
    }
  }
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(void){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(32000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(32000000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
