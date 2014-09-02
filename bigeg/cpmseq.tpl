TOP_OF_MAIN_SECTION
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(32000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(32000000);

// this is to fit the mixture of poisson model with three exons

DATA_SECTION
  init_int N  // number of exons
  init_int n  // number of observation
  init_matrix y(1,N,1,n) // y matrix with N types of reads, and n observations for each type 
  init_matrix tran(1,N,1,8)  // transcripts id that overlap each type of read
  init_vector start(1,N)    // starting positions 
  init_vector end(1,N)    //ending positions
  init_vector Ntran(1,N)  // number of transcripts overlapping in each exon
  init_int NT // number of transcripts to be solved totally
  init_vector Unitran(1,NT)  // unique transcript id that need to be solved
  // end of DATA_SECTION



PARAMETER_SECTION
  init_bounded_number lambda0(0.001,1)
  init_bounded_number lambda1(0.001,50)
  init_bounded_number lambda2(0.001,50)
  init_bounded_number lambda3(0.001,50)
  init_bounded_number lambda4(0.001,50)
  init_bounded_number lambda5(0.001,50)
  init_bounded_number lambda6(0.001,50)
  init_bounded_number lambda7(0.001,50)
  init_bounded_number lambda8(0.001,50)
  init_bounded_number lambda9(0.001,50)
  init_bounded_number lambda10(0.001,50)
  init_bounded_number lambda11(0.001,50)
  init_bounded_number lambda12(0.001,50)
  init_bounded_number lambda13(0.001,50)
  init_bounded_number lambda14(0.001,50)
  init_bounded_number lambda15(0.001,50)
  init_bounded_number pi11(0.001,0.999)
  init_bounded_number pi21(0.001,0.999)
  init_bounded_number pi31(0.001,0.999)
  init_bounded_number pi41(0.001,0.999)
  init_bounded_number pi51(0.001,0.999)
  init_bounded_number pi61(0.001,0.999)
  init_bounded_number pi71(0.001,0.999)
  init_bounded_number pi81(0.001,0.999)
  init_bounded_number pi91(0.001,0.999)
  init_bounded_number pi101(0.001,0.999)
  init_bounded_number pi111(0.001,0.999)
  init_bounded_number pi121(0.001,0.999)
  init_bounded_number pi131(0.001,0.999)
  init_bounded_number pi141(0.001,0.999)
  init_bounded_number pi151(0.001,0.999)
  matrix lamb(1,15,1,2)
  matrix pi(1,15,1,2)
  number ftmp
  objective_function_value f 
 
  
  
PROCEDURE_SECTION
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
