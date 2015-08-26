#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <time.h>       /*ʹ�õ�ǰʱ����Ϊ���������*/
#include <math.h>
#include <iostream>
#include <omp.h>
using namespace std;

#define MC_Time 20000
#define MAT 100
#define Label_NUM 10000  //C6���ӱ�Ƿ�

#define CALC_AA(Empty_NUM,C2_CON) C2_CON*10*0.16*0.443*Empty_NUM/MAT_SIZE
//#define CALC_BA(x,y) y*509.36*0.16*pow((0.443*x/MAT_SIZE),3)
#define CALC_BA(Empty_NUM,C6_CON) C6_CON*509.36*0.16*(0.443*Empty_NUM/MAT_SIZE)*(0.443*Empty_NUM/MAT_SIZE)*(0.443*Empty_NUM/MAT_SIZE)
#define CALC_AG(C2_NUM) 0.05*0.443*C2_NUM/MAT_SIZE
#define CALC_BG(C6_NUM) 0.2*0.443*C6_NUM/(MAT_SIZE*3)
#define CALC_AD(C2_NUM) 0.0386*0.443*C2_NUM/MAT_SIZE;
#define CALC_BD(C6_NUM) 0.0657*0.443*C6_NUM/(MAT_SIZE*3)
#define CALC_AB(C2_NUM,C6_NUM) 7*0.443*0.443*C2_NUM*C6_NUM/(MAT_SIZE*MAT_SIZE*3)

double MAT_SIZE = MAT*MAT;
int result[400][100][100];

/*C6�Ѹ����������Χ15��λ��*/
void C6_Desorp(int x,int y,int cur,int (&a)[MAT][MAT])
{
	int x1 = x-2, x2 = x+2;
	int y1 = y-2, y2 = y+2;
	for(int i = x1; i <= x2; i++)
		for(int j = y1; j<= y2; j++)
		{
			if(i>=0 && i<MAT && j>=0 && j<MAT)
			{
				if(a[i][j] == cur)
					a[i][j] = 0;
			}
		}
}

void reaction(double C2_CON,double C6_CON,int t)
{  /*�������*/
	unsigned seed;
	int x1,y1,x2,y2,x3,y3,x4,y4,coord[8],pc,label[Label_NUM] = {0};
	int i,j,n,nn,x,y,s,k,w,cur,tag;
	float aa,ba,ag,bg,ad,bd,ab,all;            /*ϵͳ�ڸ���Ӧ����*/
	float paa,pba,pag,pbg,pbd,pad,pab;      /*ϵͳ�ڸ���Ӧ����*/
	float R1,R3;int R2;                     /*���ɾ��������*/
	int surf[MAT][MAT]={0};                   /*��ʼ���������*/
	srand(time(NULL));                      /*��ʼ�������*/
	clock_t start,finish;
	/*��ѭ�� MC_Time �Σ�PP ��ֵΪ���ؿ���ʱ�� MCS*/
	double Empty_NUM=MAT_SIZE,C2_NUM=0,C6_NUM=0;
	seed = 25234 + 17*omp_get_thread_num();
	for(n=0;n<MC_Time;n++)
	{
		start = clock();
		/*һ�� MCS������ MAT_SIZE ��ѭ��*/
		for(nn=0;nn<MAT_SIZE;nn++)
		{
			pc = 0;
			s=0,k=0,tag=1;
			/*���� label ������δʹ�õ� tag*/
			for(i=1;i<Label_NUM;i++)
			{
				if(label[i]==0) {tag = i;break;}
			}

			/*��������ؼ���ʽ*/
			aa=CALC_AA(Empty_NUM,C2_CON);    /* C2 ��������*/
			ba=CALC_BA(Empty_NUM,C6_CON);    /* C6 ��������*/
			ag=CALC_AG(C2_NUM);              /* C2 �Ѹ�����*/
			bg=CALC_BG(C6_NUM);              /* C6 �Ѹ�����*/
			ad=CALC_AD(C2_NUM);              /* C2 ��������*/
			bd=CALC_BD(C6_NUM);              /* C6 ��������*/
			ab=CALC_AB(C2_NUM,C6_NUM);       /*˫���ӷ�Ӧ����*/
			all=aa+ba+ag+bg+ad+bd+ab;        /*������*/
			/*��Ӧ��Ӧ����*/
			paa=aa/all;pba=ba/all;pag=ag/all;pbg=bg/all;pad=ad/all;pbd=bd/all;pab=ab/all;

			/*����[0,1]���������*/
			R1=(float)rand_r(&seed)/RAND_MAX;
			R3=(float)rand_r(&seed)/RAND_MAX;

			/*�������Ӧ�����㷨��Ϊ�����߲���ʵ��*/
			/*1����������С������ C2,��-1 ��ʾ*/
			if(R1<=paa)
			{
				if(Empty_NUM<1) continue;
				do{
					R2=rand_r(&seed)%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=0);  /*�������һ�����λ��x,y��*/
				surf[x][y]=-1;    /*�ڴ˿�λ������ C2��C2 �������� 1 */
				C2_NUM += 1;
				Empty_NUM -= 1;
			}

			/*2������С���ӷ����� C6������������ʾ*/
			else if(R1<=paa+pba)
			{
				if(Empty_NUM<3) continue;
				do{
					R2=rand_r(&seed)%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=0);    /*�������һ�����λ*/

				x1=x-1,y1=y; //��
				x2=x,y2=y-1; //��
				x3=x+1,y3=y; //��
				x4=x,y4=y+1; //��
				if(x1>=0&&x1<MAT&&y1>=0&&y1<MAT&&surf[x1][y1]==0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<MAT&&y2>=0&&y2<MAT&&surf[x2][y2]==0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<MAT&&y3>=0&&y3<MAT&&surf[x3][y3]==0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<MAT&&y4>=0&&y4<MAT&&surf[x4][y4]==0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==4)
				{
					surf[coord[0]][coord[1]]=tag;
					surf[coord[2]][coord[3]]=tag;
				}
				if(pc==6)
				{
					/*�ڴ�����λλ�����ѡ������λ��˵��γ�����λ���� C6��ͬ��������C6 �� tag ��ʾ��ͬʱ C6 �������� 1*/
					if(R3<=(float)1/3){surf[coord[0]][coord[1]]=tag;surf[coord[2]][coord[3]]=tag;}
					else if(R3<=(float)2/3){surf[coord[2]][coord[3]]=tag;surf[coord[4]][coord[5]]=tag;}
					else{surf[coord[0]][coord[1]]=tag;surf[coord[4]][coord[5]]=tag;}
				}
				if(pc==8)
				{
					/*�ڴ�����λλ�����ѡ������λ��˵��γ�����λ���� C6��ͬ��������C6 �� tag ��ʾ��ͬʱ C6 �������� 1*/
					if(R3<=(float)1/4){surf[coord[0]][coord[1]]=tag;surf[coord[2]][coord[3]]=tag;}
					else if(R3<=(float)2/4){surf[coord[2]][coord[3]]=tag;surf[coord[4]][coord[5]]=tag;}
					else if(R3<=(float)3/4){surf[coord[4]][coord[5]]=tag;surf[coord[6]][coord[7]]=tag;}
					else{surf[coord[0]][coord[1]]=tag;surf[coord[6]][coord[7]]=tag;}
				}
				if(pc>=4)
				{	
					surf[x][y]=tag;
					C6_NUM += 3;
					Empty_NUM -= 3;
					label[tag] = 1;
				}
			}

			/*3������С������ C2 �Ѹ�*/
			else if(R1<=paa+pba+pag)
			{
				if(C2_NUM<1) continue;
				do{
					R2=rand_r(&seed)%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);  /*���ѡ��һ C2 ����λ*/
				surf[x][y]=0; /* C2 �Ѹ���ͬʱ C2 �Ѹ����� 1*/
				C2_NUM -= 1;
				Empty_NUM += 1;
			}

			/*4������С������ C2 ����*/
			else if(R1<=paa+pba+pag+pad)
			{
				if(C2_NUM<1) continue;
				do{
					R2=rand_r(&seed)%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);    /*���ѡ��һ C2 ����λ*/
				surf[x][y]=0;  /* C2 �Ѹ���ͬʱ C2 �������� 1*/
				C2_NUM -= 1;
				Empty_NUM += 1;
			}

			/*5��С���ӷ����� C6 �Ѹ�*/
			else if(R1<=paa+pba+pag+pad+pbg)
			{
				if(C6_NUM<3) continue;
				do{
					R2=rand_r(&seed)%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]<1);         /*���ѡ��һ C6 ����λ*/
				cur=surf[x][y];

				C6_Desorp(x,y,cur,surf);
				
				C6_NUM -= 3;
				Empty_NUM += 3;
				label[cur] = 0;
			}

			/*6��С���ӷ����� C6 ����*/
			else if(R1<=paa+pba+pag+pad+pbg+pbd)
			{
				if(C6_NUM<3) continue;
				do{
					R2=rand_r(&seed)%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]<1);      /*���ѡ��һ C6 ����λ*/
				cur=surf[x][y];

				C6_Desorp(x,y,cur,surf);
				
				C6_NUM -= 3;
				Empty_NUM += 3;
				label[cur] = 0;
			}

			/*7������С������ C2 ����ΧС���ӷ����� C6 ��Ӧ*/
			else
			{
				if(C6_NUM<3||C2_NUM<1) continue;
				do{
					R2=rand_r(&seed)%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);    /*���ѡ��һ C2 ����λ*/

				x1=x-1,y1=y;
				x2=x,y2=y-1;
				x3=x+1,y3=y;
				x4=x,y4=y+1;
				if(x1>=0&&x1<MAT&&y1>=0&&y1<MAT&&surf[x1][y1]>0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<MAT&&y2>=0&&y2<MAT&&surf[x2][y2]>0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<MAT&&y3>=0&&y3<MAT&&surf[x3][y3]>0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<MAT&&y4>=0&&y4<MAT&&surf[x4][y4]>0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==2)
				{/*�������һ�� C6 ����λ�������䷢��˫���ӷ�Ӧ*/
					cur=surf[coord[0]][coord[1]];
					C6_Desorp(coord[0],coord[1],cur,surf);
				}
				if(pc==4)
				{/*�������λ���� C6 ����*/
					/*���ѡ������һ��λ*/
					if(R3<0.5) 
					{
						cur=surf[coord[0]][coord[1]];
						C6_Desorp(coord[0],coord[1],cur,surf);
					}
					else 
					{
						cur=surf[coord[2]][coord[3]];
						C6_Desorp(coord[2],coord[3],cur,surf);
					}
				}
				if(pc==6)
				{/*�������λ���� C6 ����*/
					/*�ڴ�����λλ�����ѡ��һ��λ*/
					if(R3<=(float)1/3)
					{ 
						cur=surf[coord[0]][coord[1]];
						C6_Desorp(coord[0],coord[1],cur,surf);
					}
					else if(R3<=(float)2/3)
					{
						cur=surf[coord[2]][coord[3]];
						C6_Desorp(coord[2],coord[3],cur,surf);
					}
					else 
					{
						cur=surf[coord[4]][coord[5]];
						C6_Desorp(coord[4],coord[5],cur,surf);
					}
				}
				if(pc==8)
				{/*�������λ���� C6 ����*/
					/*�ڴ�����λλ�����ѡ��һ��λ*/
					if(R3<=(float)1/4)
					{ 
						cur=surf[coord[0]][coord[1]];
						C6_Desorp(coord[0],coord[1],cur,surf);
					}
					else if(R3<=(float)2/4)
					{
						cur=surf[coord[2]][coord[3]];
						C6_Desorp(coord[2],coord[3],cur,surf);
					}
					else if(R3<=(float)3/4)
					{
						cur=surf[coord[4]][coord[5]];
						C6_Desorp(coord[4],coord[5],cur,surf);
					}
					else 
					{
						cur=surf[coord[6]][coord[7]];
						C6_Desorp(coord[6],coord[7],cur,surf);
					}
				}
				if(pc >= 2)
				{
					surf[x][y]=0;    /*���� C2 ����λת��Ϊ��λ*/
					C6_NUM -= 3;
					C2_NUM -= 1;
					Empty_NUM += 4;
					label[cur] = 0;
				}
			}
		}
		//finish = clock();
		//cout << (double)(finish -start)/CLOCKS_PER_SEC << endl;
		//cout << n << endl;
	}

	/*for(i=0; i<MAT; i++)
		for(j=0; j<MAT; j++)
			result[t][i][j] = surf[i][j];*/

	//cout << "C2_CON="<<C2_CON<<" C6_CON="<<C6_CON<<endl;
        //cout << "C2: " << C2_NUM << "C6: " << C6_NUM << "Empty: " << Empty_NUM << endl;


	//printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",C2_NUM,C6_NUM,Empty_NUM);
}


int main()
{
	
	time_t start,finish;
	memset(result, 0, sizeof(result)); //��ʼ������������ά����
    /*����400��a,b���*/
    double C2_CON[5] = {0.1,0.2,0.3,0.4,0.5};
    double Rate[80];
    double C2C6[400][2];
    int t=0;
    for(int i=0; i<80; i++)
    {
    	Rate[i] = 0.2*i;
    }
    /*��ʼ���㷴Ӧ*/
    //omp_set_nested(1);
    //#pragma omp parallel for
    //freopen("result-powopenmp.out","w",stdout);
/*for(int i=0; i<5; i++)
    {
            for(int j=0; j<80; j++)
        {
                C2C6[t][0]=C2_CON[i];
		C2C6[t][1]=C2_CON[i]*Rate[j];
		t++;
        }
    }*/
	/*for(int i=0;i<400;++i)
		cout << "C2: " << C2C6[i][0] << "C6: " <<  C2C6[i][1] << endl;*/
	/*start = time(NULL);
#pragma omp parallel for
	for(int i=0;i<400;++i)
		reaction(C2C6[i][0],C2C6[i][1],i);
	finish = time(NULL);
	cout << "time: " << finish-start << endl;*/
    	
   start=time(NULL);


    for(int i=0; i<5; i++)
    {
    #pragma omp parallel for
    	    for(int j=0; j<80; j++)
    	{
        	reaction(C2_CON[i],C2_CON[i]*Rate[j],i*j);
    	}
	}

	finish=time(NULL);
	cout << "time: " << finish-start;

    return 0;
}
