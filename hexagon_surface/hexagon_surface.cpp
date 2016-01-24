/*����Ϊ��� R=C6_CON/C2_CON��C2_CON,C6_CON��ʱ���������� C2 �� C6 ��ָ�������ʱ��ı仯*/
/*����ͷ�ļ�*/
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>       /*ʹ�õ�ǰʱ����Ϊ���������*/
#include <math.h>
#include <iostream>
#include <omp.h>
using namespace std;
#define MC_Time 20000
#define MAT 100           /*���������л��е���������������������߳�*/
#define Labelc6_NUM 10000  //C6���ӱ�Ƿ�
#define Labelc2_NUM  10000 //C6���ӱ�Ƿ�
#define MAT_SIcur6E 10000
//double MAT_SIcur6E = MAT*MAT;    //�����С
//double result[400][100][100];

void reaction(double C2_CON,double C6_CON,int t)
{  /*�������*/
	int x1,y1,x2,y2,x3,y3,x4,y4,y5,y6,p,coord[8],pc,labelc6[Labelc6_NUM] = {0},labelc2[Labelc2_NUM] = {0};
	unsigned seed;
	int i,j,n,nn,x,y,s,k,w,cur2,cur6,c2tag,c6tag;
	int nr,nad,nbd,nax,nbx,nat,nbt;           /*ÿ��MCS��ϵͳ���淢���ĸ���Ӧ��*/
	float aa,ba,ag,bg,ad,bd,ab,all;            /*ϵͳ�ڸ���Ӧ����*/
	float paa,pba,pag,pbg,pbd,pad,pab;      /*ϵͳ�ڸ���Ӧ����*/
	float R1;int R2,R3;                     /*���ɾ��������*/
	int surf[MAT][MAT]={0};                   /*��ʼ���������*/
	/*��ÿ MCS ��ϵͳ����Ӧ���洢�ڸ�������*/
	int tor[MC_Time];int ar[MC_Time];int br[MC_Time];int ax[MC_Time];int bx[MC_Time];int at[MC_Time];int bt[MC_Time];
	int C2[MC_Time];int C6[MC_Time];int empty[MC_Time];     /*�洢ÿһ��ѭ��ǰ�ı���״̬*/
	srand(time(NULL));                      /*��ʼ�������*/
	double C2_NUM,C6_NUM,Empty_NUM;
	Empty_NUM=MAT_SIcur6E,C2_NUM=0,C6_NUM=0;
	/*��ѭ�� MC_Time �Σ�PP ��ֵΪ���ؿ���ʱ�� MCS*/
	//int count=1;
	seed = 25234 + 17*omp_get_thread_num();
	for(n=0;n<MC_Time;n++)
	{
		//printf("runing...%d\n",count);count++;

		nr=0;nad=0;nbd=0;nax=0;nbx=0;nat=0;nbt=0;    /*ÿ��ѭ��ǰ������Ӧ������*/
		for(int i=1;i<=MAT;i++)
		{
			if(i%2!=0)
			{
				for(int j=2;j<=MAT;j=j+4)
				{
					surf[i][j]=-1;
					surf[i][j+1]=-1;
				}
			}
			else
			{
				surf[i][1]=-1;
				surf[i][100]=-1;
				for(int j=4;j<=MAT-4;j=j+4)
				{
					surf[i][j]=-1;
					surf[i][j+1]=-1;
				}
				/*printf("%.lf\n",surf[i][j]);*/
			}
		}

		empty[n]=Empty_NUM;C2[n]=C2_NUM;C6[n]=C6_NUM;       /*ȡÿ MCS ����һ��ѭ����Ϊ��� MCS �ڵı���״̬���*/
		/*һ�� MCS������ MAT_SIcur6E ��ѭ��*/
		for(nn=0;nn<MAT_SIcur6E;nn++)
		{
			pc = 0;
			s=0,k=0,c2tag=0,c6tag=0;
			/*���� label ������δʹ�õ� c2tag*/
			for(i=1;i<Labelc2_NUM;i++)
			{
				if(i==1){j=1;}
				if(i!=1){j=-c2tag;}
				if(labelc2[j]==0) {c2tag = -i-1;break;}
			}
			/*���� label ������δʹ�õ� c6tag*/
			for(i=1;i<Labelc6_NUM;i++)
			{
				if(labelc6[i]==0) {c6tag = i;break;}
			}
			/*��������ؼ���ʽ*/
			aa=10*C2_CON*0.16*0.443*Empty_NUM/MAT_SIcur6E;                                          /* C2 ��������*/
			ba=509.36*C6_CON*0.16*(0.443*Empty_NUM/MAT_SIcur6E)*(0.443*Empty_NUM/MAT_SIcur6E)*(0.443*Empty_NUM/MAT_SIcur6E);    /* C6 ��������*/
			ag=0.05*0.443*C2_NUM/MAT_SIcur6E;                                               /* C2 �Ѹ�����*/
			bg=0.2*0.443*C6_NUM/(MAT_SIcur6E*3);                                            /* C6 �Ѹ�����*/
			ad=0.0386*0.443*C2_NUM/MAT_SIcur6E;                                             /* C2 ��������*/
			bd=0.0657*0.443*C6_NUM/(MAT_SIcur6E*3);                                         /* C6 ��������*/
			ab=7*0.443*0.443*C2_NUM*C6_NUM/(MAT_SIcur6E*MAT_SIcur6E*3);                     /*˫���ӷ�Ӧ����*/
			all=aa+ba+ag+bg+ad+bd+ab;                                            /*������*/

			/*��Ӧ��Ӧ����*/
			paa=aa/all;pba=ba/all;pag=ag/all;pbg=bg/all;pad=ad/all;pbd=bd/all;pab=ab/all;

			/*����[0,1]���������*/
			R1=(float)rand()/RAND_MAX;
			R3=(float)rand()/RAND_MAX;

			/*�������Ӧ�����㷨��Ϊ�����߲���ʵ��*/
			/*1����������С������ C2,��-1 ��ʾ*/
			if(R1<=paa)
			{
				if(Empty_NUM<2) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=0);  /*�������һ�����λ��p,q��*/

				if(y%2==0)
				{
					x1=x,y1=y+1;
					x2=x-1,y2=y-1;
					x3=x+1,y3=y-1;
					if(x-1==0){x2=100;}
					if(x+1==101){x3=1;}
					if(y-1==0){y2=100;y3=100;}
					if(y+1==101){y1=1;}
					if(surf[x1][y1]==0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
					if(surf[x2][y2]==0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
					if(surf[x3][y3]==0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
					if(pc==2)
					{
						surf[coord[0]][coord[1]]=c2tag;
					}
					if(pc==4)
					{
						if(R3<=(float)1/2)
						{surf[coord[0]][coord[1]]=c2tag;}
						else
						{surf[coord[2]][coord[3]]=c2tag;}
					}
					if(pc==6)
					{
						if(R3<=(float)1/3){surf[coord[0]][coord[1]]=c2tag;}
						else if(R3<=(float)2/3){surf[coord[2]][coord[3]]=c2tag;}
						else{surf[coord[4]][coord[5]]=c2tag;}
					}
					if(pc>=2)
					{	
						surf[x][y]=c2tag;p=-c2tag-1;nax=nax+1;
						C2_NUM += 2;
						Empty_NUM -= 2;
						labelc2[p] = -1;
					}
					printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",C2_NUM,C6_NUM,Empty_NUM);
				}
				else
				{
					x1=x,y1=y-1;
					x2=x-1,y2=y+1;
					x3=x+1,y3=y+1;	
					if(x-1==0){x2=100;}
					if(x+1==101){x3=1;}
					if(y-1==0){y1=100;}
					if(y+1==101){y2=1;y3=1;}
					if(surf[x1][y1]==0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
					if(surf[x2][y2]==0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
					ifsurf[x3][y3]==0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
					if(pc==2)
					{
						surf[coord[0]][coord[1]]=c2tag;
					}
					if(pc==4)
					{
						if(R3<=(float)1/2){surf[coord[0]][coord[1]]=c2tag;}
						else{surf[coord[2]][coord[3]]=c2tag;}
					}
					if(pc==6)
					{
						if(R3<=(float)1/3){surf[coord[0]][coord[1]]=c2tag;}
						else if(R3<=(float)2/3){surf[coord[2]][coord[3]]=c2tag;}
						else{surf[coord[4]][coord[5]]=c2tag;}
					}
					if(pc>=2)
					{	
						surf[x][y]=c2tag;p=-c2tag-1;nax=nax+1;
						C2_NUM += 2;
						Empty_NUM -= 2;
						labelc2[p] = -1;
					}
				}
			}


			/*2������С���ӷ����� C6������������ʾ*/
			else if(R1<=paa+pba)
			{
				if(Empty_NUM<6) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=0);    /*�������һ�����λ*/

				x1=x-2;x2=x-1;x3=x+1;x4=x+2;
				y1=y-3;y2=y-2;y3=y-1;y4=y+1;y5=y+2;y6=y+3;
				if(x-2==0){x1=100;}
				if(x-2==-1){x1=99;x2=100;}
				if(x+2==101){x4=1;}
				if(x+2==102){x3=1;x4=2;}
				if(y-3==0){y1=100;}
				if(y-3==-1){y1=99;y2=100;}
				if(y-3==-2){y1=99;y2=99;y3=100;}
				if(y+3==101){y6=1;}
				if(y+3==102){y5=1;y6=2;}
				if(y+3==103){y4=1;y5=2;y6=3;}
				if(y%2==0)
				{

					if(surf[x][y4]==0){pc += 2;}
					if(surf[x2][y3]==0){pc += 2;}
					if(surf[x3][y3]==0){pc += 2;}
					if(pc==4)
					{
						if(surf[x2][y3]!=0)
						{
							if(surf[x4][y]==0&&surf[x4][y4]==0&&surf[x3][y5]==0)
							{
								surf[x][y]=c6tag;
								surf[x][y4]=c6tag;
								surf[x3][y5]=c6tag;
								surf[x4][y4]=c6tag;
								surf[x4][y]=c6tag;
								surf[x3][y3]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
						else if(surf[x][y4]!=0)
						{
							if(surf[x2][y2]==0&&surf[x][y1]==0&&surf[x3][y2]==0)
							{
								surf[x][y]=c6tag;
								surf[x2][y3]=c6tag;
								surf[x2][y2]=c6tag;
								surf[x][y1]=c6tag;
								surf[x3][y2]=c6tag;
								surf[x3][y3]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
						else
						{
							if(surf[x-2][y]==0&&surf[x-2][y+1]==0&&surf[x-1][y+2]==0)
							{
								surf[x-1][y-1]=c6tag;
								surf[x-2][y]=c6tag;
								surf[x-2][y+1]=c6tag;
								surf[x-1][y+2]=c6tag;
								surf[x][y+1]=c6tag;
								surf[x][y]=c6tag;nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
					}
					if(pc==6)
					{
						if(R3<=(float)1/3)
						{
							if(surf[x4][y]==0&&surf[x4][y4]==0&&surf[x3][y5]==0)
							{
								surf[x][y]=c6tag;
								surf[x][y4]=c6tag;
								surf[x3][y5]=c6tag;
								surf[x4][y4]=c6tag;
								surf[x4][y]=c6tag;
								surf[x3][y3]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
						else if(R3<=(float)2/3)
						{
							if(surf[x2][y2]==0&&surf[x][y1]==0&&surf[x3][y2]==0)
							{
								surf[x][y]=c6tag;
								surf[x2][y3]=c6tag;
								surf[x2][y2]=c6tag;
								surf[x][y1]=c6tag;
								surf[x3][y2]=c6tag;
								surf[x3][y3]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
						else
						{
							if(surf[x1][y]==0&&surf[x1][y4]==0&&surf[x2][y5]==0)
							{
								surf[x][y]=c6tag;
								surf[x2][y3]=c6tag;
								surf[x1][y]=c6tag;
								surf[x1][y4]=c6tag;
								surf[x2][y5]=c6tag;
								surf[x][y4]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}  

					}
				}
				else 
				{
					if(surf[x2][y4]==0){pc += 2;}
					if(surf[x3][y4]==0){pc += 2;}
					if(surf[x][y3]==0){pc += 2;}
					if(pc==4)
					{
						if(surf[x][y3]!=0)
						{
							if(surf[x2][y5]==0&&surf[x][y6]==0&&surf[x3][y5]==0)
							{
								surf[x][y]=c6tag;
								surf[x2][y4]=c6tag;
								surf[x2][y5]=c6tag;
								surf[x][y6]=c6tag;
								surf[x3][y5]=c6tag;
								surf[x3][y4]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
						else if(surf[x2][y4]!=0)
						{
							if(surf[x3][y2]==0&&surf[x4][y3]==0&&surf[x4][y]==0)
							{
								surf[x][y]=c6tag;
								surf[x][y3]=c6tag;
								surf[x3][y2]=c6tag;
								surf[x4][y3]=c6tag;
								surf[x4][y]=c6tag;
								surf[x3][y4]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
						else
						{
							if(surf[x1][y]==0&&surf[x1][y3]==0&&surf[x2][y2]==0)
							{
								surf[x][y]=c6tag;
								surf[x2][y4]=c6tag;
								surf[x1][y]=c6tag;
								surf[x1][y3]=c6tag;
								surf[x2][y2]=c6tag;
								surf[x][y3]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
					}
					if(pc==6)
					{
						if(R3<=(float)1/3)
						{
							if(surf[x2][y5]==0&&surf[x][y6]==0&&surf[x3][y5]==0)
							{
								surf[x][y]=c6tag;
								surf[x2][y4]=c6tag;
								surf[x2][y5]=c6tag;
								surf[x][y6]=c6tag;
								surf[x3][y5]=c6tag;
								surf[x3][y4]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
						else if(R3<=(float)2/3)
						{
							if(surf[x3][y2]==0&&surf[x4][y3]==0&&surf[x4][y]==0)
							{
								surf[x][y]=c6tag;
								surf[x][y3]=c6tag;
								surf[x3][y2]=c6tag;
								surf[x4][y3]=c6tag;
								surf[x4][y]=c6tag;
								surf[x3][y4]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}
						else
						{
							if(surf[x1][y]==0&&surf[x1][y3]==0&&surf[x2][y2]==0)
							{  
								surf[x][y]=c6tag;
								surf[x2][y4]=c6tag;
								surf[x1][y]=c6tag;
								surf[x1][y3]=c6tag;
								surf[x2][y2]=c6tag;
								surf[x][y3]=c6tag;
								nbx=nbx+1;
								C6_NUM += 6;
								Empty_NUM -= 6;
								labelc6[c6tag] = 1;
							}
						}

					}

				}
			}

			/*3������С������ C2 �Ѹ�*/
			else if(R1<=paa+pba+pag)
			{
				if(C2_NUM<2) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]>=-1);  /*���ѡ��һ C2 ����λ*/

				cur2=surf[x][y];nat=nat+1;    /* C2 �Ѹ����� 1*/
				/* C2 �Ѹ����� C2 ռ�ݵ���������λ�ϵ�ֵ����Ϊ 0*/
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==cur2){surf[i][j]=0;}
					}
				}
				C2_NUM -= 2;
				Empty_NUM += 2;
				p=-cur2-1;
				labelc2[p] = 0;
			}

			/*4������С������ C2 ����*/
			else if(R1<=paa+pba+pag+pad)
			{
				if(C2_NUM<2) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]>=-1);    /*���ѡ��һ C2 ����λ*/

				cur2=surf[x][y];nad=nad+1; /* C2 �Ѹ���ͬʱ C2 �������� 1*/
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==cur2){surf[i][j]=0;}
					}
				}
				C2_NUM -= 2;
				Empty_NUM += 2;
				p=-cur2-1;
				labelc2[p] = 0;
			}

			/*5��С���ӷ����� C6 �Ѹ�*/
			else if(R1<=paa+pba+pag+pad+pbg)
			{
				if(C6_NUM<6) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]<1);         /*���ѡ��һ C6 ����λ*/
				cur6=surf[x][y];nbt=nbt+1;    /* C6 �Ѹ����� 1*/

				/* C6 �Ѹ����� C6 ռ�ݵ���������λ�ϵ�ֵ����Ϊ 0*/
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==cur6){surf[i][j]=0;}
					}
				}

				C6_NUM -= 6;
				Empty_NUM += 6;
				labelc6[cur6] = 0;
			}

			/*6��С���ӷ����� C6 ����*/
			else if(R1<=paa+pba+pag+pad+pbg+pbd)
			{
				if(C6_NUM<6) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]<1);      /*���ѡ��һ C6 ����λ*/
				cur6=surf[x][y];nbd=nbd+1;  /* C6 �������� 1*/

				/* C6 �������� C6 ռ�ݵ���������λ�ϵ�ֵ����Ϊ 0*/
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==cur6){surf[i][j]=0;}
					}
				}

				C6_NUM -= 6;
				Empty_NUM += 6;
				labelc6[cur6] = 0;
			}

			/*7������С������ C2 ����ΧС���ӷ����� C6 ��Ӧ*/
			else
			{
				if(C6_NUM<6||C2_NUM<2) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]>=-1);    /*���ѡ��һ C2 ����λ*/

				cur2=surf[x][y];
				x1=x-2;x2=x-1;x3=x+1;x4=x+2;
				y1=y-2;y2=y-1;y3=y+1;y4=y+2;
				if(x-2==0){x1=100;}
				if(x-2==-1){x1=99;x2=100;}
				if(x+2==101){x4=1;}
				if(x+2==102){x3=1;x4=2;}
				if(y-2==0){y1=0;}
				if(y-2==-1){y1=99;y2=100;}
				if(y+2==101){y4=1;}
				if(y+2==102){y3=1;y4=2;}
				if(y%2==0)
				{
					if(surf[x][y3]==cur2)
					{
						if(surf[x2][y4]>=1){pc += 2;}
						if(surf[x3][y4]>=1){pc += 2;}
						if(surf[x3][y2]>=1){pc += 2;}
						if(surf[x2][y2]>=1){pc += 2;}
						if(pc==2)
						{
							if(surf[x2][y4]>=1)
							{
								cur6=surf[x2][y4];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x3][y4]>=1)
							{
								cur6=surf[x3][y4];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x3][y2]>=1)
							{
								cur6=surf[x3][y2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x2][y2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc==4)
						{
							if(surf[x2][y4]>=1&&surf[x3][y4]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x2][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x3][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y4]>=1&&surf[x3][y2]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x3][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x3][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y2]>=1&&surf[x2][y2]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x3][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x2][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x2][y2]>=1&&surf[x2][y4]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x2][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x2][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x2][y4]>=1&&surf[x3][y2]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x2][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x3][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else
							{
								if(R3<=0.5) 
								{
									cur6=surf[x3][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x2][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}

						}
						if(pc==6)
						{


							if(surf[x2][y4]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x3][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x3][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									printf("%d\n", x2);
									cur6=surf[x2][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
									printf("%d", y2);
								}
							}
							else if(surf[x3][y4]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x3][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x2][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x2][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y2]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x2][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x2][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x3][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else 
							{
								if(R3<=1/3) 
								{
									cur6=surf[x2][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x3][y4];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x3][y2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
						}
						if(pc==8)
						{
							if(R3<=1/4) 
							{
								cur6=surf[x2][y4];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
							else if(R3<=1/2)
							{
								cur6=surf[x3][y4];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(R3<=3/4)
							{
								cur6=surf[x3][y2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x2][y2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}    
						if(pc >= 2)
						{
							surf[x][y]=0;
							surf[x][y+1]=0;
							nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
							C6_NUM -= 6;
							C2_NUM -= 2;
							Empty_NUM += 8;
							p=-cur2-1; labelc2[p] = -1;
							printf("%d", cur6);
							labelc6[cur6] = 1;
						}

					}
					else if(surf[x3][y2]==cur2)
					{
						if(surf[x][y3]>=1){pc += 2;}
						if(surf[x4][y]>=1){pc += 2;}
						if(surf[x3][y1]>=1){pc += 2;}
						if(surf[x2][y2]>=1){pc += 2;}
						if(pc==2)
						{
							if(surf[x][y3]>=1)
							{
								cur6=surf[x][y3];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x4][y]>=1)
							{
								cur6=surf[x4][y];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x3][y1]>=1)
							{
								cur6=surf[x3][y1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x-1][y-1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc==4)
						{
							if(surf[x][y3]>=1&&surf[x4][y]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x4][y]>=1&&surf[x3][y1]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y1]>=1&&surf[x2][y2]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x2][y2]>=1&&surf[x][y3]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x-1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x][y3]>=1&&surf[x3][y1]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else
							{
								if(R3<=0.5) 
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
						}
						if(pc==6)
						{


							if(surf[x][y+1]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x+2][y]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x-1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x+1][y-2]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x-1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else 
							{
								if(R3<=1/3) 
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
						}
						if(pc==8)
						{
							if(R3<=1/4) 
							{
								cur6=surf[x][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
							else if(R3<=1/2)
							{
								cur6=surf[x+2][y];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(R3<=3/4)
							{
								cur6=surf[x+1][y-2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x-1][y-1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc >= 2)
						{
							surf[x][y]=0;
							surf[x+1][y-1]=0;
							nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
							C6_NUM -= 6;
							C2_NUM -= 2;
							Empty_NUM += 8;
							p=-cur2-1; labelc2[p] = -1;
							labelc6[cur6] = 1;
						}	
					}
					else
					{
						if(surf[x1][y]>=1){pc += 2;}
						if(surf[x][y3]>=1){pc += 2;}
						if(surf[x3][y2]>=1){pc += 2;}
						if(surf[x2][y1]>=1){pc += 2;}
						if(pc==2)
						{
							if(surf[x1][y]>=1)
							{
								cur6=surf[x-2][y];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x][y3]>=1)
							{
								cur6=surf[x][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x3][y2]>=1)
							{
								cur6=surf[x+1][y-1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x-1][y-2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc==4)
						{
							if(surf[x1][y]>=1&&surf[x][y3]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x][y3]>=1&&surf[x3][y2]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y2]>=1&&surf[x2][y1]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x+1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x2][y1]>=1&&surf[x1][y]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x1][y]>=1&&surf[x3][y2]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else
							{
								if(R3<=0.5) 
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}

						}
						if(pc==6)
						{   
							if(surf[x1][y]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x+1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x][y3]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x+1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y2]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else 
							{
								if(R3<=1/3) 
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
						}
						if(pc==8)
						{
							if(R3<=1/4) 
							{
								cur6=surf[x-2][y];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
							else if(R3<=1/2)
							{
								cur6=surf[x][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(R3<=3/4)
							{
								cur6=surf[x+1][y-1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x-1][y-2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc >= 2)
						{
							surf[x][y]=0;
							surf[x-1][y-1]=0;
							nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
							C6_NUM -= 6;
							C2_NUM -= 2;
							Empty_NUM += 8;
							p=-cur2-1; labelc2[p] = -1;
							labelc6[cur6] = 1;
						}
					}

				}
				else
				{
					if(surf[x][y-1]==cur2)
					{
						if(surf[x2][y3]>=1){pc += 2;}
						if(surf[x3][y3]>=1){pc += 2;}
						if(surf[x3][y1]>=1){pc += 2;}
						if(surf[x2][y1]>=1){pc += 2;}
						if(pc==2)
						{
							if(surf[x2][y3]>=1)
							{
								cur6=surf[x-1][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x3][y3]>=1)
							{
								cur6=surf[x+1][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x3][y1]>=1)
							{
								cur6=surf[x+1][y-2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x2][y1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc==4)
						{
							if(surf[x2][y3]>=1&&surf[x3][y3]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y3]>=1&&surf[x3][y1]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y1]>=1&&surf[x2][y1]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x2][y1]>=1&&surf[x2][y3]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x2][y3]>=1&&surf[x3][y1]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else
							{
								if(R3<=0.5) 
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}

						}
						if(pc==6)
						{


							if(surf[x-1][y+1]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x+1][y+1]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x+1][y-2]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x-1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else 
							{
								if(R3<=1/3) 
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y-2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
						}
						if(pc==8)
						{
							if(R3<=1/4) 
							{
								cur6=surf[x-1][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
							else if(R3<=1/2)
							{
								cur6=surf[x+1][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(R3<=3/4)
							{
								cur6=surf[x+1][y-2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x-1][y-2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}    
						if(pc >= 2)
						{
							surf[x][y]=0;
							surf[x][y-1]=0;
							nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
							C6_NUM -= 6;
							C2_NUM -= 2;
							Empty_NUM += 8;
							//printf("%d\n", Empty_NUM);
							p=-cur2-1; labelc2[p] = -1;
							labelc6[cur6] = 1;
						}
						//printf("%d\n", Empty_NUM);

					}
					else if(surf[x-1][y+1]==cur2)
					{
						if(surf[x2][y4]>=1){pc += 2;}
						if(surf[x3][y3]>=1){pc += 2;}
						if(surf[x][y2]>=1){pc += 2;}
						if(surf[x1][y]>=1){pc += 2;}
						if(pc==2)
						{
							if(surf[x2][y4]>=1)
							{
								cur6=surf[x-1][y+2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x3][y3]>=1)
							{
								cur6=surf[x+1][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x][y2]>=1)
							{
								cur6=surf[x][y-1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x-2][y];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc==4)
						{
							if(surf[x1][y]>=1&&surf[x2][y4]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x2][y4]>=1&&surf[x3][y3]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x-1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y3]>=1&&surf[x][y2]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x][y2]>=1&&surf[x1][y]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x1][y]>=1&&surf[x3][y3]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else
							{
								if(R3<=0.5) 
								{
									cur6=surf[x-1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
						}
						if(pc==6)
						{


							if(surf[x2][y4]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y3]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x][y2]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x-2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x-1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else 
							{
								if(R3<=1/3) 
								{
									cur6=surf[x-1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x+1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
						}
						if(pc==8)
						{
							if(R3<=1/4) 
							{
								cur6=surf[x-2][y];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
							else if(R3<=1/2)
							{
								cur6=surf[x-1][y+2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(R3<=3/4)
							{
								cur6=surf[x+1][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x][y-1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc >= 2)
						{
							surf[x][y]=0;
							surf[x-1][y+1]=0;
							nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
							C6_NUM -= 6;
							C2_NUM -= 2;
							Empty_NUM += 8;
							p=-cur2-1; labelc2[p] = -1;
							labelc6[cur6] = 1;
						}	
					}
					else
					{
						if(surf[x3][y4]>=1){pc += 2;}
						if(surf[x4][y]>=1){pc += 2;}
						if(surf[x][y2]>=1){pc += 2;}
						if(surf[x2][y3]>=1){pc += 2;}
						if(pc==2)
						{
							if(surf[x3][y4]>=1)
							{
								cur6=surf[x+1][y+2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x4][y]>=1)
							{
								cur6=surf[x+2][y];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(surf[x][y2]>=1)
							{
								cur6=surf[x][y-1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x2][y3];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc==4)
						{
							if(surf[x3][y4]>=1&&surf[x4][y]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x+1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x4][y]>=1&&surf[x][y2]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x][y2]>=1&&surf[x2][y3]>=1)
							{
								if(R3<0.5) 
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x2][y3]>=1&&surf[x3][y4]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x3][y4]>=1&&surf[x][y2]>=1)
							{
								if(R3<=0.5) 
								{
									cur6=surf[x+1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else
							{
								if(R3<=0.5) 
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}

						}
						if(pc==6)
						{   
							if(surf[x3][y4]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x4][y]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else if(surf[x][y2]<=0)
							{
								if(R3<=1/3) 
								{
									cur6=surf[x-1][y+1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x+1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
							else 
							{
								if(R3<=1/3) 
								{
									cur6=surf[x+1][y+2];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
								else if(R3<=2/3)
								{
									cur6=surf[x+2][y];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
								else
								{
									cur6=surf[x][y-1];
									for(i=0;i<MAT;i++)
									{
										for(j=0;j<MAT;j++)
										{
											if(surf[i][j]==cur6){surf[i][j]=0;}
										}
									}
								}
							}
						}
						if(pc==8)
						{
							if(R3<=1/4) 
							{
								cur6=surf[x+1][y+2];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							} /*���� C6 �漰��������λ����ת��Ϊ��λ*/
							else if(R3<=1/2)
							{
								cur6=surf[x+2][y];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else if(R3<=3/4)
							{
								cur6=surf[x][y-1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
							else
							{
								cur6=surf[x-1][y+1];
								for(i=0;i<MAT;i++)
								{
									for(j=0;j<MAT;j++)
									{
										if(surf[i][j]==cur6){surf[i][j]=0;}
									}
								}
							}
						}
						if(pc >= 2)
						{
							surf[x][y]=0;
							surf[x+1][y+1]=0;
							nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
							C6_NUM -= 6;
							C2_NUM -= 2;
							Empty_NUM += 8;
							p=-cur2-1; labelc2[p] = -1;
							labelc6[cur6] = 1;
						}
					}

				}
			}
		}
		/*printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",C2_NUM,C6_NUM,Empty_NUM);*/
		/*printf("nr:%.d nad:%.d nbd:%.d nax:%.d nbx:%.d nat:%.d nbt:%.d \n",nr,nad,nbd,nax,nbx,nat,nbt);*/

		/*���� MCS ��ϵͳ����Ӧ���洢�ڸ�������*/
		tor[n]=nr;  //tor �洢˫���ӷ�Ӧ��
		ar[n]=nad;  //ar �洢С���ӳ�����
		br[n]=nbd;  //br �洢�����������
		ax[n]=nax;  //ax �洢С����������
		bx[n]=nbx;  //bx �洢�����������
		at[n]=nat;  //at �洢С�����Ѹ���
		bt[n]=nbt;  //bt �洢������Ѹ���
	}
	/*���淴Ӧ���*/
	/*for(i=0; i<MAT; i++)
	for(j=0; j<MAT; j++)
	result[t][i][j] = surf[i][j];*/
}

int main()
{
	//	memset(result, 0, 400*100*100); //��ʼ������������ά����
	/*����400��a,b���*/
	double C2_CON[5] = {0.1,0.2,0.3,0.4,0.5};
	double Rate[80];
	time_t start,finish;
	for(int i=0; i<80; i++)
	{
		Rate[i] = 0.2*i;
	}
	/*��ʼ���㷴Ӧ*/
	omp_set_num_threads(32); 
	start = time(NULL);
	for(int i=0; i<5; i++)
	{
#pragma omp parallel for
		for(int j=0; j<80; j++)
		{
			reaction(C2_CON[i],C2_CON[i]*Rate[j],i*j);
		}
	}
	finish = time(NULL);
	cout << "time:" << finish-start << endl;
	return 0;
}
