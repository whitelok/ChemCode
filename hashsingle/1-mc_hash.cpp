/*以下为求解 R=C6_CON/C2_CON（C2_CON,C6_CON）时，表面吸附 C2 和 C6 组分覆盖率随时间的变化*/
/*定义头文件*/
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>       /*使用当前时间作为随机数种子*/
#include <math.h>
using namespace std;
#define MC_Time 20000
#define MAT 100           /*方形网格行或列的网格数，代表吸附表面边长*/
#define Label_NUM 10000  //C6分子标记符

double MAT_SIZE = MAT*MAT;    //矩阵大小
double result[400][100][100];

void reaction(double C2_CON,double C6_CON,int t)
{  /*定义变量*/
	int x1,y1,x2,y2,x3,y3,x4,y4,coord[8],pc,label[Label_NUM] = {0};

	int i,j,n,nn,x,y,s,k,w,cur,tag;
	int nr,nad,nbd,nax,nbx,nat,nbt;           /*每个MCS内系统表面发生的各反应数*/
	float aa,ba,ag,bg,ad,bd,ab,all;            /*系统内各反应速率*/
	float paa,pba,pag,pbg,pbd,pad,pab;      /*系统内各反应概率*/
	float R1;int R2,R3;                     /*生成均匀随机数*/
	int surf[MAT][MAT]={0};                   /*初始化基体表面*/
	/*将每 MCS 内系统各反应数存储在各数组内*/
	int tor[MC_Time];int ar[MC_Time];int br[MC_Time];int ax[MC_Time];int bx[MC_Time];int at[MC_Time];int bt[MC_Time];
	int C2[MC_Time];int C6[MC_Time];int empty[MC_Time];     /*存储每一次循环前的表面状态*/
	srand(time(NULL));                      /*初始化随机数*/
	double C2_NUM,C6_NUM,Empty_NUM;
	Empty_NUM=MAT_SIZE,C2_NUM=0,C6_NUM=0;
	/*总循环 MC_Time 次，PP 的值为蒙特卡罗时间 MCS*/
	//int count=1;
	for(n=0;n<MC_Time;n++)
	{
	    //printf("runing...%d\n",count);count++;

		nr=0;nad=0;nbd=0;nax=0;nbx=0;nat=0;nbt=0;    /*每次循环前将各反应数归零*/
		empty[n]=Empty_NUM;C2[n]=C2_NUM;C6[n]=C6_NUM;       /*取每 MCS 的上一次循环作为这次 MCS 内的表面状态输出*/
		/*一个 MCS，包含 MAT_SIZE 次循环*/
		for(nn=0;nn<MAT_SIZE;nn++)
		{
			pc = 0;
			s=0,k=0,tag=0;
			/*查找 label 数组中未使用的 tag*/
			for(i=1;i<Label_NUM;i++)
			{
				if(label[i]==0) {tag = i;break;}
			}

			/*以下是相关计算式*/
			aa=10*C2_CON*0.16*0.443*Empty_NUM/MAT_SIZE;                                          /* C2 吸附速率*/
			ba=509.36*C6_CON*0.16*(0.443*Empty_NUM/MAT_SIZE)*(0.443*Empty_NUM/MAT_SIZE)*(0.443*Empty_NUM/MAT_SIZE);    /* C6 吸附速率*/
			ag=0.05*0.443*C2_NUM/MAT_SIZE;                                               /* C2 脱附速率*/
			bg=0.2*0.443*C6_NUM/(MAT_SIZE*3);                                            /* C6 脱附速率*/
			ad=0.0386*0.443*C2_NUM/MAT_SIZE;                                             /* C2 沉积速率*/
			bd=0.0657*0.443*C6_NUM/(MAT_SIZE*3);                                         /* C6 沉积速率*/
			ab=7*0.443*0.443*C2_NUM*C6_NUM/(MAT_SIZE*MAT_SIZE*3);                     /*双分子反应速率*/
			all=aa+ba+ag+bg+ad+bd+ab;                                            /*总速率*/

			/*对应反应概率*/
			paa=aa/all;pba=ba/all;pag=ag/all;pbg=bg/all;pad=ad/all;pbd=bd/all;pab=ab/all;

			/*产生[0,1]均匀随机数*/
			R1=(float)rand()/RAND_MAX;
			R3=(float)rand()/RAND_MAX;

			/*具体各方应过程算法分为以下七步来实现*/
			/*1、吸附线性小分子烃 C2,以-1 表示*/
			if(R1<=paa)
			{
				if(Empty_NUM<1) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=0);  /*随机生成一表面空位（p,q）*/
				surf[x][y]=-1;nax=nax+1;    /*在此空位上吸附 C2，C2 吸附数加 1 */
				C2_NUM += 1;
				Empty_NUM -= 1;
			}

			/*2、吸附小分子芳香烃 C6，以正整数表示*/
			else if(R1<=paa+pba)
			{
				if(Empty_NUM<3) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=0);    /*随机生成一表面空位*/

				x1=x-1,y1=y; //上
				x2=x,y2=y-1; //左
				x3=x+1,y3=y; //下
				x4=x,y4=y+1; //右
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
					/*在此三邻位位上随机选择两邻位与此点形成三空位吸附 C6，同样吸附的C6 以 tag 表示，同时 C6 吸附数加 1*/
					if(R3<=(float)1/3) {surf[coord[0]][coord[1]]=tag;surf[coord[2]][coord[3]]=tag;}
					else if(R3<=(float)2/3) {surf[coord[2]][coord[3]]=tag;surf[coord[4]][coord[5]]=tag;}
					else{surf[coord[0]][coord[1]]=tag;surf[coord[4]][coord[5]]=tag;}
				}
				if(pc==8)
				{
					/*在此四邻位位上随机选择两邻位与此点形成三空位吸附 C6，同样吸附的C6 以 tag 表示，同时 C6 吸附数加 1*/
					if(R3<=(float)1/4) {surf[coord[0]][coord[1]]=tag;surf[coord[2]][coord[3]]=tag;}
					else if(R3<=(float)2/4) {surf[coord[2]][coord[3]]=tag;surf[coord[4]][coord[5]]=tag;}
					else if(R3<=(float)3/4) {surf[coord[4]][coord[5]]=tag;surf[coord[6]][coord[7]]=tag;}
					else{surf[coord[0]][coord[1]]=tag;surf[coord[6]][coord[7]]=tag;}
				}
				if(pc>=2)
				{	
					surf[x][y]=tag;nbx=nbx+1;
					C6_NUM += 3;
					Empty_NUM -= 3;
					label[tag] = 1;
				}
			}

			/*3、线性小分子烃 C2 脱附*/
			else if(R1<=paa+pba+pag)
			{
				if(C2_NUM<1) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);  /*随机选择一 C2 吸附位*/
				surf[x][y]=0;nat=nat+1;  /* C2 脱附，同时 C2 脱附数加 1*/
				C2_NUM -= 1;
				Empty_NUM += 1;
			}

			/*4、线性小分子烃 C2 沉积*/
			else if(R1<=paa+pba+pag+pad)
			{
				if(C2_NUM<1) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);    /*随机选择一 C2 吸附位*/
				surf[x][y]=0;nad=nad+1;  /* C2 脱附，同时 C2 沉积数加 1*/
				C2_NUM -= 1;
				Empty_NUM += 1;
			}

			/*5、小分子芳香烃 C6 脱附*/
			else if(R1<=paa+pba+pag+pad+pbg)
			{
				if(C6_NUM<3) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]<1);         /*随机选择一 C6 吸附位*/
				cur=surf[x][y];nbt=nbt+1;    /* C6 脱附数加 1*/
				
				/* C6 脱附：将 C6 占据的三个吸附位上的值更新为 0*/
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==cur){surf[i][j]=0;}
					}
				}
				
				C6_NUM -= 3;
				Empty_NUM += 3;
				label[cur] = 0;
			}

			/*6、小分子芳香烃 C6 沉积*/
			else if(R1<=paa+pba+pag+pad+pbg+pbd)
			{
				if(C6_NUM<3) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]<1);      /*随机选择一 C6 吸附位*/
				cur=surf[x][y];nbd=nbd+1;  /* C6 沉积数加 1*/
				
				/* C6 沉积：将 C6 占据的三个吸附位上的值更新为 0*/
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==cur){surf[i][j]=0;}
					}
				}
				
				C6_NUM -= 3;
				Empty_NUM += 3;
				label[cur] = 0;
			}

			/*7、线性小分子烃 C2 与周围小分子芳香烃 C6 反应*/
			else
			{
				if(C6_NUM<3||C2_NUM<1) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);    /*随机选择一 C2 吸附位*/

				x1=x-1,y1=y;
				x2=x,y2=y-1;
				x3=x+1,y3=y;
				x4=x,y4=y+1;
				if(x1>=0&&x1<MAT&&y1>=0&&y1<MAT&&surf[x1][y1]>0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<MAT&&y2>=0&&y2<MAT&&surf[x2][y2]>0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<MAT&&y3>=0&&y3<MAT&&surf[x3][y3]>0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<MAT&&y4>=0&&y4<MAT&&surf[x4][y4]>0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==2)
				{/*如果存在一个 C6 吸附位，则与其发生双分子反应*/
					cur=surf[coord[0]][coord[1]];
					/*将此 C6 涉及的三吸附位重新转变为空位*/
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==cur){surf[i][j]=0;}
						}
					}    
				}
				if(pc==4)
				{/*如果两邻位都被 C6 吸附*/
					/*随机选择其中一邻位*/
					if(R3<=(float)1/2) cur=surf[coord[0]][coord[1]];
					else cur=surf[coord[2]][coord[3]];
					/*将此 C6 涉及的三吸附位重新转变为空位*/
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==cur){surf[i][j]=0;}
						}
					}    
				}
				if(pc==6)
				{/*如果三邻位都被 C6 吸附*/
					/*在此三邻位位上随机选择一邻位*/
					if(R3<=(float)1/3) cur=surf[coord[0]][coord[1]];
					else if(R3<=(float)2/3) cur=surf[coord[2]][coord[3]];
					else cur=surf[coord[4]][coord[5]];
					/*将此 C6 涉及的三吸附位重新转变为空位*/
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==cur){surf[i][j]=0;}
						}
					}    
				}
				if(pc==8)
				{/*如果四邻位都被 C6 吸附*/
					/*在此四邻位位上随机选择一邻位*/
					if(R3<=(float)1/4) cur=surf[coord[0]][coord[1]];
					else if(R3<=(float)2/4) cur=surf[coord[2]][coord[3]];
					else if(R3<=(float)3/4) cur=surf[coord[4]][coord[5]];
					else cur=surf[coord[6]][coord[7]];
					/*将此 C6 涉及的三吸附位重新转变为空位*/
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==cur){surf[i][j]=0;}
						}
					}    
				}
				if(pc >= 2)
				{
					surf[x][y]=0;    /*将此 C2 吸附位转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
					C6_NUM -= 3;
					C2_NUM -= 1;
					Empty_NUM += 4;
					label[cur] = 0;
				}
			}
		}
		//printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",C2_NUM,C6_NUM,Empty_NUM);

/*将此 MCS 内系统各反应数存储在各数组内*/
		tor[n]=nr;  //tor 存储双分子反应数
		ar[n]=nad;  //ar 存储小分子沉积数
		br[n]=nbd;  //br 存储大分子吸附数
		ax[n]=nax;  //ax 存储小分子吸附数
		bx[n]=nbx;  //bx 存储大分子吸附数
		at[n]=nat;  //at 存储小分子脱附数
		bt[n]=nbt;  //bt 存储大分子脱附数
	}
	/*保存反应结果*/
	for(i=0; i<MAT; i++)
		for(j=0; j<MAT; j++)
			result[t][i][j] = surf[i][j];
}

int main()
{
	memset(result, 0, 400*100*100); //初始化保存结果的三维数组
    /*生成400组a,b组合*/
    double C2_CON[5] = {0.1,0.2,0.3,0.4,0.5};
    double Rate[80];
    for(int i=0; i<80; i++)
    {
    	Rate[i] = 0.2*i;
    }
    /*开始计算反应*/
    for(int i=0; i<5; i++)
    {
    	for(int j=0; j<80; j++)
    	{
        	reaction(C2_CON[i],C2_CON[i]*Rate[j],i*j);
    	}
	}
    return 0;
}