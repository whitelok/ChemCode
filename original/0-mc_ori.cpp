/*定义头文件*/
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>       /*使用当前时间作为随机数种子*/
#include <math.h>
#include <iostream>
using namespace std;
#define MAT 100           /*方形网格行或列的网格数，代表反应表面边长*/
#define MC_Time 20000      /*蒙特卡罗总运行时间*/

int result[400][100][100];

void reaction(double C2_CON,double C6_CON,int t)
{  /*定义变量*/
	int x1,y1,x2,y2,x3,y3,x4,y4,coord[8],pc;
	
	int i,j,n,nn,p,q,n1,n2,n3,s,k,w,o,l;
	int nr,nad,nbd,nax,nbx,nat,nbt;           /*每个MCS内系统表面发生的各反应数*/
	float aa,ba,ag,bg,ad,bd,ab,all;            /*系统内各反应速率*/
	float paa,pba,pag,pbg,pbd,pad,pab;      /*系统内各反应概率*/
	float r1,r3;int r2;                     /*生成均匀随机数*/
	int surf[MAT][MAT]={0};                   /*初始化基体表面*/
	/*将每 MCS 内系统各反应数存储在各数组内*/
	int tor[MC_Time];int ar[MC_Time];int br[MC_Time];int ax[MC_Time];int bx[MC_Time];int at[MC_Time];int bt[MC_Time];
	int small[MC_Time];int large[MC_Time];int empty[MC_Time];     /*存储每一次循环前的表面状态*/
	srand(time(NULL));                      /*初始化随机数*/
	/*总循环 MC_Time 次，MC_Time 的值为蒙特卡罗时间 MCS*/
	//int count=1;
	for(n=0;n<MC_Time;n++)
	{
	    //printf("runing...%d\n",count);count++;

		nr=0;nad=0;nbd=0;nax=0;nbx=0;nat=0;nbt=0;    /*每次循环前将各反应数归零*/
		/*一个 MCS，包含 MAT*MAT 次循环*/
		for(nn=0;nn<MAT*MAT;nn++)
		{
			pc = 0;
			n1=0,n2=0;n3=0,s=0,k=0,l=1;
			/*扫描整个吸附表面，将表面的空位数记为 n1，线性小分子烃 C2 数记为 n2，小分子芳香烃 C6 数记为 n3*/
			for(i=0;i<MAT;i++)
			{
				for(j=0;j<MAT;j++)
				{
					if(surf[i][j]==0){n1=n1+1;}  /*统计基体表面空位数*/
					else if(surf[i][j]==-1){n2=n2+1;}/*统计基体表面线性小分子烃 C2 数*/
					else{n3=n3+1;}              /*统计基体表面小分子芳香烃 C6 数*/
				}
			}
			/*取每 MCS 的第一次循环作为这次 MCS 内的表面状态输出*/
			if(nn==0)
			{
				empty[n]=n1;small[n]=n2;large[n]=n3;
			}
			/*将表面吸附的 C6 以不同的正整数表示并加以区别*/
			/*查找 surf 数组中不存在的最小正整数 l*/
			for(i=0;i<MAT;i++)
			{
				for(j=0;j<MAT;j++)
				{
					if(surf[i][j]==l) {l=l+1;i=0;}
				}
			}

			/*以下是相关计算式*/
			aa=10*C2_CON*0.16*0.443*n1/2500;                                          /* C2 吸附速率*/
			ba=509.36*C6_CON*0.16*(0.443*n1/2500)*(0.443*n1/2500)*(0.443*n1/2500);    /* C6 吸附速率*/
			ag=0.05*0.443*n2/2500;                                               /* C2 脱附速率*/
			bg=0.2*0.443*n3/(2500*3);                                            /* C6 脱附速率*/
			ad=0.0386*0.443*n2/2500;                                             /* C2 沉积速率*/
			bd=0.0657*0.443*n3/(2500*3);                                         /* C6 沉积速率*/
			ab=7*0.443*0.443*n2*n3/(2500*2500*3);                                /*双分子反应速率*/
			all=aa+ba+ag+bg+ad+bd+ab;                                            /*总速率*/

			/*对应反应概率*/
			paa=aa/all;pba=ba/all;pag=ag/all;pbg=bg/all;pad=ad/all;pbd=bd/all;pab=ab/all;

			/*产生[0,1]均匀随机数*/
			r1=(float)rand()/RAND_MAX;
			r3=(float)rand()/RAND_MAX;

			/*具体各方应过程算法分为以下七步来实现*/
			/*1、吸附线性小分子烃 C2,以-1 表示*/
			if(r1<=paa)
			{
				if( n1<1 ) continue;
				do{
					r2=rand()%(MAT*MAT)+1;
					p=(r2-1)/MAT;
					q=(r2-1)%MAT;
				}while(surf[p][q]!=0);  /*随机生成一表面空位（p,q）*/
				surf[p][q]=-1;nax=nax+1;    /*在此空位上吸附 C2，C2 吸附数加 1 */
			}

			/*2、吸附小分子芳香烃 C6，以正整数表示*/
			else if(r1<=paa+pba)
			{
				if( n1<3 ) continue;
				do{
					r2=rand()%(MAT*MAT)+1;
					p=(r2-1)/MAT;q=(r2-1)%MAT;
				}while(surf[p][q]!=0);    /*随机生成一表面空位*/

				x1=p-1,y1=q;
				x2=p,y2=q-1;
				x3=p+1,y3=q;
				x4=p,y4=q+1;
				if(x1>=0&&x1<MAT&&y1>=0&&y1<MAT&&surf[x1][y1]==0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<MAT&&y2>=0&&y2<MAT&&surf[x2][y2]==0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<MAT&&y3>=0&&y3<MAT&&surf[x3][y3]==0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<MAT&&y4>=0&&y4<MAT&&surf[x4][y4]==0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==4)
				{
					surf[p][q]=l;
					surf[coord[0]][coord[1]]=l;
					surf[coord[2]][coord[3]]=l;
					nbx=nbx+1;
				}
				if(pc==6)
				{
					surf[p][q]=l;nbx=nbx+1;
					/*在此三邻位位上随机选择两邻位与此点形成三空位吸附 C6，同样吸附的C6 以 l 表示，同时 C6 吸附数加 1*/
					if(r3<=(float)1/3){surf[coord[0]][coord[1]]=l;surf[coord[2]][coord[3]]=l;}
					else if(r3<=(float)2/3){surf[coord[2]][coord[3]]=l;surf[coord[4]][coord[5]]=l;}
					else{surf[coord[0]][coord[1]]=l;surf[coord[4]][coord[5]]=l;}
				}
				if(pc==8)
				{
					surf[p][q]=l;nbx=nbx+1;
					/*在此三邻位位上随机选择两邻位与此点形成三空位吸附 C6，同样吸附的C6 以 l 表示，同时 C6 吸附数加 1*/
					if(r3<=(float)1/4){surf[coord[0]][coord[1]]=l;surf[coord[2]][coord[3]]=l;}
					else if(r3<=(float)2/4){surf[coord[2]][coord[3]]=l;surf[coord[4]][coord[5]]=l;}
					else if(r3<=(float)3/4){surf[coord[4]][coord[5]]=l;surf[coord[6]][coord[7]]=l;}
					else{surf[coord[0]][coord[1]]=l;surf[coord[6]][coord[7]]=l;}
				}
			}

			/*3、线性小分子烃 C2 脱附*/
			else if(r1<=paa+pba+pag)
			{
				if( n2<1 ) continue;
				do{
					r2=rand()%(MAT*MAT)+1;
					p=(r2-1)/MAT;q=(r2-1)%MAT;
				}while(surf[p][q]!=-1);  /*随机选择一 C2 吸附位*/
				surf[p][q]=0;nat=nat+1;  /* C2 脱附，同时 C2 脱附数加 1*/
			}

			/*4、线性小分子烃 C2 沉积*/
			else if(r1<=paa+pba+pag+pad)
			{
				if( n2<1 ) continue;
				do{
					r2=rand()%(MAT*MAT)+1;
					p=(r2-1)/MAT;q=(r2-1)%MAT;
				}while(surf[p][q]!=-1);    /*随机选择一 C2 吸附位*/
				surf[p][q]=0;nad=nad+1;  /* C2 脱附，同时 C2 沉积数加 1*/
			}

			/*5、小分子芳香烃 C6 脱附*/
			else if(r1<=paa+pba+pag+pad+pbg)
			{
				if( n3<3 ) continue;
				do{
					r2=rand()%(MAT*MAT)+1;
					p=(r2-1)/MAT;q=(r2-1)%MAT;
				}while(surf[p][q]<1);         /*随机选择一 C6 吸附位*/
				o=surf[p][q];surf[p][q]=0;nbt=nbt+1;    /* C6 脱附数加 1*/
				/* C6 脱附：将 C6 占据的三个吸附位上的值更新为 0*/
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==o){surf[i][j]=0;}
					}
				}
			}

			/*6、小分子芳香烃 C6 沉积*/
			else if(r1<=paa+pba+pag+pad+pbg+pbd)
			{
				if( n3<3 ) continue;
				do{
					r2=rand()%(MAT*MAT)+1;
					p=(r2-1)/MAT;q=(r2-1)%MAT;
				}while(surf[p][q]<1);      /*随机选择一 C6 吸附位*/
				o=surf[p][q];nbd=nbd+1;  /* C6 沉积数加 1*/
				/* C6 沉积：将 C6 占据的三个吸附位上的值更新为 0*/
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
					if(surf[i][j]==o){surf[i][j]=0;}
					}
				}
			}

			/*7、线性小分子烃 C2 与周围小分子芳香烃 C6 反应*/
			else
			{
				if( n2<1 || n3<3 ) continue;
				do{
					r2=rand()%(MAT*MAT)+1;
					p=(r2-1)/MAT;q=(r2-1)%MAT;
				}while(surf[p][q]!=-1);    /*随机选择一 C2 吸附位*/

				x1=p-1,y1=q;
				x2=p,y2=q-1;
				x3=p+1,y3=q;
				x4=p,y4=q+1;
				if(x1>=0&&x1<MAT&&y1>=0&&y1<MAT&&surf[x1][y1]>0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<MAT&&y2>=0&&y2<MAT&&surf[x2][y2]>0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<MAT&&y3>=0&&y3<MAT&&surf[x3][y3]>0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<MAT&&y4>=0&&y4<MAT&&surf[x4][y4]>0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==2)
				{/*如果存在一个 C6 吸附位，则与其发生双分子反应*/
					surf[p][q]=0;    /*将此 C2 吸附位转变为空位*/
					o=surf[coord[0]][coord[1]];
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==o){surf[i][j]=0;}
						}
					}    /*将此 C6 涉及的三吸附位重新转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
				}
				if(pc==4)
				{/*如果两邻位都被 C6 吸附*/
					surf[p][q]=0;    /*将此 C2 吸附位转变为空位*/
					/*随机选择其中一邻位*/
					if(r3<0.5) o=surf[coord[0]][coord[1]];
					else o=surf[coord[2]][coord[3]];
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==o){surf[i][j]=0;}
						}
					}    /*将此 C6 涉及的三吸附位重新转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
				}
				if(pc==6)
				{/*如果三邻位都被 C6 吸附*/
					surf[p][q]=0;    /*将此 C2 吸附位转变为空位*/
					/*在此三邻位位上随机选择一邻位*/
					if(r3<=(float)1/3) o=surf[coord[0]][coord[1]];
					else if(r3<=(float)2/3) o=surf[coord[2]][coord[3]];
					else o=surf[coord[4]][coord[5]];
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==o){surf[i][j]=0;}
						}
					}    /*将此 C6 涉及的三吸附位重新转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
				}
				if(pc==8)
				{/*如果四邻位都被 C6 吸附*/
					surf[p][q]=0;    /*将此 C2 吸附位转变为空位*/
					/*在此四邻位位上随机选择一邻位*/
					if(r3<=(float)1/4) o=surf[coord[0]][coord[1]];
					else if(r3<=(float)2/4) o=surf[coord[2]][coord[3]];
					else if(r3<=(float)3/4) o=surf[coord[4]][coord[5]];
					else o=surf[coord[6]][coord[7]];
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==o){surf[i][j]=0;}
						}
					}    /*将此 C6 涉及的三吸附位重新转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
				}
			}
		}
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
	memset(result, 0, sizeof(result)); //初始化保存结果的三维数组
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