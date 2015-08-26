#include <stdlib.h>
#include <stdio.h>
#include <time.h>       /*ʹ�õ�ǰʱ����Ϊ���������*/
#include <math.h>
#include <iostream>
#include <vector>
#include <memory.h>
using namespace std;

#define MAT 100
#define Label_NUM 10000  //C6���ӱ�Ƿ�
#define MC_Time 20000
#define MAT_SIZE 10000

#define CALC_AA(x,y) y*10*0.16*0.443*x/MAT_SIZE
#define CALC_BA(x,y) y*509.36*0.16*pow((0.443*x/MAT_SIZE), 3)
#define CALC_AG(y) 0.05*0.443*y/MAT_SIZE
#define CALC_BG(z) 0.2*0.443*z/(MAT_SIZE*3)
#define CALC_AD(y) 0.0386*0.443*y/MAT_SIZE;
#define CALC_BD(z) 0.0657*0.443*z/(MAT_SIZE*3)
#define CALC_AB(y,z) 7*0.443*0.443*y*z/(MAT_SIZE*MAT_SIZE*3)

int sum1=0;
int currentC2=0,currentC6=0;
int cC6=0,cC2=0;
vector<int> locateC2;
vector<int> locateC6;
int surface[MAT][MAT]={0};
int double_react_x, double_react_y;
int C2_C6_REACT_PATTERN[][2] = {
	{0,1},{1,0},{-1,0},{0,-1}};

inline int canC2_C6_React(int C2_NUM, int C6_NUM){
	int i,t,x,y;
	clock_t start, finish;
	int locateC2C6[4][2]={0};
	int p=0;
	//start = clock();
	//genSurface(C2_NUM-currentC2, C6_NUM-currentC6);
	t = rand() % (C2_NUM);
	t*=2;
	double_react_x = locateC2[t];
	//cout << " double_react_x: " << double_react_x + C2_C6_REACT_PATTERN[i][0] <<endl;
	double_react_y = locateC2[t+1];
	//finish = clock();
	//sum1+=(double)(finish - start);
	//cout<<"gensurface time:"<<(double)(finish - start) / CLOCKS_PER_SEC<<endl;
	for(i = 0; i < 4; i++){
		if(double_react_x + C2_C6_REACT_PATTERN[i][0] >= 0 &&
				double_react_x + C2_C6_REACT_PATTERN[i][0] < MAT &&
				double_react_y + C2_C6_REACT_PATTERN[i][1] >=0 &&
				double_react_y + C2_C6_REACT_PATTERN[i][1] < MAT &&
				surface[double_react_x + C2_C6_REACT_PATTERN[i][0]][double_react_y + C2_C6_REACT_PATTERN[i][1]]>0){
			locateC2C6[p][0] = double_react_x + C2_C6_REACT_PATTERN[i][0];
			locateC2C6[p][1] = double_react_y + C2_C6_REACT_PATTERN[i][1];
			p++;
		}
	}
	if(p==0)
		return 0;
	//cout << "OMG" << surface[double_react_x + C2_C6_REACT_PATTERN[0][0]][double_react_y + C2_C6_REACT_PATTERN[0][1]] << endl;
	surface[double_react_x][double_react_y] = 0;
	locateC2.erase(locateC2.begin()+t);
	locateC2.erase(locateC2.begin()+t);

	t = rand() % p;
	x = locateC2C6[t][0];
	//cout << "x: " << x << endl;
	y = locateC2C6[t][1];
	//cout << "y: " << y << endl;
	i = surface[x][y];
	if(i==1)
		i=0;
	for(;i<locateC6.size();i+=6){
		if(x==locateC6[i]&&y==locateC6[i+1])
			break;
	}
	/*vector<int>::iterator it,it2;
	it2=locateC6.end();
	for(it = locateC6.begin();it!=it2;it+=2){
		if(x == *it && y == *(it+1))
			break;
	}
	i = it-locateC6.begin();*/
	if(i % 6 == 0)
		t=i;
	else if(i % 6 == 2)
		t = i - 2;
	else if(i % 6 == 4)
		t = i - 4;
	else
		cout << "impossible" << endl;
	//cout << locateC6.size() << "ccc" << i << endl;
	x = locateC6[t];
	y = locateC6[t+1];
	locateC6.erase(locateC6.begin()+t);
	locateC6.erase(locateC6.begin()+t);
	surface[x][y] = 0;
	x = locateC6[t];
	y = locateC6[t+1];
	surface[x][y] = 0;
	locateC6.erase(locateC6.begin()+t);
	locateC6.erase(locateC6.begin()+t);
	x = locateC6[t];
	y = locateC6[t+1];
	surface[x][y] = 0;
	locateC6.erase(locateC6.begin()+t);
	locateC6.erase(locateC6.begin()+t);
	return 1;



}

void test1(double C2_CON,double C6_CON)
{  /*�������*/
	int x1,y1,x2,y2,x3,y3,x4,y4,coord[8],pc;
	int t;
	
	clock_t start, finish;	
	int i,j,n,nn,x,y,s,k,w,cur,tag;
	float aa,ba,ag,bg,ad,bd,ab,all;            /*ϵͳ�ڸ���Ӧ����*/
	float paa,pba,pag,pbg,pbd,pad,pab;      /*ϵͳ�ڸ���Ӧ����*/
	float R1,R3;int R2;                     /*���ɾ��������*/
	//int surf[MAT][MAT]={0};                   /*��ʼ���������*/
	srand(time(NULL));                      /*��ʼ�������*/
	/*��ѭ�� MC_Time �Σ�PP ��ֵΪ���ؿ���ʱ�� MCS*/
	int Empty_NUM,C2_NUM,C6_NUM;
	//double l1,l2,l3;
	//printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",C2_NUM,C6_NUM,Empty_NUM);
	C2_NUM=0;C6_NUM=0;Empty_NUM=MAT_SIZE;
	for(n=0;n<20000;n++)
	{
		//cout << "sum1   " << sum1;
		//cout << "C2: " << C2_NUM << "C6: " << C6_NUM << "Empty: " << Empty_NUM << endl;
		//sum1=0;
		/*һ�� MCS������ MAT*MAT ��ѭ��*/
		//start=clock();
		for(nn=0;nn<10000;nn++)
		{

						
			
			pc = 0;
			s=0,k=0;
			/*ɨ�������������棬������Ŀ�λ����Ϊ n1������С������ C2 ����Ϊ n2��С���ӷ����� C6 ����Ϊ C6_NUM*/
			//Empty_NUM=0,C2_NUM=0,C6_NUM=0;
			//for(i=0;i<MAT;i++)
			//{
				//for(j=0;j<MAT;j++)
				//{
					//if(surf[i][j]==0) {Empty_NUM=Empty_NUM+1;}  /*ͳ�ƻ�������λ��*/
					//else if(surf[i][j]==-1) {C2_NUM=C2_NUM+1;}/*ͳ�ƻ����������С������ C2 ��*/
				//	else {C6_NUM=C6_NUM+1;}              /*ͳ�ƻ������С���ӷ����� C6 ��*/
				//}
			//}
			
			//printf("%d >>>>>>C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",nn,l1,l2,l3);
			//cout << "C2: " << C2_NUM << "C6: " << C6_NUM << "Empty: " << Empty_NUM << endl;
			//if(l2!=C6_NUM)break;

			//printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",C2_NUM,C6_NUM,Empty_NUM);

			/*������������ C6 �Բ�ͬ����������ʾ����������*/
			/*���� surf �����в����ڵ���С������ tag*/
/*loop:			
			for(i=0;i<MAT;i++)
			{
				for(j=0;j<MAT;j++)
				{
					if(surf[i][j]==tag) {tag=tag+1;i=0;goto loop;}
				}
			}*/

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
			R1=(float)rand()/RAND_MAX;
			R3=(float)rand()/RAND_MAX;

			/*�������Ӧ�����㷨��Ϊ�����߲���ʵ��*/
			/*1����������С������ C2,��-1 ��ʾ*/
			if(R1<=paa)
			{
				if( Empty_NUM<1 ) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surface[x][y]!=0);  /*�������һ�����λ��x,y��*/
				surface[x][y]=-1;   /*�ڴ˿�λ������ C2��C2 �������� 1 */
				locateC2.push_back(x);
				locateC2.push_back(y);
				C2_NUM++;
				Empty_NUM--;
			}
			/*2������С���ӷ����� C6������������ʾ*/
			else if(R1<=paa+pba)
			{
				if( Empty_NUM<3 ) continue;
				do{
					R2=rand()%(MAT*MAT)+1;
					x=(R2-1)/MAT;y=(R2-1)%MAT;
				}while(surface[x][y]!=0);    /*�������һ�����λ*/

				x1=x-1,y1=y;
				x2=x,y2=y-1;
				x3=x+1,y3=y;
				x4=x,y4=y+1;
				if(x1>=0&&x1<MAT&&y1>=0&&y1<MAT&&surface[x1][y1]==0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<MAT&&y2>=0&&y2<MAT&&surface[x2][y2]==0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<MAT&&y3>=0&&y3<MAT&&surface[x3][y3]==0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<MAT&&y4>=0&&y4<MAT&&surface[x4][y4]==0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==4)
				{
					surface[x][y]=1;
					surface[coord[0]][coord[1]]=2;
					surface[coord[2]][coord[3]]=4;
					C6_NUM+=3;
					Empty_NUM-=3;
					locateC6.push_back(x);
					locateC6.push_back(y);
					locateC6.push_back(coord[0]);
					locateC6.push_back(coord[1]);
					locateC6.push_back(coord[2]);
					locateC6.push_back(coord[3]);
					
				}
				if(pc==6)
				{
					surface[x][y]=1;
					/*�ڴ�����λλ�����ѡ������λ��˵��γ�����λ���� C6��ͬ��������C6 �� tag ��ʾ��ͬʱ C6 �������� 1*/
					if(R3<=(float)1/3){surface[coord[0]][coord[1]]=2;surface[coord[2]][coord[3]]=4;locateC6.push_back(x);
					locateC6.push_back(y);
					locateC6.push_back(coord[0]);
					locateC6.push_back(coord[1]);
					locateC6.push_back(coord[2]);
					locateC6.push_back(coord[3]);}
					else if(R3<=(float)2/3){surface[coord[2]][coord[3]]=2;surface[coord[4]][coord[5]]=4;locateC6.push_back(x);
					locateC6.push_back(y);
					locateC6.push_back(coord[2]);
					locateC6.push_back(coord[3]);
					locateC6.push_back(coord[4]);
					locateC6.push_back(coord[5]);}
					else{surface[coord[0]][coord[1]]=2;surface[coord[4]][coord[5]]=4;locateC6.push_back(x);
					locateC6.push_back(y);
					locateC6.push_back(coord[0]);
					locateC6.push_back(coord[1]);
					locateC6.push_back(coord[4]);
					locateC6.push_back(coord[5]);}
					C6_NUM+=3;
					Empty_NUM-=3;
					
				}
				if(pc==8)
				{
					surface[x][y]=1;
					locateC6.push_back(x);
					locateC6.push_back(y);
					
					/*�ڴ�����λλ�����ѡ������λ��˵��γ�����λ���� C6��ͬ��������C6 �� tag ��ʾ��ͬʱ C6 �������� 1*/
					if(R3<=(float)1/4){surface[coord[0]][coord[1]]=2;surface[coord[2]][coord[3]]=4;locateC6.push_back(coord[0]);
					locateC6.push_back(coord[1]);
					locateC6.push_back(coord[2]);
					locateC6.push_back(coord[3]);}
					else if(R3<=(float)2/4){surface[coord[2]][coord[3]]=2;surface[coord[4]][coord[5]]=4;locateC6.push_back(coord[2]);
					locateC6.push_back(coord[3]);
					locateC6.push_back(coord[4]);
					locateC6.push_back(coord[5]);}
					else if(R3<=(float)3/4){surface[coord[4]][coord[5]]=2;surface[coord[6]][coord[7]]=4;locateC6.push_back(coord[4]);
					locateC6.push_back(coord[5]);
					locateC6.push_back(coord[6]);
					locateC6.push_back(coord[7]);}
					else{surface[coord[0]][coord[1]]=2;surface[coord[6]][coord[7]]=4;locateC6.push_back(coord[0]);
					locateC6.push_back(coord[1]);
					locateC6.push_back(coord[6]);
					locateC6.push_back(coord[7]);}
					C6_NUM+=3;
					Empty_NUM-=3;
				}
			}
			/*3������С������ C2 �Ѹ�*/
			//else if(R1<=paa+pba+pag)
			else if(R1<=paa+pba+pag+pad)
			{
				if( C2_NUM<1 ) continue;
				t = rand() % (C2_NUM);
				t *= 2;
				x = locateC2[t];
				y = locateC2[t+1];
				surface[x][y] = 0;
				locateC2.erase(locateC2.begin()+t);
				locateC2.erase(locateC2.begin()+t);
				C2_NUM-=1;
				Empty_NUM+=1;
			}
			/*4������С������ C2 ����*/
		       /* else if(R1<=paa+pba+pag+pad)*/
			//{
				//if( C2_NUM<1 ) continue;
				//t = rand() % (C2_NUM);
					//t *= 2;
					//x = locateC2[t];
					//y = locateC2[t+1];
					//surface[x][y] = 0;
					//locateC2.erase(locateC2.begin()+t);
					//locateC2.erase(locateC2.begin()+t);
				//C2_NUM-=1;
				//Empty_NUM+=1;
			/*}*/
			/*5��С���ӷ����� C6 �Ѹ�*/
			//else if(R1<=paa+pba+pag+pad+pbg)
			else if(R1<=paa+pba+pag+pad+pbg+pbd)
			{
				if( C6_NUM<3 ) continue;
				t = rand() % (C6_NUM/3);
					t *= 6;
					x = locateC6[t];
					y = locateC6[t+1];
					surface[x][y] = 0;
					locateC6.erase(locateC6.begin()+t);
					locateC6.erase(locateC6.begin()+t);
					x = locateC6[t];
					y = locateC6[t+1];
					surface[x][y] = 0;
					locateC6.erase(locateC6.begin()+t);
					locateC6.erase(locateC6.begin()+t);
			
					x = locateC6[t];
					y = locateC6[t+1];
					surface[x][y] = 0;
					locateC6.erase(locateC6.begin()+t);
					locateC6.erase(locateC6.begin()+t);
				C6_NUM-=3;
				Empty_NUM+=3;
			}
			/*6��С���ӷ����� C6 ����*/
		       /* else if(R1<=paa+pba+pag+pad+pbg+pbd)*/
			//{
				//if( C6_NUM<3 ) continue;
			//t = rand() % (C6_NUM/3);
					//t *= 6;
					//x = locateC6[t];
					//y = locateC6[t+1];
					//surface[x][y] = 0;
					//locateC6.erase(locateC6.begin()+t);
					//locateC6.erase(locateC6.begin()+t);
					//x = locateC6[t];
					//y = locateC6[t+1];
					//surface[x][y] = 0;
					//locateC6.erase(locateC6.begin()+t);
					//locateC6.erase(locateC6.begin()+t);
					//x = locateC6[t];
					//y = locateC6[t+1];
					//surface[x][y] = 0;
					//locateC6.erase(locateC6.begin()+t);
					//locateC6.erase(locateC6.begin()+t);
				//C6_NUM-=3;
				//Empty_NUM+=3;
			/*}*/
			/*7������С������ C2 ����ΧС���ӷ����� C6 ��Ӧ*/
			else
			{
				if( C2_NUM<1 || C6_NUM<3 ) continue;
				if(canC2_C6_React(C2_NUM,C6_NUM)){
					sum1++;
					C2_NUM-=1;
					C6_NUM-=3;
					Empty_NUM+=4;
					}
				
			}
		}
		//finish=clock();
		//cout << n << endl;
		//cout<<"time:"<<(double)(finish - start) / CLOCKS_PER_SEC<<endl;
	}

	//Empty_NUM=0,C2_NUM=0,C6_NUM=0;
	//for(i=0;i<MAT;i++)
	//{
		//for(j=0;j<MAT;j++)
		//{
			//if(surf[i][j]==0) {Empty_NUM=Empty_NUM+1;}  /*ͳ�ƻ�������λ��*/
			//else if(surf[i][j]==-1) {C2_NUM=C2_NUM+1;}/*ͳ�ƻ����������С������ C2 ��*/
		//	else {C6_NUM=C6_NUM+1;}              /*ͳ�ƻ������С���ӷ����� C6 ��*/
		//}
	//}

	//printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",C2_NUM,C6_NUM,Empty_NUM);
	cout << "C2: " << C2_NUM << "C6: " << C6_NUM << "Empty: " << Empty_NUM << endl;
}

int main()
{
	locateC6.reserve(20000);
	locateC2.reserve(10000);
	double data[]={0.1,0.1,0.1,1,0.1,15,0.3,0.1,0.3,4.2,0.3,4.4,0.3,15,0.5,0.1,0.5,11.6,0.5,11.8,0.5,15};
   	double C2_CON,C6_CON;

    //clock_t start, finish, whole_start, whole_end;
    time_t start,finish;
    printf("TEST1\n");
    //whole_start = clock();
    start = time(NULL);
    for(int i=0; i<1; i+=2)
    {
	//cout << "llll" << endl;
        C2_CON = data[i];
        C6_CON = data[i]*data[i+1];
        cout<<i/2+1<<":"<<"C2_CON="<<C2_CON<<" R="<<data[i+1]<<endl;
        //printf("%d: C2_CON=%lf R=%lf\n",i/2+1,C2_CON,data[i+1] );
        test1(C2_CON,C6_CON);
	currentC2=0,currentC6=0;
	cC6=0,cC2=0;
	locateC2.clear();
	locateC6.clear();
	memset(surface,0,sizeof(int) * MAT_SIZE);
    }
    finish = time(NULL);
    cout << "time: " << finish-start << "s" << endl;
    //whole_end = clock();
    // ���������ʱ
    //cout<<"####################"<<endl<<"Whole using time:"<<(double)(whole_end - whole_start) / CLOCKS_PER_SEC<<endl;
	cout<<"####################"<<endl<<endl;
}
