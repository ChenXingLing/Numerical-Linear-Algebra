#pragma once
#include <iostream>
#include <vector>
#include <Windows.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <chrono>
#include "Function.h"

using std::vector;
using std::cout;
using std::endl;

inline int Rand(int L,int R){
	int len=R-L+1;
	int a=rand()*rand()%len+1;
	return L+a-1;
}

int is_print=1;//是否输出解向量x
int is_print_find=0;//是否输出三分过程

inline double find(Matrix &A,Matrix &b,double &eps){
	puts("三分find中..."); if(!is_print_find)is_print_SOR=0;
	double _eps=1e-3,l=_eps*3,r=2-_eps*3; 
	while(r-l>_eps){
		double mid=(l+r)/2,f1,f2;
		//if(is_print_find)
			printf("     check(mid=%lf)\n",mid);
		GS_SOR(A,b,eps,mid-_eps),f1=SOR_cnt;
		GS_SOR(A,b,eps,mid+_eps),f2=SOR_cnt;
		if(f1==f2){
			GS_SOR(A,b,eps,mid-_eps*2),f1=SOR_cnt;
			GS_SOR(A,b,eps,mid+_eps*2),f2=SOR_cnt;
		}
		if(f1>f2)l=mid;
		else r=mid;
	}
	if(!is_print_find)is_print_SOR=1;
	return l;
}
double calc1(double x,double ep,double a){
	return (1-a)*(1-exp(-x/ep))/(1-exp(-1/ep))+a*x;
}
void sakura(double ep){
	int n=100;double a=0.5,h=1.0/n;
	Matrix A(n-1,n-1),b(n-1,1);A=0;b=0;
	for(int i=0;i<n-1;++i){
		A(i,i)=-2*ep-h;
		if(i<n-2)A(i+1,i)=ep,A(i,i+1)=ep+h;
	}
	b=a*h*h;b[n-2]-=ep+h;
	Matrix ans(n-1,1);
	for(int i=0;i<n-1;++i)ans[i]=calc1((i+1.0)/n,ep,a);
	double eps=1e-6;
	printf("ε=%.4lf a=%.2lf n=%d\n",ep,a,n);
	
	puts("\n(1).[Jacobi迭代法]");
	auto start=std::chrono::high_resolution_clock::now();
	Matrix x=Jacobi(A,b,eps);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差 |Ax-b|: "<<Norm_inf(A*x-b)<<endl;
		cout<<"误差 |x-ans|: "<<Norm_inf(x-ans)<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl<<endl;

	puts("\n(2).[G-S迭代法]");
	start=std::chrono::high_resolution_clock::now();
	x=GS_SOR(A,b,eps);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差 |Ax-b|: "<<Norm_inf(A*x-b)<<endl;
		cout<<"误差 |x-ans|: "<<Norm_inf(x-ans)<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl<<endl;

	double w=find(A,b,eps);
	printf("\n(3).[SOR迭代法]  (最佳松弛因子:%lf)\n",w);
	start=std::chrono::high_resolution_clock::now();
	x=GS_SOR(A,b,eps,w);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差 |Ax-b|: "<<Norm_inf(A*x-b)<<endl;
		cout<<"误差 |x-ans|: "<<Norm_inf(x-ans)<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl<<endl;

}
void exercise_1(){
	puts("————————————————————\n【ecercise_1(1)】");
	sakura(1);
	puts("————————————————————\n【ecercise_1(2)】");
	sakura(0.1);
	puts("————————————————————\n【ecercise_1(3)】");
	sakura(0.01);
	puts("————————————————————\n【ecercise_1(4)】");
	sakura(0.0001);
	return;
}