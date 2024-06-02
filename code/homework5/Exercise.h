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

inline double f(double x,double y){return sin(x*y);}
inline double phi(double x,double y){return x*x+y*y;}
#define Poi(x,y) ((x)*(n-1)+(y))
inline void Conj(Matrix A,Matrix b,double eps){
	puts("[实用共轭梯度法]");
	auto start=std::chrono::high_resolution_clock::now();
	Matrix x=ConjGrad(A,b,eps);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差 |Ax-b|: "<<Norm_inf(A*x-b)<<endl;
	}
	else cout<<"求解失败！！！"<<endl<<endl;
}
void exercise_1(){
	puts("————————————————————\n【ecercise_1】");
	
	int n=20;double h=1.0/n;double eps=1e-7;
	Matrix A((n-1)*(n-1),(n-1)*(n-1));A=0;
	for(int i=0;i<(n-1)*(n-1);++i){
		A(i,i)=1+h*h/4;
		if((i+1)%(n-1))A(i+1,i)=A(i,i+1)=-1/4.0;
	}
	for(int i=0;i<(n-2)*(n-1);++i)A(i+n-1,i)=A(i,i+n-1)=-1/4.0;
	Matrix b((n-1)*(n-1),1);b=0;
	for(int i=0;i<n-1;++i)
		for(int j=0;j<n-1;++j)
			b[Poi(i,j)]=h*h*f((i+1)*h,(j+1)*h)/4;
	for(int j=0;j<n-1;++j)
		b[Poi(0,j)]+=phi((0+1 -1)*h,(j+1)*h)/4,
		b[Poi(n-2,j)]+=phi((n-2+1 +1)*h,(j+1)*h)/4;
	for(int i=0;i<n-1;++i)
		b[Poi(i,0)]+=phi((i+1)*h,(0+1 -1)*h)/4,
		b[Poi(i,n-2)]+=phi((i+1)*h,(n-2+1 +1)*h)/4;
	puts("");
	Conj(A,b,eps);

	return;
}
inline void sakura(int n){
	double eps=1e-7;
	Matrix A(n,n),b(n,1);A.Hilbert();b=0;
	for(int i=0;i<n;b[i]/=3,++i)
		for(int j=0;j<n;++j)
			b[i]+=A(i,j);
	printf("\n[ %d 阶Hilbert]\n",n);
	Conj(A,b,eps);
}
void exercise_2(){
	puts("————————————————————\n【ecercise_2】");
	sakura(20),sakura(40),sakura(60),sakura(80);
}
#define Vec vector<double>
void exercise_3(){
	puts("————————————————————\n【ecercise_3】");
	Matrix A(vector<Vec >{Vec{10,1,2,3,4},Vec{1,9,-1,2,-3},Vec{2,-1,7,3,-5},Vec{3,2,3,12,-1},Vec{4,-3,-5,-1,15}});
	Matrix b(vector<Vec >{Vec{12},Vec{-27},Vec{14},Vec{-17},Vec{12}});
	double eps=1e-7;
	puts("\n(1).[Jacobi迭代法]");
	auto start=std::chrono::high_resolution_clock::now();
	Matrix x=Jacobi(A,b,eps);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差 |Ax-b|: "<<Norm_inf(A*x-b)<<endl;
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
	}
	else cout<<"求解失败！！！"<<endl<<endl;

	printf("\n(3).");Conj(A,b,eps);
}