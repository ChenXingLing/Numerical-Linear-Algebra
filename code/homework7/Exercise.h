#pragma once
#include <algorithm>
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

int is_print_A=1;//是否输出A
int is_print_Q=1;//是否输出Q

#define Vec vector<double>
void exercise_1_(int id,int n,int ed=10){

	Matrix A(n,n);A=0;
	for(int i=0;i<n;++i){
		A(i,i)=4;
		if(i<n-1)A(i+1,i)=A(i,i+1)=1;
	}

	printf("\n(%d).[n=%d, 扫描次数=%d]\n",id,n,ed);
	auto start=std::chrono::high_resolution_clock::now();
	Matrix Q=Jacobi_Threshold(A,ed);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	cout<<"  Time: "<<(duration)<<" us"<<endl;
	Vec ans;ans.clear();
	for(int i=0;i<n;++i)ans.push_back(A(i,i));
	std::sort(ans.begin(),ans.end());
	cout<<"  全部特征值："<<endl;
	for(int i=0;i<ans.size();++i)cout<<setw(9)<<ans[i]<<" ";puts("");
	if(is_print_A&&id==1)puts("  Ak:"),A.print();
	if(is_print_Q&&id==1)puts("  Qk:"),Q.print();
}
void exercise_1(){
	puts("————————————————————\n【ecercise_1】");
	for(int i=1;i<=6;++i)exercise_1_(i,40+10*i);
}

void exercise_2(){
	puts("————————————————————\n【ecercise_2】");
	int n=100;Matrix A(n,n),x(n,1),y(n,1);A=0,x=2,y=-1;y[0]=0;
	for(int i=0;i<n;++i){
		A(i,i)=2;
		if(i<n-1)A(i+1,i)=A(i,i+1)=-1;
	}

	puts("\n(1).[最小特征值]");
	auto start=std::chrono::high_resolution_clock::now();
	double lambda=Bisection(1,x,y);
	Matrix v=Inverse_Power(A,lambda);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	cout<<"  Time: "<<(duration)<<" us"<<endl;
	cout<<"  最小特征值："<<lambda<<endl;
	cout<<"  对应特征向量："<<endl;v.Trans().print();

	puts("\n(2).[最大特征值]");
	start=std::chrono::high_resolution_clock::now();
	lambda=Bisection(n,x,y);
	v=Inverse_Power(A,lambda);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	cout<<"  Time: "<<(duration)<<" us"<<endl;
	cout<<"  最大特征值："<<lambda<<endl;
	cout<<"  对应特征向量："<<endl;v.Trans().print();
}