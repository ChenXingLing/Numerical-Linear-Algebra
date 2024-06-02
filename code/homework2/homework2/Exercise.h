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

void exercise_1(){
	puts("————————————————————\n【ecercise (1)】");
	cout<<"【 5~20 阶 Hilbert 矩阵的无穷范数条件数】"<<endl;
	for(int n=5;n<=20;++n){
		Matrix A(n,n);A.Hilbert();
		cout<<"【 n = "<<setw(2)<<n<<" 】  "<<setw(11)<<Cond_inf(A)<<endl;
	}
}

void exercise_2(){
	puts("————————————————————\n【ecercise (2)】");

	for(int n=5;n<=30;++n){
		Matrix A(n,n),b(n,1),xr(n,1);A=0,b=0;
		for(int i=0;i<n;++i){
			A(i,i)=A(i,n-1)=1;
			for(int j=0;j<i;++j)A(i,j)=-1;
			xr[i]=Rand(-100,100)/10.0;
		}
		b=A*xr; Matrix x=Solve_CPLU(A,b);
		cout<<"【 n = "<<setw(2)<<n<<" 】";
		cout<<"  估计相对误差： "<<setw(11)<<Cond_inf(A)*Norm_inf(b-A*x)/Norm_inf(b)<<"   | ";
		cout<<"  真实相对误差："<<setw(11)<<Norm_inf(xr-x)/Norm_inf(xr)<<"   | ";
		cout<<"  约为 "<<setw(7)<<Cond_inf(A)*Norm_inf(b-A*x)/Norm_inf(b)/Norm_inf(xr-x)*Norm_inf(xr)<<" 倍"<<endl;
	}
}