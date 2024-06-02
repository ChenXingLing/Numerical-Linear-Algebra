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

#define Vec vector<double>
#define Vec_C vector<complex<double>>
void exercise_1(){
	puts("————————————————————\n【ecercise_1】");

	double eps=1e-6;

	puts("\n(1).[x^3+x^2-5x+3=0]");
	auto start=std::chrono::high_resolution_clock::now();
	cout<<"  模最大根: "<<MaxRoot(Vec{1,1,-5,3},eps)<<endl;
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	cout<<"  Time: "<<(duration)<<" us"<<endl;

	puts("\n(2).[x^3-3x-1=0]");
	start=std::chrono::high_resolution_clock::now();
	cout<<"  模最大根: "<<MaxRoot(Vec{1,0,-3,-1},eps)<<endl;
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	cout<<"  Time: "<<(duration)<<" us"<<endl;

	puts("\n(3).[x^8+101x^7+208.01x^6+10891.01x^5+9802.08x^4+79108.9x^3-99902x^2+790x-1000=0]");
	start=std::chrono::high_resolution_clock::now();
	cout<<"  模最大根: "<<MaxRoot(Vec{1,101,208.01,10891.01,9802.08,79108.9,-99902,790,-1000},eps)<<endl;
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	cout<<"  Time: "<<(duration)<<" us"<<endl;
}

void exercise_2_2(){
	puts("————————————————————\n【ecercise_2(2)】");
	
	double eps=1e-6;

	puts("\n(1).[x^{41}+x^3+1=0]");
	Vec a(42,0);a[0]=a[38]=a[41]=1;
	auto start=std::chrono::high_resolution_clock::now();
	vector<CP> x=ImplicitQR(Friend(a));
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.size()){
		cout<<"  全部根为: "<<endl;
		for(int i=0;i<x.size();++i)
			if(x[i].imag()==0)cout<<x[i].real()<<endl;
			else cout<<x[i].real()<<(x[i].imag()>0?" +":" ")<<x[i].imag()<<"i"<<endl;
		cout<<"  Time: "<<(duration)<<" us"<<endl;
	}
	else cout<<"  求解失败！！！"<<endl<<endl;
}

void exercise_2_3(){
	puts("————————————————————\n【ecercise_2(3)】");
	
	double eps=1e-6;

	Matrix A=Matrix(vector<Vec>{Vec{9.1,3.0,2.6,4.0},Vec{4.2,5.3,4.7,1.6},Vec{3.2,1.7,9.4,0},Vec{6.1,4.9,3.5,6.2}});

	puts("\n(1).[x=0.9]");
	A(2,3)=0.9;
	auto start=std::chrono::high_resolution_clock::now();
	vector<CP> x=ImplicitQR(A);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.size()){
		cout<<"  全部特征值为: "<<endl;
		for(int i=0;i<x.size();++i)
			if(x[i].imag()==0)cout<<x[i].real()<<endl;
			else cout<<x[i].real()<<(x[i].imag()>0?" +":" ")<<x[i].imag()<<"i"<<endl;
		cout<<"  Time: "<<(duration)<<" us"<<endl;
	}
	else cout<<"  求解失败！！！"<<endl<<endl;

	puts("\n(2).[x=1.0]");
	A(2,3)=1.0;
	start=std::chrono::high_resolution_clock::now();
	x=ImplicitQR(A);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.size()){
		cout<<"  全部特征值为: "<<endl;
		for(int i=0;i<x.size();++i)
			if(x[i].imag()==0)cout<<x[i].real()<<endl;
			else cout<<x[i].real()<<(x[i].imag()>0?" +":" ")<<x[i].imag()<<"i"<<endl;
		cout<<"  Time: "<<(duration)<<" us"<<endl;
	}
	else cout<<"  求解失败！！！"<<endl<<endl;

	puts("\n(3).[x=1.1]");
	A(2,3)=1.1;
	start=std::chrono::high_resolution_clock::now();
	x=ImplicitQR(A);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.size()){
		cout<<"  全部特征值为: "<<endl;
		for(int i=0;i<x.size();++i)
			if(x[i].imag()==0)cout<<x[i].real()<<endl;
			else cout<<x[i].real()<<(x[i].imag()>0?" +":" ")<<x[i].imag()<<"i"<<endl;
		cout<<"  Time: "<<(duration)<<" us"<<endl;
	}
	else cout<<"  求解失败！！！"<<endl<<endl;
}