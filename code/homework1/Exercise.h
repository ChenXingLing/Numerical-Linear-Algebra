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

void LU(Matrix A,Matrix b){//����LU�ֽⷨ
	puts("\n(1).[LU�ֽⷨ]");
	auto start=std::chrono::high_resolution_clock::now();
	Matrix x=Solve_LU(A,b);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"��Ϊ: X^T="<<endl; (x.Trans()).print();
		cout<<"����ʱ��: "<<(duration)<<" us"<<endl;
		cout<<"���(Ax-b������Ԫ��ƽ����): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"���ʧ�ܣ�����"<<endl;

	puts("\n(2).[ȫ��ԪLU�ֽⷨ]");
	start=std::chrono::high_resolution_clock::now();
	x=Solve_FPLU(A,b);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"��Ϊ: X^T="<<endl; (x.Trans()).print();
		cout<<"����ʱ��: "<<(duration)<<" us"<<endl;
		cout<<"���(Ax-b������Ԫ��ƽ����): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"���ʧ�ܣ�����"<<endl;

	puts("\n(3).[����ԪLU�ֽⷨ]");
	start=std::chrono::high_resolution_clock::now();
	x=Solve_CPLU(A,b);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"��Ϊ: X^T="<<endl; (x.Trans()).print();
		cout<<"����ʱ��: "<<(duration)<<" us"<<endl;
		cout<<"���(Ax-b������Ԫ��ƽ����): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"���ʧ�ܣ�����"<<endl;
}

void Cholesky(Matrix A,Matrix b){//����Cholesky�ֽⷨ
	puts("\n(1).[Cholesky�ֽ�]");
	auto start=std::chrono::high_resolution_clock::now();
	Matrix x=Solve_LL(A,b);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"��Ϊ: X^T="<<endl; (x.Trans()).print();
		cout<<"����ʱ��: "<<(duration)<<" us"<<endl;
		cout<<"���(Ax-b������Ԫ��ƽ����): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"���ʧ�ܣ�����"<<endl;

	puts("\n(2).[Cholesky+�ֽ�]");
	start=std::chrono::high_resolution_clock::now();
	x=Solve_LDL(A,b);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"��Ϊ: X^T="<<endl; (x.Trans()).print();
		cout<<"����ʱ��: "<<(duration)<<" us"<<endl;
		cout<<"���(Ax-b������Ԫ��ƽ����): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"���ʧ�ܣ�����"<<endl;
}


void exercise_1(){
	puts("����������������������������������������\n��ecercise_1��");

	int n=84;
	Matrix A(n,n),b(n,1);A=0,b=0;
	for(int i=0;i<n;++i){
		A(i,i)=6;
		if(i<n-1)A(i,i+1)=1,A(i+1,i)=8;
	}
	b[0]=7,b[n-1]=14;
	for(int i=1;i<n-1;++i)b[i]=15;

	LU(A,b);
}

void exercise_23_1(){
	puts("����������������������������������������\n��ecercise_2(1)��");

	int n=100;
	Matrix A(n,n),b(n,1);A=0,b=0;
	for(int i=0;i<n;++i){
		A(i,i)=10;
		if(i<n-1)A(i,i+1)=A(i+1,i)=1;
	}
	for(int i=0;i<n;++i)b[i]=Rand(1,100);
	cout<<"�����: b^T="<<endl; (b.Trans()).print();

	Cholesky(A,b);

	puts("����������������������������������������\n��ecercise_3(1)��");

	LU(A,b);
}
void exercise_23_2(){
	puts("����������������������������������������\n��ecercise_2(2)��");

	int n=40;
	Matrix A(n,n),b(n,1);A=0,b=0;
	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
			b[i]+=(A(i,j)=1.0/(i+j+1));

	Cholesky(A,b);

	puts("����������������������������������������\n��ecercise_3(2)��");

	LU(A,b);
}