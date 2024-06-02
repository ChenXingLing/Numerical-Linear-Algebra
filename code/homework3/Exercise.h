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

void LU(Matrix A,Matrix b){//三种LU分解法
	puts("\n(1).[LU分解法]");
	auto start=std::chrono::high_resolution_clock::now();
	Matrix x=Solve_LU(A,b);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差(Ax-b中所有元的平方和): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl;

	puts("\n(2).[全主元LU分解法]");
	start=std::chrono::high_resolution_clock::now();
	x=Solve_FPLU(A,b);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差(Ax-b中所有元的平方和): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl;

	puts("\n(3).[列主元LU分解法]");
	start=std::chrono::high_resolution_clock::now();
	x=Solve_CPLU(A,b);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差(Ax-b中所有元的平方和): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl<<endl;
}

void Cholesky(Matrix A,Matrix b){//两种Cholesky分解法
	puts("\n(1).[Cholesky分解]");
	auto start=std::chrono::high_resolution_clock::now();
	Matrix x=Solve_LL(A,b);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差(Ax-b中所有元的平方和): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl;

	puts("\n(2).[Cholesky+分解]");
	start=std::chrono::high_resolution_clock::now();
	x=Solve_LDL(A,b);
	end=std::chrono::high_resolution_clock::now();
	duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差(Ax-b中所有元的平方和): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl<<endl;
}

void QR(Matrix A,Matrix b){//QR分解法
	puts("\n(1).[QR分解法]");
	auto start=std::chrono::high_resolution_clock::now();
	Matrix x=Solve_QR(A,b);
	auto end=std::chrono::high_resolution_clock::now();
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if(x.nrow()){
		cout<<"解为: X^T="<<endl; if(is_print)(x.Trans()).print();
		cout<<"运行时间: "<<(duration)<<" us"<<endl;
		cout<<"误差(Ax-b中所有元的平方和): "<<(A*x-b).sum2()<<endl<<endl;
	}
	else cout<<"求解失败！！！"<<endl<<endl;
}



void exercise_1_1(){
	puts("————————————————————\n【ecercise_1(1)】");

	int n=84;
	Matrix A(n,n),b(n,1);A=0,b=0;
	for(int i=0;i<n;++i){
		A(i,i)=6;
		if(i<n-1)A(i,i+1)=1,A(i+1,i)=8;
	}
	b[0]=7,b[n-1]=14;
	for(int i=1;i<n-1;++i)b[i]=15;

	QR(A,b);

}

void exercise_1_2(){
	puts("————————————————————\n【ecercise_1(2)】");

	int n=100;
	Matrix A(n,n),b(n,1);A=0,b=0;
	for(int i=0;i<n;++i){
		A(i,i)=10;
		if(i<n-1)A(i,i+1)=A(i+1,i)=1;
	}
	for(int i=0;i<n;++i)b[i]=Rand(1,100);
	cout<<"随机得: b^T="<<endl; (b.Trans()).print();

	LU(A,b);
	Cholesky(A,b);
	QR(A,b);
}

void exercise_1_3(){
	puts("————————————————————\n【ecercise_1(3)】");

	int n=40;
	Matrix A(n,n),b(n,1);A=0,b=0;
	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
			b[i]+=(A(i,j)=1.0/(i+j+1));

	LU(A,b);
	Cholesky(A,b);
	QR(A,b);
}

void exercise_2(){
	puts("————————————————————\n【ecercise_2】");
	vector<double>t={-1,-0.75,-0.5,0,0.25,0.5,0.75};
	vector<double>y={1,0.8125,0.75,1,1.3125,1.75,2.3125};
	int m=7,n=3;
	Matrix A(m,n),b(m,1);
	for(int i=0;i<m;++i){
		for(int j=0;j<n;++j)A(i,j)=pow(t[i],n-1-j);
		b[i]=y[i];
	}
	Matrix x=Solve_QR(A,b);
	cout<<"拟合多项式系数: "; (x.Trans()).print();
	cout<<"最小残向量 2 范数: "<<sqrt((A*x-b).sum2())<<endl<<endl;
}

void exercise_3(){
	puts("————————————————————\n【ecercise_3】");
	vector<vector<double>> _A =
		{ {1,4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0},
		{1,5.0208, 1, 3.531, 1.5, 2, 7, 4, 62, 1, 1, 0},
		{1,4.5429, 1, 2.275, 1.175, 1, 6, 3, 40,  2, 1, 0},
		{1,4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1, 0},
		{1,5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0},
		{1,3.891, 1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0},
		{1,5.898, 1, 5.85, 1.24, 1, 7, 3, 51, 2, 1,  1},
		{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
		{1,15.4202, 2.5,  9.8, 3.42, 2, 10, 5, 42, 2, 1, 1},
		{1,14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1},
		{1,5.8282, 1, 6.435, 1.225, 2, 6, 3, 32, 1, 1, 0},
		{1,5.3003, 1, 4.9883, 1.552, 1, 6, 3, 30, 1, 2, 0},
		{1,6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2, 0},
		{1,5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0},
		{1,5.05, 1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1},
		{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
		{1,8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4, 1, 0},
		{1,6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1},
		{1,7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0},
		{1,9.0384, 1, 7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0},
		{1,5.9894, 1, 5.52, 1.256, 2, 6, 3, 40, 4, 1, 1},
		{1,7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1, 0},
		{1,8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1},
		{1,6.0931, 1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0},
		{1,8.3607, 1.5, 9.15, 1.777, 2., 8, 4, 48, 1, 1, 1},
		{1,8.14, 1, 8, 1.504, 2, 7, 3, 3, 1, 3, 0},
		{1,9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1, 0},
		{1,12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
	vector<double> _b =
		{ 25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
		28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
		30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
		37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0 };
	int m=28,n=12;
	Matrix A(m,n),b(m,1);
	for(int i=0;i<m;++i){
		for(int j=0;j<n;++j)A(i,j)=_A[i][j];
		b[i]=_b[i];
	}
	Matrix x=Solve_QR(A,b);
	cout<<"拟合多项式系数: "; (x.Trans()).print();
	cout<<"最小残向量 2 范数: "<<sqrt((A*x-b).sum2())<<endl<<endl;
}