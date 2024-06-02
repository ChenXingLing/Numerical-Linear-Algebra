#pragma once
#include <iostream>
#include <complex>
#include <iomanip>
#include <vector>
#include <cmath>
using std::setprecision;
using std::complex;
using std::vector;
using std::cout;
using std::setw;
using std::pair;
#define mp std::make_pair

template<typename T>
class _Matrix {
public:
	_Matrix() {}
	_Matrix(int N, int M) {new_(N,M);}
	_Matrix(vector<vector<T>> _A){
		new_(_A.size(),(_A.begin())->size());
		auto I=_A.begin();
		for(int i=0;i<rows;++i,++I){
			auto J=I->begin();
			for(int j=0;j<cols;++j,++J)
				data[i*cols+j]=*J;
		}
	}
	_Matrix(vector<T> _A){
		new_(_A.size(),1);
		auto I=_A.begin();
		for(int i=0;i<rows;++i,++I){
			data[i]=*I;
		}
	}
	~_Matrix() {delete_();}

	void new_(int N,int M){data=new T[(rows=N)*(cols=M)];}
	void delete_(){if(data!=NULL){delete[] data;data=NULL;}rows=cols=0;}

	int nrow() const {return rows;}
	int ncol() const {return cols;}

	T a(int i,int j)const {//超出的部分用0补全
		return (i<rows&&j<cols)?data[i*cols+j]:0;
	}

	T& operator()(int i, int j) {return data[i*cols+j];}
	T& operator[](int k) {return data[k];}

	_Matrix col(int j) const{
		_Matrix ans(rows,1);
		for(int i=0;i<rows;++i)ans.data[i]=a(i,j);
		return ans;
	}
	_Matrix row(int i) const{
		_Matrix ans(1,cols);
		for(int j=0;j<cols;++j)ans.data[j]=a(i,j);
		return ans;
	}
	_Matrix submat(int a, int b, int c, int d) const {//(a~c,b~d)
		a--,b--,c--,d--;
		_Matrix ans(c-a+1,d-b+1);
		for(int i=0;i<ans.rows;++i)
			for(int j=0;j<ans.cols;++j)
				ans.data[i*ans.cols+j]=data[(a+i)*cols+b+j];
		return ans;
	}

	_Matrix(const _Matrix& rhs){//赋值构造函数
		new_(rhs.rows,rhs.cols);
		for(int i=0;i<rows*cols;++i)data[i]=rhs.data[i];
	}
	_Matrix& operator= (const _Matrix& rhs) {
		delete_(),new_(rhs.rows,rhs.cols);
		for(int i=0;i<rows*cols;++i)data[i]=rhs.data[i];
		return *this;
	}
	_Matrix& operator= (T v) {
		for(int i=0;i<rows*cols;++i)data[i]=v;
		return *this;
	}

	_Matrix operator+ (const _Matrix& rhs) {
		_Matrix ans(Max(rows,rhs.rows),Max(cols,rhs.cols));
		for(int i=0;i<ans.rows;++i)
			for(int j=0;j<ans.cols;++j)
				ans.data[i*ans.cols+j]=a(i,j)+rhs.a(i,j);
		return ans;
	}
	_Matrix operator- (const _Matrix& rhs) {
		_Matrix ans(Max(rows,rhs.rows),Max(cols,rhs.cols));
		for(int i=0;i<ans.rows;++i)
			for(int j=0;j<ans.cols;++j)
				ans.data[i*ans.cols+j]=a(i,j)-rhs.a(i,j);
		return ans;
	}
	_Matrix operator* (const _Matrix& rhs) {
		_Matrix ans(rows,rhs.cols);
		for(int i=0;i<ans.rows;++i)
			for(int j=0;j<ans.cols;++j){
				ans.data[i*ans.cols+j]=0;
				for(int k=0;k<Min(cols,rhs.rows);++k)
					ans.data[i*ans.cols+j]+=data[i*cols+k]*rhs.data[k*rhs.cols+j];
			}
		return ans;
	}

	_Matrix& operator+= (const _Matrix& rhs) {return *this=*this+rhs;}
	_Matrix& operator-= (const _Matrix& rhs) {return *this=*this-rhs;}
	_Matrix& operator*= (const _Matrix& rhs) {return *this=*this*rhs;}

	_Matrix operator+ (T v) {
		_Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=data[i]+v;
		return ans;
	}
	_Matrix operator- (T v) {
		_Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=data[i]-v;
		return ans;
	}
	_Matrix operator* (T v) {
		_Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=data[i]*v;
		return ans;
	}
	_Matrix operator/ (T v) {
		_Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=data[i]/v;
		return ans;
	}

	_Matrix& operator+= (T v) {return *this=*this+v;}
	_Matrix& operator-= (T v) {return *this=*this-v;}
	_Matrix& operator*= (T v) {return *this=*this*v;}
	_Matrix& operator/= (T v) {return *this=*this/v;}

	void print () const {
		printf("  [Matrix size = (%d x %d)]\n", rows, cols);
		for(int i=0;i<rows;puts(""),++i)
			for(int j=0;j<cols;++j)
				//printf("%8.4lf",data[i*cols+j]);
				//cout<<setw(9)<<setprecision(4)<<data[i*cols+j]<<" ";
				cout<<setw(14)<<data[i*cols+j]<<" ";
		puts("");
	}
	
	_Matrix Identity () {//把自己填充成单位矩阵
		for(int i=0;i<rows;++i)
			for(int j=0;j<cols;++j)
				data[i*cols+j]=(i==j);
		return *this;
	}

	_Matrix Hilbert () {//把自己填充成Hilbert矩阵
		for(int i=0;i<rows;++i)
			for(int j=0;j<cols;++j)
				data[i*cols+j]=1.0/(i+j+1);
		return *this;
	}

	_Matrix Trans () const {//返回转置矩阵
		_Matrix ans(cols,rows);
		for(int i=0;i<rows;++i)
			for(int j=0;j<cols;++j)
				ans.data[j*rows+i]=data[i*cols+j];
		return ans;
	}

	_Matrix Sign () const {//提取符号
		_Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=(data[i]>0?1:(data[i]<0?-1:0));
		return ans;
	}

	T sum2(){//矩阵所有元素平方和
		T ans=0;
		for(int i=0;i<rows*cols;++i)ans+=data[i]*data[i];
		return ans;
	}
	T abs_max(){//矩阵所有元素绝对值的最大值
		T ans=0;
		for(int i=0;i<rows*cols;++i)ans=max(ans,abs(data[i]));
		return ans;
	}

	void swap_row(int x,int y){//交换第x行和第y行
		for(int j=0;j<cols;++j)swap(data[x*cols+j],data[y*cols+j]);
	}
	void swap_col(int x,int y){//交换第x列和第y列
		for(int i=0;i<rows;++i)swap(data[i*cols+x],data[i*cols+y]);
	}
	void swap_rows(vector<int>&O){
		for(int i=0;i<O.size();++i)swap_row(i,O[i]);
	}
	void swap_rows_(vector<int>&O){
		for(int i=O.size()-1;i>=0;--i)swap_row(i,O[i]);
	}

	/*分解成功返回1，否则返回0*/
	bool LU(){//LU分解
		if(rows!=cols){puts("Error!!! LU:非方阵");return 0;}
		int n=rows;
		for(int i=0;i<n-1;++i){
			if(!data[i*cols+i]){puts("Error!!! LU:除0");return 0;}
			for(int j=i+1;j<n;++j){
				data[j*cols+i]/=data[i*cols+i];
				for(int k=i+1;k<n;++k)
					data[j*cols+k]-=data[j*cols+i]*data[i*cols+k];
			}
		}
		return 1;
	}
	bool FPLU(vector<int>&P,vector<int>&Q){//全主元LU分解
		if(rows!=cols){puts("Error!!! FPLU:非方阵");return 0;}
		int n=rows;P.clear(),Q.clear();
		for(int i=0;i<n-1;++i){
			int p=i,q=i;
			for(int j=i;j<n;++j)
				for(int k=i;k<n;++k)
					if(Abs(data[j*cols+k])>Abs(data[p*cols+q]))
						p=j,q=k;
			if(!data[p*cols+q]){puts("Error!!! FPLU:除0");return 0;}
			swap_row(i,p),P.push_back(p);
			swap_col(i,q),Q.push_back(q);
			for(int j=i+1;j<n;++j){
				data[j*cols+i]/=data[i*cols+i];
				for(int k=i+1;k<n;++k)
					data[j*cols+k]-=data[j*cols+i]*data[i*cols+k];
			}
		}
		return 1;
	}
	bool CPLU(vector<int>&P){//列主元LU分解
		if(rows!=cols){puts("Error!!! CPLU:非方阵");return 0;}
		int n=rows;P.clear();
		for(int i=0;i<n-1;++i){
			int p=i;
			for(int j=i;j<n;++j)
				if(Abs(data[j*cols+i])>Abs(data[p*cols+i]))
					p=j;
			if(!data[p*cols+i]){puts("Error!!! CPLU:除0");return 0;}
			swap_row(i,p),P.push_back(p);
			for(int j=i+1;j<n;++j){
				data[j*cols+i]/=data[i*cols+i];
				for(int k=i+1;k<n;++k)
					data[j*cols+k]-=data[j*cols+i]*data[i*cols+k];
			}
		}
		return 1;
	}
	bool LL(){//LL分解
		if(rows!=cols){puts("Error!!! LL:非方阵");return 0;}
		int n=rows;
		for(int i=0;i<n;++i){
			if(data[i*cols+i]<0){puts("Error!!! LL:无法求平方根");return 0;}
			if(!data[i*cols+i]){puts("Error!!! LL:除0");return 0;}
			data[i*cols+i]=sqrt(data[i*cols+i]);
			for(int j=i+1;j<n;++j)data[j*cols+i]/=data[i*cols+i];
			for(int j=i+1;j<n;++j)
				for(int k=j;k<n;++k)
					data[k*cols+j]-=data[k*cols+i]*data[j*cols+i];
		}
		return 1;
	}
	bool LDL(){//LDL分解
		if(rows!=cols){puts("Error!!! LDL:非方阵");return 0;}
		int n=rows;vector<T>v;
		for(int j=0;j<n;++j){
			v.clear();
			for(int i=0;i<j;++i)v.push_back(data[j*cols+i]*data[i*cols+i]);
			for(int i=0;i<j;++i)
				data[j*cols+j]-=data[j*cols+i]*v[i];
			if(!data[j*cols+j]){puts("Error!!! LDL:除0");return 0;}
			for(int i=j+1;i<n;++i){
				for(int k=0;k<j;++k)
					data[i*cols+j]-=data[i*cols+k]*v[k];
				data[i*cols+j]/=data[j*cols+j];
			}
		}
		return 1;
	}
	
	/*求解1范数*/
	T Norm_1(_Matrix A){
		T ans=0;
		for(int j=0;j<A.ncol();++j){
			T sum=0;
			for(int i=0;i<A.nrow();++i)sum+=abs(A(i,j));
			ans=max(ans,sum);
		}
		return ans;
	}
	/*求解∞范数*/
	T Norm_inf(){
		return Norm_1(Trans());
	}

	_Matrix house(_Matrix x,T &beta){//计算householder变换向量
		int n=x.nrow();
		x/=x.Norm_inf();
		_Matrix v=x;v[0]=0;//v(2:n)=x(2:n)
		T sigma=v.sum2();//sigma= x(2)^2 + x(3)^2 +...+ x(n)^2
		if(sigma==0)beta=2,v[0]=1;//beta=0;//另一种写法：beta=2,v[0]=1;
		else{
			T alpha=sqrt(x[0]*x[0]+sigma);//alpha= x(1)^2 + sigma
			if(x[0]<=0)v[0]=x[0]-alpha;//v(1)= x(1)-alpha
			else v[0]=-sigma/(x[0]+alpha);//v(1)= -sigma / (x(1)+alpha)
			beta=2*v[0]*v[0]/(sigma+v[0]*v[0]);//beta= 2 * v(1)^2 / (sigma + v(1)^2)
			v/=v[0];//v= v / v(1)
		}
		return v;
	}
	_Matrix QR(){//QR分解(HouseHoder法)
		int m=rows,n=cols;_Matrix d=_Matrix(n,1);
		for(int j=0;j<n;++j){
			T beta; _Matrix v=house(submat(j+1,j+1,m,j+1),beta);
			_Matrix tmp=(v*beta) * (v.Trans()*submat(j+1,j+1,m,n));
			for(int k=j;k<m;++k)
				for(int l=j;l<n;++l)
					data[k*cols+l]-=tmp(k-j,l-j);  //A(j:m,j:n)=(I-beta* v*v^T) * A(j:m,j:n)
			d[j]=beta;
			for(int k=j+1;k<m;++k)data[k*cols+j]=v[k-j];
		}
		return d;
	}

	_Matrix Hessenberg(){//上Hessenberg化
		if(rows!=cols){puts("Error!!! Hessenberg:矩阵大小不合");return _Matrix(0,0);}
		int n=rows;
		for(int k=0;k<n-2;++k){
			T beta;_Matrix v=house(submat(k+2,k+1,n,k+1),beta),tmp;
			tmp=(v*beta) * (v.Trans()*submat(k+2,k+1,n,n));
			for(int i=k+1;i<n;++i)
				for(int j=k;j<n;++j)
					data[i*cols+j]-=tmp(i-k-1,j-k);
			tmp=(submat(1,k+2,n,n) * (v*beta)) * v.Trans();
			for(int i=0;i<n;++i)
				for(int j=k+1;j<n;++j)
					data[i*cols+j]-=tmp(i,j-k-1);
		}
		return *this;
	}
	_Matrix DoubleStepQR(){//双重步位移QR迭代
		_Matrix A=*this;int n=rows,m=n-1;
		if(n<=2)return *this;
		T s=A(m-1,m-1)+A(n-1,n-1);
		T t=A(m-1,m-1)*A(n-1,n-1)-A(m-1,n-1)*A(n-1,m-1);
		T x=A(0,0)*A(0,0)+A(0,1)*A(1,0)-s*A(0,0)+t;
		T y=A(1,0)*(A(0,0)+A(1,1)-s);
		T z=A(1,0)*A(2,1);
		for(int k=0;k<n-2;++k){
			_Matrix xyz(3,1);xyz[0]=x,xyz[1]=y,xyz[2]=z;
			T beta;_Matrix v=house(xyz,beta),tmp;
			int q=max(1,k);tmp=(v*beta) * (v.Trans()*A.submat(k+1,q,k+3,n));
			for(int i=k;i<k+3;++i)
				for(int j=q-1;j<n;++j)
					A(i,j)-=tmp(i-k,j-q+1);
			int r=min(k+4,n);tmp=(A.submat(1,k+1,r,k+3) * (v*beta)) * v.Trans();
			for(int i=0;i<r;++i)
				for(int j=k;j<k+3;++j)
					A(i,j)-=tmp(i,j-k);
			x=A(k+1,k),y=A(k+2,k);if(k<n-3)z=A(k+3,k);
		}
		_Matrix xy(2,1);xy[0]=x,xy[1]=y;
		T beta;_Matrix v=house(xy,beta),tmp;
		tmp=(v*beta) * (v.Trans()*A.submat(n-1,n-2,n,n));
		for(int i=n-1 -1;i<n;++i)
			for(int j=n-2 -1;j<n;++j)
				A(i,j)-=tmp(i-n+2,j-n+3);
		tmp=(A.submat(1,n-1,n,n) * (v*beta)) * v.Trans();
		for(int i=0;i<n;++i)
			for(int j=n-1 -1;j<n;++j)
				A(i,j)-=tmp(i,j-n+2);
		return *this=A;
	}

	inline void house_left(_Matrix &A,_Matrix &v,T &beta,int a,int b,int c,int d){
		_Matrix tmp=(v*beta) * (v.Trans()*A.submat(a,b,c,d));
		for(int i=a-1;i<c;++i)
			for(int j=b-1;j<d;++j)
				A(i,j)-=tmp(i-a+1,j-b+1);
	}
	inline void house_right(_Matrix &A,_Matrix &v,T &beta,int a,int b,int c,int d){
		_Matrix tmp=(A.submat(a,b,c,d)*v) * (v.Trans()*beta);
		for(int i=a-1;i<c;++i)
			for(int j=b-1;j<d;++j)
				A(i,j)-=tmp(i-a+1,j-b+1);
	}
	inline void givens_left(T &c,T &s,int p,int k){
		for(int j=0;j<cols;++j){
			double Apj=(*this)(p,j),Akj=(*this)(k,j);
			(*this)(p,j)=Apj*c-Akj*s,(*this)(k,j)=Apj*s+Akj*c;
		}
	}
	inline void givens_right(T &c,T &s,int p,int k){
		for(int i=0;i<cols;++i){
			double Aip=(*this)(i,p),Aik=(*this)(i,k);
			(*this)(i,p)=Aip*c-Aik*s,(*this)(i,k)=Aip*s+Aik*c;
		}
	}

	_Matrix TwoDiag(_Matrix &P,_Matrix &Q){//二对角化
		_Matrix A=*this;int m=A.nrow(),n=A.ncol();
		for(int k=0;k<n;++k){
			T beta;_Matrix v=house(A.submat(k+1,k+1,m,k+1),beta);
			house_left(A,v,beta,k+1,k+1,m,n);
			house_left(P,v,beta,k+1,1,m,m);
			if(k<n-1){
				v=house((A.submat(k+1,k+2,k+1,n)).Trans(),beta);
				house_right(A,v,beta,k+1,k+2,m,n);
				house_right(Q,v,beta,1,k+2,n,n);
			}
		}
		return *this=A;
	}
	_Matrix TwoDiag(_Matrix &P,_Matrix &Q,int p,int q){
		_Matrix A=*this;int m=A.nrow(),n=A.ncol();
		if(p>=n)p=0;if(q>=n)q=0;
		for(int k=p;k<q;++k){
			T beta;_Matrix v=house(A.submat(k+1,k+1,m,k+1),beta);
			house_left(A,v,beta,k+1,k+1,m,n);
			house_left(P,v,beta,k+1,1,m,m);
			if(k<q-1){
				v=house((A.submat(k+2,k+2,k+2,n)).Trans(),beta);
				house_right(A,v,beta,k+2,k+2,m,n);
				house_right(Q,v,beta,1,k+2,n,n);
			}
		}
		return *this=A;
	}

private:
	int rows, cols;
	T *data;
	int Max(int a,int b){return a>b?a:b;}
	int Min(int a,int b){return a<b?a:b;}
	T Abs(T x){return x<0?-x:x;}
	void swap(T &x,T &y){T c=x;x=y,y=c;}
};
#define CP complex<double>
#define Matrix _Matrix<double>
#define Matrix_C _Matrix<CP>
Matrix_C Mat_to_C(Matrix A){//实矩阵转换成复矩阵
	Matrix_C A_C(A.nrow(),A.ncol());
	for(int i=0;i<A.nrow();++i)
		for(int j=0;j<A.ncol();++j)
			A_C(i,j)=A(i,j);
	return A_C;
}
/*求解成功返回1，否则返回0*/
bool Sovle_Lxb(Matrix& L,Matrix& b,bool flag=0) {//前代法解下三角方程(flag=1时表示单位下三角)
	if(L.nrow()!=L.ncol()||L.nrow()!=b.nrow()||b.ncol()!=1){
		puts("Error!!! Sovle_Uxb:矩阵大小不合");return 0;
	}
	int n=L.nrow();
	for(int j=0;j<n;++j){
		if(!flag){//非单位
			if(!L(j,j)){puts("Error!!! Sovle_Lxb:除0");return 0;}
			b[j]/=L(j,j);
		}
		for(int k=j+1;k<n;++k)
			b[k]-=b[j]*L(k,j);
	}
	return 1;
}
bool Sovle_Uxb(Matrix& U,Matrix& b) {//回代法解上三角方程
	if(U.nrow()!=U.ncol()||U.nrow()!=b.nrow()||b.ncol()!=1){
		puts("Error!!! Sovle_Uxb:矩阵大小不合");return 0;
	}
	int n=U.nrow();
	for(int j=n-1;j>=0;--j){
		if(!U(j,j)){puts("Error!!! Sovle_Uxb:除0");printf("U(%d,%d)\n",j,j);return 0;}
		b[j]/=U(j,j);
		for(int k=j-1;k>=0;--k)
			b[k]-=b[j]*U(k,j);
	}
	return 1;
}

/*求解成功返回答案矩阵，否则返回空矩阵*/
Matrix Solve_LU(Matrix A, Matrix b) {//Guass
	if(!A.LU())return Matrix(0,0);
	if(!Sovle_Lxb(A,b,1))return Matrix(0,0);//单位下三角
	if(!Sovle_Uxb(A,b))return Matrix(0,0);//上三角
	return b;
}
Matrix Solve_FPLU(Matrix A, Matrix b) {//全主元Guass
	vector<int>P,Q;
	if(!A.FPLU(P,Q))return Matrix(0,0);
	b.swap_rows(P);
	if(!Sovle_Lxb(A,b,1))return Matrix(0,0);//单位下三角
	if(!Sovle_Uxb(A,b))return Matrix(0,0);//上三角
	b.swap_rows_(Q);
	return b;
}
Matrix Solve_CPLU(Matrix A, Matrix b) {//列主元Guass
	vector<int>P;
	if(!A.CPLU(P))return Matrix(0,0);
	b.swap_rows(P);
	if(!Sovle_Lxb(A,b,1))return Matrix(0,0);//单位下三角
	if(!Sovle_Uxb(A,b))return Matrix(0,0);//上三角
	return b;
}
Matrix Solve_LL(Matrix A,Matrix b){//Cholesky
	if(!A.LL())return Matrix(0,0);
	if(!Sovle_Lxb(A,b))return Matrix(0,0);//下三角
	A=A.Trans();
	if(!Sovle_Uxb(A,b))return Matrix(0,0);//上三角
	return b;
}
Matrix Solve_LDL(Matrix A,Matrix b){//Cholesky+
	if(!A.LDL())return Matrix(0,0);
	Matrix AT=A.Trans();int n=AT.nrow();
	for(int i=0;i<n;++i)
		for(int j=i+1;j<n;++j)
			AT(i,j)*=AT(i,i);
	if(!Sovle_Lxb(A,b,1))return Matrix(0,0);//单位下三角
	if(!Sovle_Uxb(AT,b))return Matrix(0,0);//上三角
	return b;
}

/*求解1范数*/
double Norm_1(Matrix A){
	double ans=0;
	for(int j=0;j<A.ncol();++j){
		double sum=0;
		for(int i=0;i<A.nrow();++i)sum+=abs(A(i,j));
		ans=max(ans,sum);
	}
	return ans;
}
/*求解1范数及最大值对应列*/
double Norm_1(Matrix A,int &p){
	double ans=0;p=0;
	for(int j=0;j<A.ncol();++j){
		double sum=0;
		for(int i=0;i<A.nrow();++i)sum+=abs(A(i,j));
		if(sum>ans)ans=sum,p=j;
	}
	return ans;
}
/*求解∞范数*/
double Norm_inf(Matrix A){
	return Norm_1(A.Trans());
}
/*求解∞范数及最大值对应行*/
double Norm_inf(Matrix A,int &p){
	return Norm_1(A.Trans(),p);
}
/*求解矩阵逆的∞范数*/
double Norm_invinf(Matrix A){
	int n=A.nrow(),p=0;
	Matrix x(n,1),w,v,z;x=(1.0/n);
	while(1){
		w=Solve_CPLU(A.Trans(),x);
		v=w.Sign();
		z=Solve_CPLU(A,v);
		double tmp=Norm_inf(z,p);
		if(tmp<=(z.Trans()*x)(0,0))return Norm_1(w);
		else x=0,x(p,0)=1;
	}
}
/*求解∞范数条件数*/
double Cond_inf(Matrix A){
	return Norm_invinf(A)*Norm_inf(A);
}

/*求解成功返回答案矩阵，否则返回空矩阵*/
Matrix Solve_QR(Matrix A, Matrix b){//QR分解法
	int m=A.nrow(),n=A.ncol();
	Matrix d=A.QR();
	for(int i=0;i<n;++i){//计算(Q^T *b)
		double beta=d[i];
		if(beta==0)continue;
		Matrix v(m,1);v=0;
		v[i+0]=1;//v(1)=1
		for(int j=i+1;j<m;++j)v[j]=A(j,i);
		b=b-(v*beta)*(v.Trans()*b);
	}
	Matrix R=A.submat(1,1,n,n),b_=b.submat(1,1,n,1);
	if(!Sovle_Uxb(R,b_))return Matrix(0,0);//上三角 (解 Rx = Q^T *b)
	return b_;
}

/*迭代至满足精度要求返回答案矩阵，否则返回空矩阵*/
Matrix Jacobi(Matrix A,Matrix b,double eps){//Jacobi迭代法
	printf("(Jacobi: 终止条件:eps=%.7lf)\n",eps);
	if(A.nrow()!=A.ncol()){puts("Error!!! Jacobi:矩阵大小不合");return Matrix(0,0);}
	int n=A.nrow();Matrix D(n,1);
	for(int i=0;i<n;++i)D[i]=A(i,i),A(i,i)=0;
	Matrix x(n,1),last(n,1);last=1;
	int ed=2e8/n/n; //最大迭代次数
	for(int O=0;O<ed;++O){
		x=b-A*last;
		for(int i=0;i<n;++i)x[i]/=D[i];
		//printf("  Jacobi: 迭代 %d 次，误差：%lf\n",O,Norm_inf(x-last));
		if(Norm_inf(x-last)<eps){printf("  Jacobi: 迭代 %d 次\n",O);return x;}
		last=x;
	}
	puts("Error!!! Jacobi:迭代超时");return Matrix(0,0);
}
int SOR_cnt;bool is_print_SOR=1;
Matrix GS_SOR(Matrix A,Matrix b,double eps,double w=0){//GS迭代法 / SOR迭代法
	if(w){if(is_print_SOR)printf("(SOR: 终止条件:eps=%.6lf 松弛因子:w=%.4lf)\n",eps,w);}
	else printf("(GS: 终止条件:eps=%.7lf)\n",eps);
	if(A.nrow()!=A.ncol()){puts("Error!!! GS_SOR:矩阵大小不合");return Matrix(0,0);}
	int n=A.nrow();
	Matrix x(n,1),last(n,1);last=1;
	int ed=2e8/n/n; //最大迭代次数
	for(int O=0;O<ed;++O){
		for(int i=0;i<n;++i){
			x[i]=0;
			for(int j=0;j<i;++j)x[i]-=A(i,j)*x[j];//L*x_k
			for(int j=i+1;j<n;++j)x[i]-=A(i,j)*last[j]; //U*x_{k-1}
			(x[i]+=b[i])/=A(i,i);
			if(w)x[i]=(1-w)*last[i]+w*x[i];
		}
		//printf("  Jacobi: 迭代 %d 次，误差：%lf\n",O,Norm_inf(x-last));
		if(Norm_inf(x-last)<eps){if(is_print_SOR)printf("  GS_SOR: 迭代 %d 次\n",O);SOR_cnt=O;return x;}
		last=x;
	}
	if(is_print_SOR)puts("Error!!! GS_SOR:迭代超时");SOR_cnt=2e8;return Matrix(0,0);
}

/*迭代至满足精度要求返回答案矩阵，否则返回空矩阵*/
Matrix ConjGrad(Matrix A,Matrix b,double eps){//实用共轭梯度法
	if(A.nrow()!=A.ncol()){puts("Error!!! ConjGrad:矩阵大小不合");return Matrix(0,0);}
	int n=A.nrow();Matrix x(n,1);x=1;
	Matrix p,r=b-A*x;double P=r.sum2(),P_;
	int ed=2e8/n/n; //最大迭代次数
	for(int O=1;O<ed;++O){
		p=(O==1?r:r+p*(P/P_));
		double alpha=P/(p.Trans()*A*p)[0];x+=p*alpha;
		r-=A*p*(alpha),P_=P,P=r.sum2();
		if(Norm_inf(p*alpha)<eps){printf("  ConjGrad: 迭代 %d 次\n",O);return x;}
	}
	puts("Error!!! ConjGrad:迭代超时");return Matrix(0,0);
}

/*迭代至满足精度要求返回答案矩阵，否则返回空矩阵*/
Matrix Power(Matrix A,double eps){//幂法
	if(A.nrow()!=A.ncol()){puts("Error!!! Power:矩阵大小不合");return Matrix(0,0);}
	int n=A.nrow();
	Matrix u(n,1),v(n,1);u=0,v=0;
	int ok=0;int ed=2e8/n/n/n; //最大迭代次数
	for(int i=1;i<=n;++i){
		u=0;u(i,1)=1;
		for(int O=1;O<=ed&&!ok;++O){
			int p;v=A*u;
			if(Norm_inf(v,p)<eps)break;
			v/=v[p];
			if(Norm_inf(v-u)<eps){printf("  Power: 迭代 %d 次\n",O);return v;}
			u=v;
		}
	}
	puts("Error!!! Power: 迭代超时");return Matrix(0,0);
}
Matrix Friend(vector<double> a){//根据多项式系数构造友方阵
	int n=a.size()-1;
	if(n<=0){puts("Error!!! Friend:次数过小");return Matrix(0,0);}
	if(!a[0]){puts("Error!!! Friend:首项为0");return Matrix(0,0);}
	Matrix A(n,n);A=0;
	for(int i=0;i<n;++i){
		A(i,n-1)=-a[n-i]/a[0];
		if(i<n-1)A(i+1,i)=1;
	}
	return A;
}
double MaxRoot(vector<double>a,double eps){//模最大根
	if(a.size()<=1){puts("Error!!! MaxRoot:次数过小");return 0;}
	if(!a[0]){puts("Error!!! MaxRoot:首项为0");return 0;}
	Matrix A=Friend(a),x=Power(A,eps),y=A*x;
	int p;Norm_inf(y,p);return y[p];
}
inline void vieta_root(Matrix A,CP &x1,CP &x2){//解二阶矩阵特征根
	double sum=A(0,0)+A(1,1),prod=A(0,0)*A(1,1)-A(0,1)*A(1,0);
	double tmp=sum*sum-4*prod;
	if(tmp<0)x1=CP(sum/2,sqrt(-tmp)/2),x2=CP(sum/2,-sqrt(-tmp)/2);//虚根
	else x1=(sum+sqrt(tmp))/2,x2=(sum-sqrt(tmp))/2;//实根
	return;
}
/*将特征根存入vector*/
void Schur(Matrix A,vector<CP> &ans){//Schur分解求特征值
	if(A.nrow()!=A.ncol()){puts("Error!!! Schur:矩阵大小不合");return;}
	int n=A.nrow();
	for(int i=0;i<n;++i)
		if(i==n-1||!A(i+1,i))ans.push_back(A(i,i));//1阶子矩阵
		else{//2阶子矩阵
			CP x1,x2;vieta_root(A.submat(i+1,i+1,i+2,i+2),x1,x2);
			ans.push_back(x1),ans.push_back(x2),++i;
		}
	return;
}
void ImplicitQR_(Matrix A,vector<CP> &ans,int &deep){//隐式QR递归迭代函数
	int n=A.nrow();
	if(n<=2){Schur(A,ans);return;}
	for(int i=1;i<n;++i)
		if(abs(A(i,i-1))<(abs(A(i,i))+abs(A(i-1,i-1)))*DBL_EPSILON)
			A(i,i-1)=0;
	int m=n,l=0;
	for(int i=n-1;i>=1;--i)if(A(i,i-1)!=0&&(i==1||A(i-1,i-2)!=0)){m=i+1;break;}//拟上三角 A(m:n-1) (m+1~n)
	for(int i=m-1;i>=1;--i)if(A(i,i-1)==0){l=i;break;}//不可约的上Hessenberg A(l:m-1) (l+1~m)
	if(m<n)Schur(A.submat(m+1,m+1,n,n),ans);//如果有H33
	if(l<m)ImplicitQR_((A.submat(l+1,l+1,m,m)).DoubleStepQR(),ans,++deep);//如果有H22
	if(l>0)ImplicitQR_(A.submat(1,1,l,l),ans,deep);//如果有H11
	return;
}
/*返回所有特征根*/
vector<CP> ImplicitQR(Matrix A){//隐式QR算法
	vector<CP>ans;ans.clear();
	if(A.nrow()!=A.ncol()){puts("Error!!! ImplicitQR:矩阵大小不合");return ans;}
	int cnt=0; ImplicitQR_(A.Hessenberg(),ans,cnt);
	printf("  ImplicitQR: 迭代 %d 次\n",cnt);
	return ans;
}

/*进行ed次扫描，返回变换矩阵*/
Matrix Jacobi_Threshold(Matrix &A,int ed=10){//过关Jacobi法
	if(A.nrow()!=A.ncol()){puts("Error!!! Jacobi_Threshold:矩阵大小不合");return Matrix(0,0);}
	int n=A.nrow();Matrix Q(n,n);Q.Identity();
	if(n<=1)return Q;
	double delta=0;
	for(int i=0;i<n-1;++i)
		for(int j=i+1;j<n;++j)
			delta+=A(i,j)*A(i,j);
	int O=0;
	for(int T=1;T<=ed;++T){//扫描次数
		if(delta==0)break;
		delta/=n;
		for(int i=0;i<n-1;++i)
			for(int j=i+1;j<n;++j)
				if(delta<abs(A(i,j))){
					++O;
					double tmp=(A(j,j)-A(i,i))/2/A(i,j);
					double t=(tmp>=0?1:-1)/(abs(tmp)+sqrt(1+tmp*tmp));
					double c=1/sqrt(1+t*t),s=c*t;
					for(int k=0;k<n;++k){
						double Aki=A(k,i),Akj=A(k,j);
						A(k,i)=Aki*c-Akj*s,A(k,j)=Aki*s+Akj*c;
					}
					for(int l=0;l<n;++l){
						double Ail=A(i,l),Ajl=A(j,l);
						A(i,l)=Ail*c-Ajl*s,A(j,l)=Ail*s+Ajl*c;
					}

					for(int k=0;k<n;++k){
						double Qki=Q(k,i),Qkj=Q(k,j);
						Q(k,i)=Qki*c-Qkj*s,Q(k,j)=Qki*s+Qkj*c;
					}
				}
	}
	printf("  Jacobi_Threshold: 迭代 %d 次\n",O);
	return Q.Trans();
}
int Bisection_Sn(double lambda,Matrix x,Matrix y){
	int n=x.nrow(),s=0;double q=x[0]-lambda;
	if(n==1)return x[0]<lambda;
	for(int i=0;i<n;++i){
		if(q<0)++s;
		if(i<n-1){
			if(y[i+1]==0)return s+Bisection_Sn(lambda,x.submat(i+2,1,n,1),y.submat(i+2,1,n,1));
			if(q==0)q=abs(y[i+1])*DBL_EPSILON;
			q=x[i+1]-lambda-y[i+1]*y[i+1]/q;
		}
	}
	return s;
}
/*二分迭代至满足精度要求，返回第K大特征值*/
double Bisection(int K,Matrix x,Matrix y,double eps=1e-10){//二分法
	if(x.nrow()!=y.nrow()){puts("Error!!! Bisection:矩阵大小不合");return 0;}
	int n=x.nrow();double tmp=abs(x[n-1])+abs(y[n-1]);
	for(int i=0;i<n-1;++i)
		tmp=max(tmp,abs(x[i])+abs(y[i])+abs(y[i+1]));
	double l=-tmp,r=tmp;int O=0;
	while(r-l>eps){
		double mid=(l+r)/2;++O;
		if(Bisection_Sn(mid,x,y)<K)l=mid;
		else r=mid;
	}
	printf("  Bisection: 迭代 %d 次\n",O);
	return l;
}
/*迭代至满足精度要求返回答案矩阵，否则返回空矩阵*/
Matrix Inverse_Power(Matrix A,double lambda,double eps=1e-10){//反幂法
	if(A.nrow()!=A.ncol()){puts("Error!!! Inverse_Power:矩阵大小不合");return Matrix(0,0);}
	int n=A.nrow();
	Matrix B,I(n,n),u(n,1),v(n,1),y;I.Identity(),u=1;
	B=A-I*lambda,B.LU();int ed=2e8/n/n;
	for(int O=1;O<=ed;++O){
		y=u,Sovle_Lxb(B,y,1);
		v=y,Sovle_Uxb(B,v);
		u=v/sqrt(v.sum2());
		if(Norm_inf(A*u-u*lambda)<eps){printf("  Inverse_Power: 迭代 %d 次\n",O);return u;}
	}
	puts("Error!!! Inverse_Power: 迭代超时");return Matrix(0,0);
}


inline void SVD_(Matrix &A,Matrix &P,Matrix &Q,int p,int q){
	--q;
	if(q>p+1){
		double alpha=A(q,q)*A(q,q)+A(q-1,q)*A(q-1,q);
		double delta=(A(q-1,q-1)*A(q-1,q-1)+A(q-2,q-1)*A(q-2,q-1)-alpha)/2;
		double beta=A(q-1,q-1)*A(q-1,q);
		double mu=alpha-beta*beta/(delta+(delta>=0?1:-1)*sqrt(delta*delta+beta*beta));
		double y=A(p,p)*A(p,p)-mu,z=A(p,p)*A(p,p+1);
		for(int k=p;k<q;++k){
			double c=-y/sqrt(y*y+z*z),s=z/sqrt(y*y+z*z);
			y=c*A(k,k)-s*A(k,k+1),z=-s*A(k+1,k+1);
			A.givens_right(c,s,k,k+1),Q.givens_right(c,s,k,k+1);
			c=-y/sqrt(y*y+z*z),s=z/sqrt(y*y+z*z);
			if(k<q-1)y=c*A(k,k+1)-s*A(k+1,k+1),z=-s*A(k+1,k+2);
			A.givens_left(c,s,k,k+1),P.givens_left(c,s,k,k+1);
		}
	}
	else if(q==p+1){
		double mu=(A(p,p)*A(p,p)+A(p,q)*A(p,q)-A(q,q)*A(q,q))/A(p,q)/A(q,q);
		double t=(mu+sqrt(mu*mu+4))/2;
		double c=1/sqrt(1+t*t),s=t*c;
		double s_=s*A(p,p),c_=s*A(p,q)+c*A(q,q);
		A.givens_left(c,s,p,q),P.givens_left(c,s,p,q);
		c=c_/sqrt(c_*c_+s_*s_),s=s_/sqrt(c_*c_+s_*s_);
		A.givens_right(c,s,p,q),Q.givens_right(c,s,p,q);
	}
}
/*迭代至满足要求返回D=PAQ*/
Matrix SVD(Matrix A,Matrix &P,Matrix &Q,double eps=1e-10){//SVD迭代
	int m=A.nrow(),n=A.ncol();
	if(m<n){
		Matrix D=SVD(A.Trans(),Q,P);
		P=P.Trans(),Q=Q.Trans();
		return D.Trans();
	}
	P=Matrix(m,m),Q=Matrix(n,n);P.Identity(),Q.Identity();
	A.TwoDiag(P,Q);
	int p,q=n,O=0,ed=2e8/n/m/n;//最大迭代次数
	while(O<ed){
		for(int i=0;i<n;++i){
			if(i<n-1&&abs(A(i,i+1))<(abs(A(i,i))+abs(A(i+1,i+1)))*eps)A(i,i+1)=0;
			if(abs(A(i,i))<Norm_inf(A)*eps)A(i,i)=0;
		}
		while(q>1&&A(q-2,q-1)==0)--q;
		if(q==1){
			for(int i=0;i<n;++i)if(A(i,i)<0){
				for(int j=0;j<n;++j)Q(j,i)=-Q(j,i);
				A(i,i)=-A(i,i);
			}
			break;
		}
		p=q-2;
		while(p&&A(p-1,p)!=0)--p;
		int ok=0;
		for(int i=1;i<q-p;++i)
			if(A(q-i-1,q-i-1)==0){
				for(int k=q-i+1;k<q;++k){//对应行化零
					double t=sqrt(A(k,k)*A(k,k)+A(q-i,k)*A(q-i,k));
					double c=A(k,k)/t,s=-A(q-i,k)/t;
					A.givens_left(c,s,q-i,k),P.givens_left(c,s,q-i,k);
				}
				A.TwoDiag(P,Q,p,q); ok=1;
				break;
			}
		if(!ok)SVD_(A,P,Q,p,q),++O;
		
	}
	printf("  SVD: 迭代 %d 次\n",O); return A;
}