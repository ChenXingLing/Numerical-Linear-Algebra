#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using std::setprecision;
using std::vector;
using std::cout;
using std::setw;
using std::pair;
#define mp std::make_pair

class Matrix {
public:
	Matrix() {}
	Matrix(int N, int M) {new_(N,M);}
	Matrix(vector<vector<double>> _A){
		new_(_A.size(),(_A.begin())->size());
		auto I=_A.begin();
		for(int i=0;i<rows;++i,++I){
			auto J=I->begin();
			for(int j=0;j<cols;++j,++J)
				data[i*cols+j]=*J;
		}
	}
	Matrix(vector<double> _A){
		new_(_A.size(),1);
		auto I=_A.begin();
		for(int i=0;i<rows;++i,++I){
			data[i]=*I;
		}
	}
	~Matrix() {delete_();}

	void new_(int N,int M){data=new double[(rows=N)*(cols=M)];}
	void delete_(){if(data!=NULL){delete[] data;data=NULL;}rows=cols=0;}

	int nrow() const {return rows;}
	int ncol() const {return cols;}

	double a(int i,int j)const {//超出的部分用0补全
		return (i<rows&&j<cols)?data[i*cols+j]:0;
	}

	double& operator()(int i, int j) {return data[i*cols+j];}
	double& operator[](int k) {return data[k];}

	Matrix col(int j) const{
		Matrix ans(rows,1);
		for(int i=0;i<rows;++i)ans.data[i]=a(i,j);
		return ans;
	}
	Matrix row(int i) const{
		Matrix ans(1,cols);
		for(int j=0;j<cols;++j)ans.data[j]=a(i,j);
		return ans;
	}
	Matrix submat(int a, int b, int c, int d) const {//(a~c,b~d)
		a--,b--,c--,d--;
		Matrix ans(c-a+1,d-b+1);
		for(int i=0;i<ans.rows;++i)
			for(int j=0;j<ans.cols;++j)
				ans.data[i*ans.cols+j]=data[(a+i)*cols+b+j];
		return ans;
	}

	Matrix(const Matrix& rhs){//赋值构造函数
		new_(rhs.rows,rhs.cols);
		for(int i=0;i<rows*cols;++i)data[i]=rhs.data[i];
	}
	Matrix& operator= (const Matrix& rhs) {
		delete_(),new_(rhs.rows,rhs.cols);
		for(int i=0;i<rows*cols;++i)data[i]=rhs.data[i];
		return *this;
	}
	Matrix& operator= (double v) {
		for(int i=0;i<rows*cols;++i)data[i]=v;
		return *this;
	}

	Matrix operator+ (const Matrix& rhs) {
		Matrix ans(Max(rows,rhs.rows),Max(cols,rhs.cols));
		for(int i=0;i<ans.rows;++i)
			for(int j=0;j<ans.cols;++j)
				ans.data[i*ans.cols+j]=a(i,j)+rhs.a(i,j);
		return ans;
	}
	Matrix operator- (const Matrix& rhs) {
		Matrix ans(Max(rows,rhs.rows),Max(cols,rhs.cols));
		for(int i=0;i<ans.rows;++i)
			for(int j=0;j<ans.cols;++j)
				ans.data[i*ans.cols+j]=a(i,j)-rhs.a(i,j);
		return ans;
	}
	Matrix operator* (const Matrix& rhs) {
		Matrix ans(rows,rhs.cols);
		for(int i=0;i<ans.rows;++i)
			for(int j=0;j<ans.cols;++j){
				ans.data[i*ans.cols+j]=0;
				for(int k=0;k<Min(cols,rhs.rows);++k)
					ans.data[i*ans.cols+j]+=data[i*cols+k]*rhs.data[k*rhs.cols+j];
			}
		return ans;
	}

	Matrix& operator+= (const Matrix& rhs) {return *this=*this+rhs;}
	Matrix& operator-= (const Matrix& rhs) {return *this=*this-rhs;}
	Matrix& operator*= (const Matrix& rhs) {return *this=*this*rhs;}

	Matrix operator+ (double v) {
		Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=data[i]+v;
		return ans;
	}
	Matrix operator- (double v) {
		Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=data[i]-v;
		return ans;
	}
	Matrix operator* (double v) {
		Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=data[i]*v;
		return ans;
	}
	Matrix operator/ (double v) {
		Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=data[i]/v;
		return ans;
	}

	Matrix& operator+= (double v) {return *this=*this+v;}
	Matrix& operator-= (double v) {return *this=*this-v;}
	Matrix& operator*= (double v) {return *this=*this*v;}
	Matrix& operator/= (double v) {return *this=*this/v;}

	void print () const {
		printf("  [Matrix size = (%d x %d)]\n", rows, cols);
		for(int i=0;i<rows;puts(""),++i)
			for(int j=0;j<cols;++j)
				//printf("%8.4lf",data[i*cols+j]);
				//cout<<setw(9)<<setprecision(4)<<data[i*cols+j]<<" ";
				cout<<setw(9)<<data[i*cols+j]<<" ";
		puts("");
	}
	
	Matrix Hilbert () {//把自己填充成Hilbert矩阵
		for(int i=0;i<rows;++i)
			for(int j=0;j<cols;++j)
				data[i*cols+j]=1.0/(i+j+1);
		return *this;
	}

	Matrix Trans () const {//返回转置矩阵
		Matrix ans(cols,rows);
		for(int i=0;i<rows;++i)
			for(int j=0;j<cols;++j)
				ans.data[j*rows+i]=data[i*cols+j];
		return ans;
	}

	Matrix Sign () const {//提取符号
		Matrix ans(rows,cols);
		for(int i=0;i<rows*cols;++i)ans.data[i]=(data[i]>0?1:(data[i]<0?-1:0));
		return ans;
	}

	double sum2(){//矩阵所有元素平方和
		double ans=0;
		for(int i=0;i<rows*cols;++i)ans+=data[i]*data[i];
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
		int n=rows;vector<double>v;
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
	double Norm_1(Matrix A){
		double ans=0;
		for(int j=0;j<A.ncol();++j){
			double sum=0;
			for(int i=0;i<A.nrow();++i)sum+=abs(A(i,j));
			ans=max(ans,sum);
		}
		return ans;
	}
	/*求解∞范数*/
	double Norm_inf(){
		return Norm_1(Trans());
	}

	Matrix house(Matrix x,double &beta){//计算householder变换向量
		int n=x.nrow();
		x/=x.Norm_inf();
		Matrix v=x;v[0]=0;//v(2:n)=x(2:n)
		double sigma=v.sum2();//sigma= x(2)^2 + x(3)^2 +...+ x(n)^2
		if(sigma==0)beta=0;
		else{
			double alpha=sqrt(x[0]*x[0]+sigma);//alpha= x(1)^2 + sigma
			if(x[0]<=0)v[0]=x[0]-alpha;//v(1)= x(1)-alpha
			else v[0]=-sigma/(x[0]+alpha);//v(1)= -sigma / (x(1)+alpha)
			beta=2*v[0]*v[0]/(sigma+v[0]*v[0]);//beta= 2 * v(1)^2 / (sigma + v(1)^2)
			v/=v[0];//v= v / v(1)
		}
		return v;
	}
	Matrix QR(){//QR分解(HouseHoder法)
		int m=rows,n=cols;Matrix d=Matrix(n,1);
		for(int j=0;j<n;++j){
			double beta; Matrix v=house(submat(j+1,j+1,m,j+1),beta);
			Matrix tmp=(v*beta) * (v.Trans()*submat(j+1,j+1,m,n));
			for(int k=j;k<m;++k)
				for(int l=j;l<n;++l)
					data[k*cols+l]-=tmp(k-j,l-j);  //A(j:m,j:n)=(I-beta* v*v^T) * A(j:m,j:n)
			d[j]=beta;
			for(int k=j+1;k<m;++k)data[k*cols+j]=v[k-j];
		}
		return d;
	}

private:
	int rows, cols;
	double *data;
	int Max(int a,int b){return a>b?a:b;}
	int Min(int a,int b){return a<b?a:b;}
	double Abs(double x){return x<0?-x:x;}
	void swap(double &x,double &y){double c=x;x=y,y=c;}
};

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
Matrix Jacobi(Matrix A,Matrix b,double eps){
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
Matrix GS_SOR(Matrix A,Matrix b,double eps,double w=0){
	if(w){if(is_print_SOR)printf("(SOR: 终止条件:eps=%.6lf 松弛因子:w=%.4lf)\n",eps,w);}
	else printf("(GS: 终止条件:eps=%.7lf)\n",eps);
	if(A.nrow()!=A.ncol()){puts("Error!!! Jacobi:矩阵大小不合");return Matrix(0,0);}
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
		if(Norm_inf(x-last)<eps){if(is_print_SOR)printf("  Jacobi: 迭代 %d 次\n",O);SOR_cnt=O;return x;}
		last=x;
	}
	if(is_print_SOR)puts("Error!!! Jacobi:迭代超时");SOR_cnt=2e8;return Matrix(0,0);
}

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