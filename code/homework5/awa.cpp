#include<bits/stdc++.h>
using namespace std; 
#define LL long long
string A[400][400];
int main()
{
	int n=20;
	for(int i=0;i<(n-1)*(n-1);++i){
		A[i][i]="(1+h*h/4)";
		if((i+1)%(n-1))A[i+1][i]=A[i][i+1]="(-1/4.0)";
	}
	for(int i=0;i<(n-2)*(n-1);++i)A[i+n-1][i]=A[i][i+n-1]="(-1/4.0)";
	for(int i=0;i<(n-1);++i){
		for(int j=0;j<(n-1);++j){
			printf("b([%d,%d],1)=: ",i+1,j+1);
			for(int k=0;k<n-1;++k)
				for(int l=0;l<n-1;++l)
					if(A[i*(n-1)+j][k*(n-1)+l].size()){
						cout<<A[i*(n-1)+j][k*(n-1)+l];printf("*u[%d,%d]) + ",k+1,l+1);
					}
			puts("");
		}
		puts("");
	}
}
