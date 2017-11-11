#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>
#include<queue>
#include<algorithm>
#include<map>
#include<iomanip>
#define INF 99999999
using namespace std;

const int MAX=110000+10;
char s[MAX*2];
int p[MAX*2];

int main(){
	while(scanf("%s",s)!=EOF){
		int len=strlen(s),id=0,maxlen=0;
		for(int i=len;i>=0;--i){//插入'#'
			s[i+i+2]=s[i];
			s[i+i+1]='#';
		}//插入了len+1个'#',最终的s长度是1~len+len+1即2*len+1,首尾s[0]和s[2*len+2]要插入不同的字符
		s[0]='*';//s[0]='*',s[len+len+2]='\0',防止在while时p[i]越界
		for(int i=2;i<2*len+1;++i){
			if(p[id]+id>i)p[i]=min(p[2*id-i],p[id]+id-i);
			else p[i]=1;
			while(s[i-p[i]] == s[i+p[i]])++p[i];
			if(id+p[id]<i+p[i])id=i;
			if(maxlen<p[i])maxlen=p[i];
		}
		cout<<maxlen-1<<endl;
	}
	return 0;
}
