#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <algorithm>
using namespace std;
typedef long long ll;
const int inf=0x3f3f3f3f;
const int MAXN=100010;
int sa[MAXN];
int t1[MAXN],t2[MAXN],c[MAXN];
int Rank[MAXN],height[MAXN];
void construct_sa(int s[],int n,int m){
    int i,j,p,*x=t1,*y=t2;
    for(i=0;i<m;i++)c[i]=0;
    for(i=0;i<n;i++)c[x[i]=s[i]]++;
    for(i=1;i<m;i++)c[i]+=c[i-1];
    for(i=n-1;i>=0;i--)sa[--c[x[i]]]=i;
    for(j=1;j<=n;j<<=1){
        p=0;
        for(i=n-j;i<n;i++)y[p++]=i;
        for(i=0;i<n;i++)if(sa[i]>=j)y[p++]=sa[i]-j;
        for(i=0;i<m;i++)c[i]=0;
        for(i=0;i<n;i++)c[x[y[i]]]++;
        for(i=1;i<m;i++)c[i]+=c[i-1];
        for(i=n-1;i>=0;i--)sa[--c[x[y[i]]]]=y[i];
        swap(x,y);
        p=1;x[sa[0]]=0;
        for(i=1;i<n;i++)
            x[sa[i]]=y[sa[i-1]]==y[sa[i]] && y[sa[i-1]+j]==y[sa[i]+j]?p-1:p++;
        if(p>=n)break;
        m=p;
    }
}
void construct_lcp(int s[],int n){
    int k=0;
    for(int i=0;i<=n;i++) Rank[sa[i]]=i;
    for(int i=0;i<n;i++){
        if(k)k--;
        int j=sa[Rank[i]-1];
        while(s[i+k]==s[j+k])k++;
        height[Rank[i]]=k;
    }
}
int fvis[MAXN],pr[MAXN],s[MAXN];
char str[MAXN];
int main(){
    int T,cas=1;
    char ch[10];
    scanf("%d",&T);
    while(T--){
        scanf("%s",ch);
        scanf("%s",str);
        int n=strlen(str),flag=0,fpr;
        for(int i=0;i<=n;i++) s[i]=str[i];
        memset(fvis,0,sizeof(fvis));
        memset(pr,0,sizeof(pr));
        for(int i=n-1;i>=0;i--){
            if(str[i]==ch[0]){
                flag=1;fpr=i;
            }
            if(flag) fvis[i]=1,pr[i]=fpr;
        }
        ll ans=0;
        construct_sa(s,n+1,128);
        construct_lcp(s,n);
        for(int i=0;i<n;i++){
            if(fvis[sa[i]]) ans=ans+n-max((sa[i]+height[i]),pr[sa[i]]);
        }
        printf("Case #%d: %I64d\n",cas++,ans);
    }
    return 0;
}
