#include<bits/stdc++.h>
#define ll long long
using namespace std;
const int maxn=100005;
char ch[10];
char str[maxn];
bool vis[maxn];
int pos[maxn];

int s[maxn];
int sa[maxn];
int t1[maxn],t2[maxn],c[maxn];
int Rank[maxn],height[maxn];
void construct_sa(int s[],int n,int m)
{
    int i,j,p,*x=t1,*y=t2;
    for(i=0;i<m;i++)c[i]=0;
    for(i=0;i<n;i++)c[x[i]=s[i]]++;
    for(i=1;i<m;i++)c[i]+=c[i-1];
    for(i=n-1;i>=0;i--)sa[--c[x[i]]]=i;
    for(j=1;j<=n;j<<=1)
    {
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
void construct_lcp(int s[],int n)
{
    int k=0;
    for(int i=0;i<=n;i++)Rank[sa[i]]=i;
    for(int i=0;i<n;i++)
    {
        if(k)k--;
        int j=sa[Rank[i]-1];
        while(s[i+k]==s[j+k])k++;
        height[Rank[i]]=k;
    }
}
int main()
{
//    freopen("input.txt","r",stdin);
    int T;
    cin>>T;
    for(int Case=1;Case<=T;Case++)
    {
        scanf("%s",ch);scanf("%s",str);
        int n=strlen(str);
        for(int i=0;i<=n;i++)s[i]=str[i];
        int flag=0,max_c_pos;
        memset(vis,0,sizeof(vis[0])*n);
        for(int i=n-1;i>=0;i--)
        {
            if(s[i]==ch[0]){flag=1;max_c_pos=i;}
            if(flag){vis[i]=1;pos[i]=max_c_pos;}
        }
        ll ans=0;
        construct_sa(s,n+1,128);
        construct_lcp(s,n);
        for(int i=0;i<=n;i++)
            if(vis[sa[i]])ans+=n-max(sa[i]+height[i],pos[sa[i]]);
        printf("Case #%d: %I64d\n",Case,ans);
    }
    return 0;
}
