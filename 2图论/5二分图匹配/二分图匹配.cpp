无权：
二分图：把一个图的顶点划分为两个不相交集U和V ，使得每一条边都分别连接U、V中的顶点。如果存在这样的划分，则此图为一个二分图。
匹配：一个「匹配」（matching）是一个边的集合，其中任意两条边都没有公共顶点。
最大匹配：一个图所有匹配中，所含匹配边数最多的匹配，称为这个图的最大匹配。
完美匹配：如果一个图的某个匹配中，所有的顶点都是匹配点，那么它就是一个完美匹配。

（点）覆盖集：一个点集，满足原图所有边至少有一个端点在集合里。
极小（点）覆盖集(minimal vertex covering)：本身为点覆盖，其真子集都不是。
最小（点）覆盖集(minimum vertex covering)：点数最少的点覆盖。
点覆盖数(vertex covering number)：最小点覆盖的点数。

（点）独立集：一个点集，满足集合里任两个结点在原图里不相邻。
极大（点）独立集(maximal independent set)：本身为独立集，再加入任何点都不是。
最大（点）独立集(maximum independent set)：点数最多的独立集。
（点）独立数(independent number)：最大独立集的点。

最小点覆盖集跟最大点独立集互补（如果点覆盖集的某条边两个端点都在最小点覆盖集里，补集肯定不在，如果只有一个顶点在集合里，反过来也是只有一个顶点在，故不相邻）

边覆盖集：一个边集，满足原图所有点都至少与集合里的一条边邻接。
极小边覆盖(minimal edge covering)：本身是边覆盖，其真子集都不是。
最小边覆盖(minimum edge covering)：边数最少的边覆盖。
边覆盖数(edge covering number)：最小边覆盖的边数。

边独立集：一个边集，满足边集中的任两边不邻接。
极大边独立集(maximal edge independent set)：本身为边独立集，再加入任何边都不是。
最大边独立集(maximum edge independent set)：边数最多的边独立集。
边独立数(edge independent number)：最大边独立集的边数。
边独立集又称匹配(matching)，相应的有极大匹配(maximal matching)，最大匹配(maximum matching)，匹配数(matching number)。

const int maxn = 1000 + 5; // 单侧顶点的最大数目
// 二分图最大基数匹配
struct BPM {
  int n, m;               // 左右顶点个数
  vector<int> G[maxn];    // 邻接表
  int left[maxn];         // left[i]为右边第i个点的匹配点编号，-1表示不存在
  bool T[maxn];           // T[i]为右边第i个点是否已标记
  int right[maxn];        // 求最小覆盖用
  bool S[maxn];           // 求最小覆盖用
  void init(int n, int m) {
    this->n = n;
    this->m = m;
    for(int i = 0; i < n; i++) G[i].clear();
  }
  void AddEdge(int u, int v) {
    G[u].push_back(v);
  }
  bool match(int u){
    S[u] = true;
    for(int i = 0; i < G[u].size(); i++) {
      int v = G[u][i];
      if (!T[v]){
        T[v] = true;
        if (left[v] == -1 || match(left[v])){
          left[v] = u;
          right[u] = v;
          return true;
        }
      }
    }
    return false;
  }
  // 求最大匹配
  int solve() {
    memset(left, -1, sizeof(left));
    memset(right, -1, sizeof(right));
    int ans = 0;
    for(int u = 0; u < n; u++) { // 从左边结点u开始增广
      memset(S, 0, sizeof(S));
      memset(T, 0, sizeof(T));
      if(match(u)) ans++;
    }
    return ans;
  }
  // 求最小覆盖。X和Y为最小覆盖中的点集
  int mincover(vector<int>& X, vector<int>& Y) {
    int ans = solve();
    memset(S, 0, sizeof(S));
    memset(T, 0, sizeof(T));
    for(int u = 0; u < n; u++)
      if(right[u] == -1) match(u); // 从所有X未盖点出发增广
    for(int u = 0; u < n; u++)
      if(!S[u]) X.push_back(u); // X中的未标记点
    for(int v = 0; v < m; v++)
      if(T[v]) Y.push_back(v);  // Y中的已标记点
   return ans;
  }
};
//poj_1469
/*==================================================*\
| 二分图匹配（Hopcroft－Carp 的算法）
| INIT: g[][]邻接矩阵;
| CALL: res = MaxMatch(); Nx, Ny要初始化！！！
| Mx，My为match
| 时间复杂度为O（V^0.5 E）
\*==================================================*/
/***********************Hopcroft－Carp 算法****************************************/
#include <cstdio>
#include <queue>
#include <cstring>
using namespace std;

const int MAXN = 310;
const int INF = 1 << 28;
bool flag;
int p,n;
int  Mx[MAXN], My[MAXN], Nx, Ny;
int dx[MAXN], dy[MAXN], dis;
bool vst[MAXN],g[110][310];
bool searchP(void)    //BFS
{
    queue <int> Q;
    dis = INF;
    memset(dx, -1, sizeof(dx));
    memset(dy, -1, sizeof(dy));
    for (int i = 1; i <= Nx; i++)
    if (Mx[i] == -1){
       Q.push(i); dx[i] = 0;
    }
    while (!Q.empty()) {
        int u = Q.front(); Q.pop();
        if (dx[u] > dis) break;        //说明该增广路径长度大于dis还没有结束，等待下一次BFS在扩充
           for (int v = 1; v <= Ny; v++)
               if (g[u][v] && dy[v] == -1) {        //v是未匹配点
                  dy[v] = dx[u]+1;
                if (My[v] == -1) dis = dy[v];    //得到本次BFS的最大遍历层次
                else{
                     dx[My[v]] = dy[v]+1;         //v是匹配点，继续延伸
                     Q.push(My[v]);
                     }
                }
    }
    return dis != INF;
}

bool DFS(int u){
    for (int v = 1; v <= Ny; v++)
    if (!vst[v] && g[u][v] && dy[v] == dx[u]+1) {
       vst[v] = 1;
       if (My[v] != -1 && dy[v] == dis) continue;   //层次（也就是增广路径的长度）大于本次查找的dis，是searchP被break的情况，也就是还不确定是否是增广路径，只有等再次调用searchP()在判断。
       if (My[v] == -1 || DFS(My[v])) {     //是增广路径，更新匹配集
       My[v] = u; Mx[u] = v;
       return 1;
       }
    }
 return 0;
}

int MaxMatch(void){
    int res = 0;
    memset(Mx, -1, sizeof(Mx));
    memset(My, -1, sizeof(My));
    while (searchP()) {
          memset(vst, 0, sizeof(vst));
          for (int i = 1; i <= Nx; i++)
              if (Mx[i] == -1 && DFS(i)) res++;   //查找到一个增广路径，匹配数res++
    }
    return res;
}

/**********************************************************************/
int main()
{
    int i,j,k,t,v,cnt;
    scanf("%d",&t);
    while (t--)
    {
          scanf("%d %d", &p, &n);
          for (i = 1; i <= p; i++)
              for (j = 1; j <= n; j++)
                  g[i][j] = false;
          flag = true;
          for (i = 1; i <= p; i++)
          {
              scanf("%d",&k);
              if (k == 0)
                 flag = false;
              while (k--)
              {
                    scanf("%d",&v);
                    g[i][v]  = true;
              }
          }
          Nx = p; Ny = n;
          if (flag)
          {
               cnt = MaxMatch();
               if (cnt == p)
                  printf("YES\n");
               else printf("NO\n");
          }
          else printf("NO\n");
    }

    return 0;
}
带权：
二分图带权最优匹配：对于二分图的每条边都有一个权（非负），要求一种完美匹配方案，
使得所有匹配边的权和最大，记做最优完美匹配。（特殊的，当所有边的权为1时，就是最大完美匹配问题）
const int maxn=505;
const int INF=0x3f3f3f3f;
int n;
int W[maxn][maxn];
int Left[maxn];
int Lx[maxn],Ly[maxn];
int slack[maxn];
bool S[maxn],T[maxn];
bool match(int i)
{
	S[i]=true;
	for(int j=1;j<=n;j++)if(!T[j])
		if(Lx[i]+Ly[j]==W[i][j])
		{
            T[j]=true;
            if(!Left[j]||match(Left[j]))
			{
				Left[j]=i;
				return true;
			}
		}
		else slack[j]=min(slack[j],Lx[i]+Ly[j]-W[i][j]);
	return false;
}
void update()
{
	int a=INF;
	for(int i=1;i<=n;i++)if(!T[i])a=min(a,slack[i]);
	for(int i=1;i<=n;i++)
	{
		if(S[i])Lx[i]-=a;
		if(T[i])Ly[i]+=a;
	}
}
void KM()
{
	for(int i=1;i<=n;i++)
	{
		Left[i]=Lx[i]=Ly[i]=0;
		for(int j=1;j<=n;j++)Lx[i]=max(Lx[i],W[i][j]);
	}
    for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)slack[j]=INF;
		for(;;)
		{
			for(int j=1;j<=n;j++)S[j]=T[j]=0;
			if(match(i))break;else update();
		}
	}
}
