��Ȩ��
����ͼ����һ��ͼ�Ķ��㻮��Ϊ�������ཻ��U��V ��ʹ��ÿһ���߶��ֱ�����U��V�еĶ��㡣������������Ļ��֣����ͼΪһ������ͼ��
ƥ�䣺һ����ƥ�䡹��matching����һ���ߵļ��ϣ��������������߶�û�й������㡣
���ƥ�䣺һ��ͼ����ƥ���У�����ƥ���������ƥ�䣬��Ϊ���ͼ�����ƥ�䡣
����ƥ�䣺���һ��ͼ��ĳ��ƥ���У����еĶ��㶼��ƥ��㣬��ô������һ������ƥ�䡣

���㣩���Ǽ���һ���㼯������ԭͼ���б�������һ���˵��ڼ����
��С���㣩���Ǽ�(minimal vertex covering)������Ϊ�㸲�ǣ������Ӽ������ǡ�
��С���㣩���Ǽ�(minimum vertex covering)���������ٵĵ㸲�ǡ�
�㸲����(vertex covering number)����С�㸲�ǵĵ�����

���㣩��������һ���㼯�����㼯���������������ԭͼ�ﲻ���ڡ�
���󣨵㣩������(maximal independent set)������Ϊ���������ټ����κε㶼���ǡ�
��󣨵㣩������(maximum independent set)���������Ķ�������
���㣩������(independent number)�����������ĵ㡣

��С�㸲�Ǽ����������������������㸲�Ǽ���ĳ���������˵㶼����С�㸲�Ǽ�������϶����ڣ����ֻ��һ�������ڼ����������Ҳ��ֻ��һ�������ڣ��ʲ����ڣ�

�߸��Ǽ���һ���߼�������ԭͼ���е㶼�����뼯�����һ�����ڽӡ�
��С�߸���(minimal edge covering)�������Ǳ߸��ǣ������Ӽ������ǡ�
��С�߸���(minimum edge covering)���������ٵı߸��ǡ�
�߸�����(edge covering number)����С�߸��ǵı�����

�߶�������һ���߼�������߼��е������߲��ڽӡ�
����߶�����(maximal edge independent set)������Ϊ�߶��������ټ����κα߶����ǡ�
���߶�����(maximum edge independent set)���������ı߶�������
�߶�����(edge independent number)�����߶������ı�����
�߶������ֳ�ƥ��(matching)����Ӧ���м���ƥ��(maximal matching)�����ƥ��(maximum matching)��ƥ����(matching number)��

const int maxn = 1000 + 5; // ���ඥ��������Ŀ
// ����ͼ������ƥ��
struct BPM {
  int n, m;               // ���Ҷ������
  vector<int> G[maxn];    // �ڽӱ�
  int left[maxn];         // left[i]Ϊ�ұߵ�i�����ƥ����ţ�-1��ʾ������
  bool T[maxn];           // T[i]Ϊ�ұߵ�i�����Ƿ��ѱ��
  int right[maxn];        // ����С������
  bool S[maxn];           // ����С������
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
  // �����ƥ��
  int solve() {
    memset(left, -1, sizeof(left));
    memset(right, -1, sizeof(right));
    int ans = 0;
    for(int u = 0; u < n; u++) { // ����߽��u��ʼ����
      memset(S, 0, sizeof(S));
      memset(T, 0, sizeof(T));
      if(match(u)) ans++;
    }
    return ans;
  }
  // ����С���ǡ�X��YΪ��С�����еĵ㼯
  int mincover(vector<int>& X, vector<int>& Y) {
    int ans = solve();
    memset(S, 0, sizeof(S));
    memset(T, 0, sizeof(T));
    for(int u = 0; u < n; u++)
      if(right[u] == -1) match(u); // ������Xδ�ǵ��������
    for(int u = 0; u < n; u++)
      if(!S[u]) X.push_back(u); // X�е�δ��ǵ�
    for(int v = 0; v < m; v++)
      if(T[v]) Y.push_back(v);  // Y�е��ѱ�ǵ�
   return ans;
  }
};
//poj_1469
/*==================================================*\
| ����ͼƥ�䣨Hopcroft��Carp ���㷨��
| INIT: g[][]�ڽӾ���;
| CALL: res = MaxMatch(); Nx, NyҪ��ʼ��������
| Mx��MyΪmatch
| ʱ�临�Ӷ�ΪO��V^0.5 E��
\*==================================================*/
/***********************Hopcroft��Carp �㷨****************************************/
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
        if (dx[u] > dis) break;        //˵��������·�����ȴ���dis��û�н������ȴ���һ��BFS������
           for (int v = 1; v <= Ny; v++)
               if (g[u][v] && dy[v] == -1) {        //v��δƥ���
                  dy[v] = dx[u]+1;
                if (My[v] == -1) dis = dy[v];    //�õ�����BFS�����������
                else{
                     dx[My[v]] = dy[v]+1;         //v��ƥ��㣬��������
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
       if (My[v] != -1 && dy[v] == dis) continue;   //��Σ�Ҳ��������·���ĳ��ȣ����ڱ��β��ҵ�dis����searchP��break�������Ҳ���ǻ���ȷ���Ƿ�������·����ֻ�е��ٴε���searchP()���жϡ�
       if (My[v] == -1 || DFS(My[v])) {     //������·��������ƥ�伯
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
              if (Mx[i] == -1 && DFS(i)) res++;   //���ҵ�һ������·����ƥ����res++
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
��Ȩ��
����ͼ��Ȩ����ƥ�䣺���ڶ���ͼ��ÿ���߶���һ��Ȩ���Ǹ�����Ҫ��һ������ƥ�䷽����
ʹ������ƥ��ߵ�Ȩ����󣬼�����������ƥ�䡣������ģ������бߵ�ȨΪ1ʱ�������������ƥ�����⣩
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
