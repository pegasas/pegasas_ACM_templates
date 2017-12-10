const int maxn = 100 + 10;
const int INF = 1000000000;

struct Edge {
  int from, to, cap, flow;
};

bool operator < (const Edge& a, const Edge& b) {
  return a.from < b.from || (a.from == b.from && a.to < b.to);
}

struct ISAP {
  int n, m, s, t;
  vector<Edge> edges;
  vector<int> G[maxn];   // �ڽӱ�G[i][j]��ʾ���i�ĵ�j������e�����е����
  bool vis[maxn];        // BFSʹ��
  int d[maxn];           // ����㵽i�ľ���
  int cur[maxn];        // ��ǰ��ָ��
  int p[maxn];          // ������·�ϵ���һ����
  int num[maxn];        // �����ż���

  void AddEdge(int from, int to, int cap) {
    edges.push_back((Edge){from, to, cap, 0});
    edges.push_back((Edge){to, from, 0, 0});
    m = edges.size();
    G[from].push_back(m-2);
    G[to].push_back(m-1);
  }

  bool BFS() {
    memset(vis, 0, sizeof(vis));
    queue<int> Q;
    Q.push(t);
    vis[t] = 1;
    d[t] = 0;
    while(!Q.empty()) {
      int x = Q.front(); Q.pop();
      for(int i = 0; i < G[x].size(); i++) {
        Edge& e = edges[G[x][i]^1];
        if(!vis[e.from] && e.cap > e.flow) {
          vis[e.from] = 1;
          d[e.from] = d[x] + 1;
          Q.push(e.from);
        }
      }
    }
    return vis[s];
  }

  void ClearAll(int n) {
    this->n = n;
    for(int i = 0; i < n; i++) G[i].clear();
    edges.clear();
  }

  void ClearFlow() {
    for(int i = 0; i < edges.size(); i++) edges[i].flow = 0;
  }

  int Augment() {
    int x = t, a = INF;
    while(x != s) {
      Edge& e = edges[p[x]];
      a = min(a, e.cap-e.flow);
      x = edges[p[x]].from;
    }
    x = t;
    while(x != s) {
      edges[p[x]].flow += a;
      edges[p[x]^1].flow -= a;
      x = edges[p[x]].from;
    }
    return a;
  }

  int Maxflow(int s, int t, int need) {
    this->s = s; this->t = t;
    int flow = 0;
    BFS();
    memset(num, 0, sizeof(num));
    for(int i = 0; i < n; i++) num[d[i]]++;
    int x = s;
    memset(cur, 0, sizeof(cur));
    while(d[s] < n) {
      if(x == t) {
        flow += Augment();
        if(flow >= need) return flow;
        x = s;
      }
      int ok = 0;
      for(int i = cur[x]; i < G[x].size(); i++) {
        Edge& e = edges[G[x][i]];
        if(e.cap > e.flow && d[x] == d[e.to] + 1) { // Advance
          ok = 1;
          p[e.to] = G[x][i];
          cur[x] = i; // ע��
          x = e.to;
          break;
        }
      }
      if(!ok) { // Retreat
        int m = n-1; // ��ֵע��
        for(int i = 0; i < G[x].size(); i++) {
          Edge& e = edges[G[x][i]];
          if(e.cap > e.flow) m = min(m, d[e.to]);
        }
        if(--num[d[x]] == 0) break;
        num[d[x] = m+1]++;
        cur[x] = 0; // ע��
        if(x != s) x = edges[p[x]].from;
      }
    }
    return flow;
  }

  vector<int> Mincut() { // call this after maxflow
    BFS();
    vector<int> ans;
    for(int i = 0; i < edges.size(); i++) {
      Edge& e = edges[i];
      if(!vis[e.from] && vis[e.to] && e.cap > 0) ans.push_back(i);
    }
    return ans;
  }

  void Reduce() {
    for(int i = 0; i < edges.size(); i++) edges[i].cap -= edges[i].flow;
  }

  void print() {
    printf("Graph:\n");
    for(int i = 0; i < edges.size(); i++)
      printf("%d->%d, %d, %d\n", edges[i].from, edges[i].to , edges[i].cap, edges[i].flow);
  }
};

