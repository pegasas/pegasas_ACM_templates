����һ������ͼG������ÿ������xi���2�����2i��2i+1���ֱ��ʾxiΪ�ٺ�xiΪ��
����xiΪ�ٻ���xjΪ���������������2i+1->2j 2j+1->2i
��ʵ�ʲ����У�����һ������(i,j)�����Ǿ͵���solver.add_clause(i^1,j^1)��ʾҪôiΪ�ǣ�ҪôjΪ��
const int maxn = 2000 + 10;

struct TwoSAT {
  int n;
  vector<int> G[maxn*2];
  bool mark[maxn*2];
  int S[maxn*2], c;

  bool dfs(int x) {
    if (mark[x^1]) return false;
    if (mark[x]) return true;
    mark[x] = true;
    S[c++] = x;
    for (int i = 0; i < G[x].size(); i++)
      if (!dfs(G[x][i])) return false;
    return true;
  }

  void init(int n) {
    this->n = n;
    for (int i = 0; i < n*2; i++) G[i].clear();
    memset(mark, 0, sizeof(mark));
  }

  // x = xval or y = yval
  void add_clause(int x, int xval, int y, int yval) {
    x = x * 2 + xval;
    y = y * 2 + yval;
    G[x^1].push_back(y);
    G[y^1].push_back(x);
  }

  bool solve() {
    for(int i = 0; i < n*2; i += 2)
      if(!mark[i] && !mark[i+1]) {
        c = 0;
        if(!dfs(i)) {
          while(c > 0) mark[S[--c]] = false;
          if(!dfs(i+1)) return false;
        }
      }
    return true;
  }
};
