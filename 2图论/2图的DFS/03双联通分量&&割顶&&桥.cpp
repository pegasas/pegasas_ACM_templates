对于无向图G,如果删除某个点u后，联通分量数目增加，称u为图的割顶(cut vertex)
对于无向图G,如果删除某条边e后，联通分量数目增加，称u为图的桥(bridge)
一个无向图中的每一个极大点（边）双连通子图称作此无向图的点-双连通分量(边-双连通分量)。
双连通分量可细分为：
(1)点-双连通分量:在无向图中两点间至少有两条路径，且路径中（不算头尾）的点不同。不同的点-双连通分量最多有一个公共点，这个点必定是“割顶”。
(2)边-双连通分量:在一个无向图中两点间至少有两条路径，且路径中的边不同。边-双连通分量中一定没有桥。
为什么点-双连通分量必须存边
证明如下：
首先要明确边-双连通分量和点-双连通分量的区别与联系
1.二者都是基于无向图
2.边双连通分量是删边后还连通，而后者是删点
3.点双连通分量一定是边双连通分量（除两点一线的特殊情况），反之不一定
4.点双连通分量可以有公共点，而边双连通分量不能有公共边
由于4，显然，求解边双连通分量只需先一遍dfs求桥，在一遍dfs求点（不经过桥即可）
但如果求点双连通分量，就要更复杂：　　
1.如果存边
根据dfs的性质，每条边都有且只有一次入栈，而由于性质3和性质4，点双连通分量没有公共边，
所以弹出这个点双连通分量里的所有边就一定包含这里面的所有点，而且一定不含其他点双连通分量的边。
因此求解时只需弹出这个点双连通分量里的所有边，并记录这些边的点即可（要判重，一个点可出现多次），正确。　　
2.如果存点
根据dfs的性质，每个点同样有且只有一次入栈。但注意，由于性质4，你将一个点出栈后，
还可能有别的点双连通分量包含它，错误。

点-双连通分量
int pre[maxn],iscut[maxn],bccno[maxn],dfs_clock,bcc_cnt; // 割顶的bccno无意义
vector<int> G[maxn],bcc[maxn];
stack<Edge> S;
int dfs(int u,int fa)
{
    int lowu=pre[u]=++dfs_clock;
    int child=0;
    for(int i=0;i<G[u].size();i++)
    {
        int v=G[u][i];
        Edge e=(Edge){u,v};
        if(!pre[v])// 没有访问过v
        {
            S.push(e);
            child++;
            int lowv=dfs(v,u);
            lowu=min(lowu,lowv); // 用后代的low函数更新自己
            if(lowv>=pre[u])
            {
                iscut[u]=true;
                bcc_cnt++;bcc[bcc_cnt].clear();
                for(;;)
                {
                    Edge x=S.top();S.pop();
                    if(bccno[x.u]!=bcc_cnt){bcc[bcc_cnt].push_back(x.u);bccno[x.u]=bcc_cnt;}
                    if(bccno[x.v]!=bcc_cnt){bcc[bcc_cnt].push_back(x.v);bccno[x.v]=bcc_cnt;}
                    if(x.u==u&&x.v==v)break;
                }
            }
        }
        else if(pre[v]<pre[u]&&v!=fa)
        {
            S.push(e);
            lowu = min(lowu, pre[v]); // 用反向边更新自己
        }
    }
    if(fa<0&&child==1)iscut[u]=0;
    return lowu;
}
边-双联通分量
void dfs(int u,int fa)
{
    int lowu=pre[u]=++dfs_clock;
    S.push(u);
    instack[u]=1;
    for(int i=0;i<G[u].size();i++)
    {
        int v=edges[G[u][i]].v;
        if(G[u][i]==(fa^1))continue;
        if(!pre[v])
        {
            int lowv=dfs(v,i);
            if(lowv>pre[u])isbridge[G[u][i]]=true;
            lowu=min(lowu,lowv);
        }
        else if(instack[v])
            lowu=min(lowu,pre[v]);
    }
    if(pre[u]==lowu)
    {
        bcc_cnt++;
        while(1)
        {
            int v=S.top();S.pop();
            instack[v]=0;
            bccno[v]=bcc_cnt;
            if(v==u)break;
        }
    }
}
