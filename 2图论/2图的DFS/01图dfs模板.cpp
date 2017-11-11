void dfs(int u)
{
    vis[u]=1;
    pre_deal_with(u);//处理往下dfs处理节点u的相关信息
    for(int i=0;i<G[u].size();i++)//访问u的邻接表
    {
        int v=G[u][i];
        if(vis[v]==1)continue;
        dfs(v);
    }
    after_deal_with(u);//处理往下dfs返回后节点u的相关信息
}
