void dfs(int u)
{
    vis[u]=1;
    pre_deal_with(u);//��������dfs����ڵ�u�������Ϣ
    for(int i=0;i<G[u].size();i++)//����u���ڽӱ�
    {
        int v=G[u][i];
        if(vis[v]==1)continue;
        dfs(v);
    }
    after_deal_with(u);//��������dfs���غ�ڵ�u�������Ϣ
}
