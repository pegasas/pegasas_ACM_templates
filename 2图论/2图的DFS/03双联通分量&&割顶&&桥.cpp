��������ͼG,���ɾ��ĳ����u����ͨ������Ŀ���ӣ���uΪͼ�ĸ(cut vertex)
��������ͼG,���ɾ��ĳ����e����ͨ������Ŀ���ӣ���uΪͼ����(bridge)
һ������ͼ�е�ÿһ������㣨�ߣ�˫��ͨ��ͼ����������ͼ�ĵ�-˫��ͨ����(��-˫��ͨ����)��
˫��ͨ������ϸ��Ϊ��
(1)��-˫��ͨ����:������ͼ�����������������·������·���У�����ͷβ���ĵ㲻ͬ����ͬ�ĵ�-˫��ͨ���������һ�������㣬�����ض��ǡ������
(2)��-˫��ͨ����:��һ������ͼ�����������������·������·���еı߲�ͬ����-˫��ͨ������һ��û���š�
Ϊʲô��-˫��ͨ����������
֤�����£�
����Ҫ��ȷ��-˫��ͨ�����͵�-˫��ͨ��������������ϵ
1.���߶��ǻ�������ͼ
2.��˫��ͨ������ɾ�ߺ���ͨ����������ɾ��
3.��˫��ͨ����һ���Ǳ�˫��ͨ������������һ�ߵ��������������֮��һ��
4.��˫��ͨ���������й����㣬����˫��ͨ���������й�����
����4����Ȼ������˫��ͨ����ֻ����һ��dfs���ţ���һ��dfs��㣨�������ż��ɣ�
��������˫��ͨ��������Ҫ�����ӣ�����
1.������
����dfs�����ʣ�ÿ���߶�����ֻ��һ����ջ������������3������4����˫��ͨ����û�й����ߣ�
���Ե��������˫��ͨ����������б߾�һ����������������е㣬����һ������������˫��ͨ�����ıߡ�
������ʱֻ�赯�������˫��ͨ����������бߣ�����¼��Щ�ߵĵ㼴�ɣ�Ҫ���أ�һ����ɳ��ֶ�Σ�����ȷ������
2.������
����dfs�����ʣ�ÿ����ͬ������ֻ��һ����ջ����ע�⣬��������4���㽫һ�����ջ��
�������б�ĵ�˫��ͨ����������������

��-˫��ͨ����
int pre[maxn],iscut[maxn],bccno[maxn],dfs_clock,bcc_cnt; // ���bccno������
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
        if(!pre[v])// û�з��ʹ�v
        {
            S.push(e);
            child++;
            int lowv=dfs(v,u);
            lowu=min(lowu,lowv); // �ú����low���������Լ�
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
            lowu = min(lowu, pre[v]); // �÷���߸����Լ�
        }
    }
    if(fa<0&&child==1)iscut[u]=0;
    return lowu;
}
��-˫��ͨ����
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
