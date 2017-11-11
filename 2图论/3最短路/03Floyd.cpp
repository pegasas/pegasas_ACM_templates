void Floyd()
{
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            dist[i][j]=M[i][j];
    for(int k=0;k<n;k++)
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                if(dist[i][j]<dist[i][k]+dist[k][j])
                    dist[i][j]=min(dist[i][j],dist[i][k]+dist[k][j]);
}
