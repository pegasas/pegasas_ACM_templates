KMP：给出两个字符串A（称为模板串）和B（称为子串），长度分别为lenA和lenB，要求在线性时间内，
对于每个A[i]（0<=i<lenA)，求出A[i]往前和B的前缀匹配的最大匹配长度，
记为ex[i]（或者说，ex[i]为满足A[i-z+1..i]==B[0..z-1]的最大的z值）。
KMP的主要目的是求B是不是A的子串，以及若是，B在A中所有出现的位置（当ex[i]=lenB时）。
【算法】
设next[i]为满足B[i-z+1..i]==B[0..z-1]的最大的z值（也就是B的自身匹配）。设目前next[0..lenB-1]与ex[0..i-1]均已求出，要用它们来求ex[i]的值。
根据ex的定义，有A[i-1-ex[i-1]+1..i-1]==B[0..ex[i-1]-1]，这时，若有A[i]==B[ex[i-1]]，则可以直接得到ex[i]=ex[i-1]+1
（因为i-1-ex[i-1]+1即i-ex[i-1]，现在由于A[i]==B[ex[i-1]]，可得A[i-ex[i-1]..i]==B[0..ex[i-1]]，即A[i-ex[i-1]+1-1..i]==B[0..ex[i-1]+1-1]，
所以ex[i]=ex[i-1]+1）。若A[i]!=B[ex[i-1]]？设j=next[ex[i-1]-1]，则根据next定义得B[ex[i-1]-j..ex[i-1]-1]==B[0..j-1]，
又因为A[i-ex[i-1]..i-1]==B[0..ex[i-1]-1]得A[i-j..i-1]==B[ex[i-1]-j..ex[i-1]-1]，这样有A[i-j..i-1]==B[0..j-1]！
也就是此时只需再比较A[i]与B[j]的值是否相等即可，若相等，可得ex[i]=j+1，若仍不相等，则更新j为next[j-1]，
继续比较A[i]与B[j]是否相等……直到A[i]与B[j]相等或直到j==0时，A[i]仍不等于B[j]，此时ex[i]=0。边界：求ex[0]时，初始j（用来代替ex[i-1]）为0。
现在还有一个问题，如何求next？显然next就是以B自身为模板串，B为子串的“自身匹配”，用类似的办法即可，唯一不同的是next[0]=lenB可以直接得到，
求next[1]时，初始j（代替next[i-1]）为0。

void getNext(const char x[],int len)
{
    int i,j;
    j = Next[0] = -1;
    i = 0;
    while(i < len)
    {
        while(-1 != j &&  x[i] != x[j]) j = Next[j];
        ++i,++j;
        if(x[i] == x[j]) Next[i] = Next[j];
        else Next[i] = j;
    }
}
int KMP(const char p[],int plen,const char s[],int slen)
{
    int i,j;
    int ans = 0;
    getNext(p,plen);
    i = j = 0;
    while(i < slen)
    {
        while(-1 != j && s[i] != p[j]) j = Next[j];
        ++i,++j;
        if(j >= plen)
        {
            ans++;
            j = 0;
        }
    }
    return ans;
}
