扩展KMP：给出模板串A和子串B，长度分别为lenA和lenB，要求在线性时间内，对于每个A[i]（0<=i<lenA)，
求出A[i..lenA-1]与B的最长公共前缀长度，记为ex[i]（或者说，ex[i]为满足A[i..i+z-1]==B[0..z-1]的最大的z值）。
扩展KMP可以用来解决很多字符串问题，如求一个字符串的最长回文子串和最长重复子串。
【算法】
设next[i]为满足B[i..i+z-1]==B[0..z-1]的最大的z值（也就是B的自身匹配）。设目前next[0..lenB-1]与ex[0..i-1]均已求出，要用它们来求ex[i]的值。
设p为目前A串中匹配到的最远位置，k为让其匹配到最远位置的值（或者说，k是在0<=i0<i的所有i0值中，
使i0+ex[i0]-1的值最大的一个，p为这个最大值，即k+ex[k]-1），显然，p之后的所有位都是未知的，
也就是目前还无法知道A[p+1..lenA-1]中的任何一位和B的任何一位是否相等。
根据ex的定义可得，A[k..p]==B[0..p-k]，因为i>k，所以又有A[i..p]==B[i-k..p-k]，设L=next[i-k]，
则根据next的定义有B[0..L-1]==B[i-k..i-k+L-1]。考虑i-k+L-1与p-k的关系：
（1）i-k+L-1<p-k，即i+L<=p。这时，由A[i..p]==B[i-k..p-k]可以得到A[i..i+L-1]==B[i-k..i-k+L-1]，
又因为B[0..L-1]==B[i-k..i-k+L-1]所以A[i..i+L-1]==B[0..L-1]，这就说明ex[i]>=L。又由于next的定义可得，A[i+L]必然不等于B[L]
（否则A[i..i+L]==B[0..L]，因为i+L<=p，所以A[i..i+L]==B[i-k..i-k+L]，这样B[0..L]==B[i-k..i-k+L]，故next[i-k]的值应为L+1或更大），
这样，可以直接得到ex[i]=L！
（2）i+k-L+1>=p-k，即i+L>p。这时，首先可以知道A[i..p]和B[0..p-i]是相等的（因为A[i..p]==B[i-k..p-k]，而i+k-L+1>=p-k，
由B[0..L-1]==B[i-k..i-k+L-1]可得B[0..p-i]==B[i-k..p-k]，即A[i..p]==B[0..p-i]），然后，对于A[p+1]和B[p-i+1]是否相等，
目前是不知道的（因为前面已经说过，p是目前A串中匹配到的最远位置，在p之后无法知道任何一位的匹配信息），
因此，要从A[p+1]与B[p-i+1]开始往后继续匹配（设j为目前B的匹配位置的下标，一开始j=p-i+1，每次比较A[i+j]与B[j]是否相等，
直到不相等或者越界为止，此时的j值就是ex[i]的值）。在这种情况下，p的值必然会得到延伸，因此更新k和p的值。
边界：ex[0]的值需要预先求出，然后将初始的k设为0，p设为ex[0]-1。
对于求next数组，也是“自身匹配”，类似KMP的方法处理即可。唯一的不同点也在边界上：可以直接知道next[0]=lenB，next[1]的值预先求出，然后初始k=1，p=ex[1]。
需要严重注意的是，在上述的情况（2）中，本该从A[p+1]与B[p-i+1]开始匹配，但是，若p+1<i，也就是p-i+1<0
（这种情况是有可能发生的，当ex[i-1]=0，且前面的ex值都没有延伸到i及以后的时候）的话，需要将A、B的下标都加1
（因为此时p必然等于i-2，如果A、B的下标用两个变量x、y控制的话，x和y都要加1）！！
void GetNext(char *P,int plen,int *Next)//计算P自身的Next[i]，Next[i]即满足P[0..z]==P[i..i+z-1]的最大z值，即每个位置i与自身从头开始的最长前缀匹配长度
{
    int a=0;
    Next[0]=plen;//初始化,显然0和0的最长前缀匹配长度是Plen
    while(a<plen-1&&P[a]==P[a+1])a++;
    Next[1]=a;//计算Next[1],即位置1与自身从头开始的最长前缀匹配长度a，并且Next[1]=a
    a=1;
    for(int k=2;k<plen;k++)
    {
        int p=a+Next[a]-1,L=Next[k-a];
        /*
        p为之前求出位置i与自身从头开始的最长前缀匹配长度的匹配到的最远位置，a为取到这个最远位置的起点，即a+Next[a]-1==p,
        那么我们知道P[0..Next[a]-1]==P[a..a+Next[a]-1].当前计算的位置k>a,所以P[k-a..Next[a]-1]==P[k..a+Next[a]-1],对于k-a这个
        位置，又有P[0..Next[k-a]-1]==P[k-a..k-a+Next[k-a]-1].如果k-a+Next[k-a]-1<Next[a]-1,如果令L=Next[k-a],就变成k+L-1<P的话，那么
        P[k-a..k-a+Next[k-a]-1]==P[k..k+Next[k-a]-1]==P[0..Next[k-a]-1]，即表示从k开始Next[k-a]长度的前缀与P从头开始Next[k-a]前缀
        相等，并且不可能更长，否则违反Next[k-a]的"最长"前缀匹配长度的定义。如果k-a+Next[k-a]-1<Next[a]-1,如果令L=Next[k-a],就变成k+L-1>=P的话,
        当前已经匹配到的长度是p-k+1，之后我们就暴力匹配并且此时取到匹配最远位置变成了k,并且利用计算的结果更新Next[k]值。
        P[0..L-1]==P[k-a..k-a+L-1]当前要计算位置k为起点与P开头匹配的最长前缀长度。
        (1)
        */
        if((k-1)+L>=p)
        {
            int j=(p-k+1)>0? p-k+1:0;
            while(k+j<plen&&P[k+j]==P[j])j++;
            Next[k]=j;
            a=k;
        }
        else Next[k]=L;
    }
}
//预处理T的Next数组
void GetExtend(char *T,char *P,int tlen,int plen,int *Next,int *extend)//计算T和P的extend[i],extend[i]表示T[i..i+len-1]==P[0..len-1]的最大len值，即T的位置i开始的和P从头开始的最长匹配前缀长度
{
    int a=0;
    GetNext(P,plen,Next);
    int MinLen=tlen<plen?tlen:plen;
    while(a<MinLen&&T[a]==P[a])a++;
    extend[0]=a;//初始化extend[0]
    a=0;
    for(int k=1;k<tlen;k++)
    {
        int p=a+extend[a]-1,L=Next[k-a];
        if((k-1)+L>=p)//利用Next和刚才相同的过程求extend[k]
        {
            int j=(p-k+1)>0?p-k+1:0;
            while(k+j<tlen&&j<plen&&T[k+j]==P[j]) j++;
            extend[k]=j;
            a=k;
        }
        else extend[k]=L;
    }
}
//求模式串S与匹配串T的扩展kmp
