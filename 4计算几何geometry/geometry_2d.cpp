#include<bits/stdc++.h>
using namespace std;
//判断x大于，小于或者等于0
const double eps = 1e-10;
int dcmp(double x)
{
    if(fabs(x)<eps)return 0;
    else return x<0?-1:1;
}
//PI的值
const double PI=acos(-1.0);
//角度deg换成弧度
double torad(double deg){return deg/180*PI;}


//点定义
struct Point
{
    double x,y;
    Point(double x=0, double y=0):x(x),y(y){}
};
//向量定义
typedef Point Vector;
Vector operator + (Vector A, Vector B) { return Vector(A.x+B.x, A.y+B.y); }
Vector operator - (Point A, Point B) { return Vector(A.x-B.x, A.y-B.y); }
Vector operator * (Vector A, double p) { return Vector(A.x*p, A.y*p); }
Vector operator / (Vector A, double p) { return Vector(A.x/p, A.y/p); }
double Dot(const Vector& A, const Vector& B) {return A.x*B.x+A.y*B.y;}
double Length(const Vector& A){return sqrt(Dot(A, A));}
double Angle(const Vector& A,const Vector& B){return acos(Dot(A, B)/Length(A)/Length(B));}
double Cross(const Vector& A,const Vector& B){return A.x*B.y-A.y*B.x;}
Vector Normal(Vector A)
{
    double L = Length(A);
    return Vector(-A.y/L, A.x/L);
}
Point GetLineIntersection(const Point& P,const Point& v,const Point& Q,const Point& w)
{
    Vector u=P-Q;
    double t=Cross(w,u)/Cross(v,w);
    return P+v*t;
}
Vector Rotate(const Vector& A,double rad)
{
    return Vector(A.x*cos(rad)-A.y*sin(rad),A.x*sin(rad)+A.y*cos(rad));
}
// 理论上这个“小于”运算符是错的，因为可能有三个点a, b, c, a和b很接近（即a<b好b<a都不成立），b和c很接近，但a和c不接近
// 所以使用这种“小于”运算符的前提是能排除上述情况
bool operator < (const Point& a,const Point& b)
{
    return a.x<b.x||(a.x==b.x&&a.y<b.y);
}
bool operator == (const Point& a, const Point &b)
{
    return dcmp(a.x-b.x)==0&&dcmp(a.y-b.y)==0;
}
bool SegmentProperIntersection(const Point& a1, const Point& a2, const Point& b1, const Point& b2)
{
    double c1=Cross(a2-a1,b1-a1),c2=Cross(a2-a1,b2-a1),
    c3=Cross(b2-b1,a1-b1),c4=Cross(b2-b1,a2-b1);
    return dcmp(c1)*dcmp(c2)<0&&dcmp(c3)*dcmp(c4)<0;
}
bool OnSegment(const Point& p, const Point& a1, const Point& a2)
{
    return dcmp(Cross(a1-p, a2-p))==0&&dcmp(Dot(a1-p, a2-p))<0;
}
double DistanceToSegment(const Point& P, const Point& A, const Point& B)
{
    if(A==B)return Length(P-A);
    Vector v1=B-A,v2=P-A,v3=P-B;
    if(dcmp(Dot(v1,v2))<0)return Length(v2);
    else if(dcmp(Dot(v1,v3))>0)return Length(v3);
    else return fabs(Cross(v1,v2))/Length(v1);
}
//得到点在向量AB上的投影
Point GetLineProjection(Point P, Point A, Point B)
{
    Vector v=B-A;
    return A+v*(Dot(v,P-A)/Dot(v, v));
}
//点P到直线AB的距离
double DistanceToLine(Point P,Point A,Point B)
{
    Vector v1=B-A,v2=P-A;
    return fabs(Cross(v1,v2))/Length(v1);//如果不取绝对值，得到的是有向距离
}
//求向量v的角度
double angle(Vector v)
{
    return atan2(v.y,v.x);
}
//线段AB的长度
double Dist2(const Point& A, const Point& B)
{
    return (A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y);
}
//直线定义
struct Line
{
    Point p;//直线上任意一点
    Vector v;//方向向量
    double ang;//极角，即从x正半轴旋转到向量v所需要的角（弧度）
    Line(){}
    Line(Point p,Vector v):p(p),v(v){ang=atan2(v.y,v.x);}
    Point point(double t){return p+v*t;}
    Line move(double d){return Line(p+Normal(v)*d,v);}
    bool operator < (const Line& L) const
    {
        return ang < L.ang;
    }
};
Line getLine(double x1,double y1,double x2,double y2)
{
    Point p1(x1,y1);
    Point p2(x2,y2);
    return Line(p1,p2-p1);
}
//二直线交点，假定交点惟一存在
Point GetLineIntersection(const Line& a,const Line& b)
{
    Vector u = a.p-b.p;
    double t = Cross(b.v, u) / Cross(a.v, b.v);
    return a.p+a.v*t;
}
//圆定义
struct Circle
{
    Point c;
    double r;
    Circle(Point c,double r):c(c),r(r){}
    Point point(double a){return Point(c.x+cos(a)*r,c.y+sin(a)*r);}
};
//读一个点
Point read_point()
{
    double x,y;
    scanf("%lf%lf",&x,&y);
    return Point(x,y);
};

//求直线L和圆C的交点
int getLineCircleIntersection(Line L,Circle C,double& t1,double& t2,vector<Point>& sol)
{
    double a=L.v.x,b=L.p.x-C.c.x,c=L.v.y,d=L.p.y-C.c.y;
    double e=a*a+c*c,f=2*(a*b+c*d),g=b*b+d*d-C.r*C.r;
    double delta=f*f-4*e*g;//判别式
    if(dcmp(delta)<0)return 0;//相离
    if(dcmp(delta)==0)//相切
    {
        t1=t2=-f/(2*e);
        sol.push_back(L.point(t1));
        return 1;
    }
    //相交
    t1 = (-f - sqrt(delta)) / (2 * e); sol.push_back(L.point(t1));
    t2 = (-f + sqrt(delta)) / (2 * e); sol.push_back(L.point(t2));
    return 2;
}
//求圆c1,c2的全部交点
int getCircleCircleIntersection(Circle C1, Circle C2, vector<Point>& sol)
{
    double d=Length(C1.c-C2.c);
    if(dcmp(d)==0)
    {
        if(dcmp(C1.r - C2.r) == 0)return -1;//重合，无穷多交点
        return 0;//同心圆半径不同，交点数为0
    }
    if(dcmp(C1.r+C2.r-d)<0)return 0;//外离
    if(dcmp(fabs(C1.r-C2.r)-d)>0)return 0;//内含
    double a=angle(C2.c - C1.c);
    double da=acos((C1.r*C1.r+d*d-C2.r*C2.r)/(2*C1.r*d));
    Point p1=C1.point(a-da),p2=C1.point(a+da);
    sol.push_back(p1);
    if(p1==p2)return 1;
    sol.push_back(p2);
    return 2;
}
const double TWO_PI=PI*2;
double NormalizeAngle(double rad,double center=PI)
{
    return rad-TWO_PI*floor((rad+PI-center)/TWO_PI);
}
//求圆c1,c2的全部交点,交点相对于圆1的极角保存在rad中
void getCircleCircleIntersection(Point c1, double r1, Point c2, double r2, vector<double>& rad)
{
    double d=Length(c1-c2);
    if(dcmp(d)==0)return;//不管是内含还是重合，都不相交
    if(dcmp(r1+r2-d)<0)return;
    if(dcmp(fabs(r1-r2)-d)>0)return;
    double a=angle(c2 - c1);
    double da=acos((r1*r1 + d*d - r2*r2) / (2*r1*d));
    rad.push_back(NormalizeAngle(a-da));
    rad.push_back(NormalizeAngle(a+da));
}
//求三角形p1,p2,p3的外接圆
Circle CircumscribedCircle(Point p1, Point p2, Point p3)
{
    double Bx=p2.x-p1.x,By=p2.y-p1.y;
    double Cx=p3.x-p1.x,Cy=p3.y-p1.y;
    double D=2*(Bx*Cy-By*Cx);
    double cx=(Cy*(Bx*Bx+By*By)-By*(Cx*Cx+Cy*Cy))/D+p1.x;
    double cy=(Bx*(Cx*Cx+Cy*Cy)-Cx*(Bx*Bx+By*By))/D+p1.y;
    Point p=Point(cx,cy);
    return Circle(p,Length(p1-p));
}
//求三角形p1,p2,p3的内切圆
Circle InscribedCircle(Point p1, Point p2, Point p3)
{
    double a=Length(p2-p3);
    double b=Length(p3-p1);
    double c=Length(p1-p2);
    Point p=(p1*a+p2*b+p3*c)/(a+b+c);
    return Circle(p,DistanceToLine(p,p1,p2));
}
//过点p到圆C的切线。v[i]是第i条切线的向量。返回切线条数
int getTangents(Point p, Circle C, Vector* v)
{
    Vector u=C.c-p;
    double dist=Length(u);
    if(dist<C.r)return 0;
    else if(dcmp(dist - C.r)==0) //p在圆上，只有一条切线
    {
        v[0] = Rotate(u, PI/2);
        return 1;
    }
    else
    {
        double ang=asin(C.r/dist);
        v[0]=Rotate(u,-ang);
        v[1]=Rotate(u,+ang);
        return 2;
    }
}
//给出所有经过点p并与直线L相切的半径为r的圆
vector<Point> CircleThroughPointTangentToLineGivenRadius(Point p, Line L, double r)
{
    vector<Point> ans;
    double t1, t2;
    getLineCircleIntersection(L.move(-r), Circle(p, r), t1, t2, ans);
    getLineCircleIntersection(L.move(r), Circle(p, r), t1, t2, ans);
    return ans;
}
//给出两条不平行直线(x1,y1)-(x2,y2)与(x3,y3)-(x4,y4),求所有半径为r并且同时与这两条直线相切的圆（圆心）
vector<Point> CircleTangentToLinesGivenRadius(Line a, Line b, double r)
{
    vector<Point> ans;
    Line L1 = a.move(-r), L2 = a.move(r);
    Line L3 = b.move(-r), L4 = b.move(r);
    ans.push_back(GetLineIntersection(L1, L3));
    ans.push_back(GetLineIntersection(L1, L4));
    ans.push_back(GetLineIntersection(L2, L3));
    ans.push_back(GetLineIntersection(L2, L4));
    return ans;
}
//求出所有与圆c1,c2外切的半径为r的圆
vector<Point> CircleTangentToTwoDisjointCirclesWithRadius(Circle c1, Circle c2, double r)
{
    vector<Point> ans;
    Vector v = c2.c-c1.c;
    double dist=Length(v);
    int d=dcmp(dist-c1.r-c2.r-r*2);
    if(d>0)return ans;
    getCircleCircleIntersection(Circle(c1.c, c1.r+r),Circle(c2.c, c2.r+r), ans);
    return ans;
}
// formatting
double lineAngleDegree(Vector v)
{
    double ang=angle(v)*180.0/PI;
    while(dcmp(ang)<0)ang+=360.0;
    while(dcmp(ang-180)>=0)ang-=180.0;
    return ang;
}
typedef vector<Point> Polygon;
// 点集凸包
// 如果不希望在凸包的边上有输入点，把两个 <= 改成 <
// 如果不介意点集被修改，可以改成传递引用
vector<Point> ConvexHull(vector<Point> p)
{
    // 预处理，删除重复点
    sort(p.begin(),p.end());
    p.erase(unique(p.begin(),p.end()),p.end());

    int n=p.size();
    int m=0;
    vector<Point> ch(n+1);
    for(int i=0;i<n;i++)
    {
        while(m>1&&Cross(ch[m-1]-ch[m-2],p[i]-ch[m-2])<=0)m--;
        ch[m++]=p[i];
    }
    int k=m;
    for(int i=n-2;i>=0;i--)
    {
        while(m>k&&Cross(ch[m-1]-ch[m-2],p[i]-ch[m-2])<=0)m--;
        ch[m++]=p[i];
    }
    if(n>1)m--;
    ch.resize(m);
    return ch;
}
//多边形的有向面积
double PolygonArea(vector<Point> p)
{
    double area = 0;
    int n=p.size();
    for(int i=1;i<n-1;i++)
    area+=Cross(p[i]-p[0],p[i+1]-p[0]);
    return area/2;
}
// 过两点p1, p2的直线一般方程ax+by+c=0
// (x2-x1)(y-y1) = (y2-y1)(x-x1)
void getLineGeneralEquation(const Point& p1,const Point& p2,double& a,double& b,double &c)
{
    a=p2.y-p1.y;
    b=p1.x-p2.x;
    c=-a*p1.x-b*p1.y;
}
//判断点p是否在多边形poly内
int IsPointInPolygon(const Point& p, const vector<Point>& poly)
{
    int wn = 0;
    int n = poly.size();
    for(int i = 0; i < n; i++)
    {
        const Point& p1 = poly[i];
        const Point& p2 = poly[(i+1)%n];
        if(p1 == p || p2 == p || OnSegment(p, p1, p2)) return -1; // 在边界上
        int k = dcmp(Cross(p2-p1, p-p1));
        int d1 = dcmp(p1.y - p.y);
        int d2 = dcmp(p2.y - p.y);
        if(k > 0 && d1 <= 0 && d2 > 0) wn++;
        if(k < 0 && d2 <= 0 && d1 > 0) wn--;
    }
    if (wn != 0) return 1; // 内部
    return 0; // 外部
}
//判断两个多边形c1,c2是否不相交
bool ConvexPolygonDisjoint(const vector<Point> ch1, const vector<Point> ch2)
{
    int c1 = ch1.size();
    int c2 = ch2.size();
    for(int i = 0; i < c1; i++)
    if(IsPointInPolygon(ch1[i], ch2) != 0) return false; // 内部或边界上
    for(int i = 0; i < c2; i++)
    if(IsPointInPolygon(ch2[i], ch1) != 0) return false; // 内部或边界上
    for(int i = 0; i < c1; i++)
    for(int j = 0; j < c2; j++)
      if(SegmentProperIntersection(ch1[i], ch1[(i+1)%c1], ch2[j], ch2[(j+1)%c2])) return false;
    return true;
}
//返回点集直径的平方（对踵点对）（旋转卡壳）
double diameter2(vector<Point>& points)
{
    vector<Point> p = ConvexHull(points);
    int n = p.size();
    if(n == 1) return 0;
    if(n == 2) return Dist2(p[0], p[1]);
    p.push_back(p[0]); // 免得取模
    double ans = 0;
    for(int u = 0, v = 1; u < n; u++)
    {
        // 一条直线贴住边p[u]-p[u+1]
        for(;;)
        {
            // 当Area(p[u], p[u+1], p[v+1]) <= Area(p[u], p[u+1], p[v])时停止旋转
            // 即Cross(p[u+1]-p[u], p[v+1]-p[u]) - Cross(p[u+1]-p[u], p[v]-p[u]) <= 0
            // 根据Cross(A,B) - Cross(A,C) = Cross(A,B-C)
            // 化简得Cross(p[u+1]-p[u], p[v+1]-p[v]) <= 0
            double diff = Cross(p[u+1]-p[u], p[v+1]-p[v]);
            if(diff<=0)
            {
                ans=max(ans,Dist2(p[u],p[v]));//u和v是对踵点
                if(diff==0)ans=max(ans,Dist2(p[u],p[v+1]));//diff == 0时u和v+1也是对踵点
                break;
            }
            v = (v + 1) % n;
        }
    }
    return ans;
}
//点p在有向直线L的左边（线上不算）
bool OnLeft(const Line& L,const Point& p)
{
    return Cross(L.v,p-L.p) > 0;
}

// 半平面交主过程
vector<Point> HalfplaneIntersection(vector<Line> L)
{
    int n=L.size();
    sort(L.begin(),L.end());//按极角排序

    int first, last;         // 双端队列的第一个元素和最后一个元素的下标
    vector<Point> p(n);      // p[i]为q[i]和q[i+1]的交点
    vector<Line> q(n);       // 双端队列
    vector<Point> ans;       // 结果

    q[first=last=0] = L[0];  // 双端队列初始化为只有一个半平面L[0]
    for(int i = 1; i < n; i++)
    {
        while(first < last && !OnLeft(L[i], p[last-1])) last--;
        while(first < last && !OnLeft(L[i], p[first])) first++;
        q[++last] = L[i];
        if(fabs(Cross(q[last].v, q[last-1].v)) < eps) // 两向量平行且同向，取内侧的一个
        {
            last--;
            if(OnLeft(q[last], L[i].p)) q[last] = L[i];
        }
        if(first < last) p[last-1] = GetLineIntersection(q[last-1], q[last]);
    }
    while(first < last && !OnLeft(q[first], p[last-1])) last--; // 删除无用平面
    if(last - first <= 1) return ans; // 空集
    p[last] = GetLineIntersection(q[last], q[first]); // 计算首尾两个半平面的交点

    // 从deque复制到输出中
    for(int i = first; i <= last; i++) ans.push_back(p[i]);
    return ans;
}

struct Edge
{
    int from, to; // 起点，终点，左边的面编号
    double ang;
};
const int maxn = 10000 + 10; // 最大边数
// 平面直线图（PSGL）实现
struct PSLG {
  int n, m, face_cnt;
  double x[maxn], y[maxn];
  vector<Edge> edges;
  vector<int> G[maxn];
  int vis[maxn*2];  // 每条边是否已经访问过
  int left[maxn*2]; // 左面的编号
  int prev[maxn*2]; // 相同起点的上一条边（即顺时针旋转碰到的下一条边）的编号

  vector<Polygon> faces;
  double area[maxn]; // 每个polygon的面积

  void init(int n) {
    this->n = n;
    for(int i = 0; i < n; i++) G[i].clear();
    edges.clear();
    faces.clear();
  }

  // 有向线段from->to的极角
  double getAngle(int from, int to) {
    return atan2(y[to]-y[from], x[to]-x[from]);
  }

  void AddEdge(int from, int to) {
    edges.push_back((Edge){from, to, getAngle(from, to)});
    edges.push_back((Edge){to, from, getAngle(to, from)});
    m = edges.size();
    G[from].push_back(m-2);
    G[to].push_back(m-1);
  }

  // 找出faces并计算面积
  void Build() {
    for(int u = 0; u < n; u++) {
      // 给从u出发的各条边按极角排序
      int d = G[u].size();
      for(int i = 0; i < d; i++)
        for(int j = i+1; j < d; j++) // 这里偷个懒，假设从每个点出发的线段不会太多
          if(edges[G[u][i]].ang > edges[G[u][j]].ang) swap(G[u][i], G[u][j]);
      for(int i = 0; i < d; i++)
        prev[G[u][(i+1)%d]] = G[u][i];
    }

    memset(vis, 0, sizeof(vis));
    face_cnt = 0;
    for(int u = 0; u < n; u++)
      for(int i = 0; i < G[u].size(); i++) {
        int e = G[u][i];
        if(!vis[e]) { // 逆时针找圈
          face_cnt++;
          Polygon poly;
          for(;;) {
            vis[e] = 1; left[e] = face_cnt;
            int from = edges[e].from;
            poly.push_back(Point(x[from], y[from]));
            e = prev[e^1];
            if(e == G[u][i]) break;
            assert(vis[e] == 0);
          }
          faces.push_back(poly);
        }
      }

    for(int i = 0; i < faces.size(); i++) {
      area[i] = PolygonArea(faces[i]);
    }
  }
};
// 假定poly没有相邻点重合的情况，只需要删除三点共线的情况
Polygon simplify(const Polygon& poly)
{
    Polygon ans;
    int n = poly.size();
    for(int i = 0; i < n; i++) {
    Point a = poly[i];
    Point b = poly[(i+1)%n];
    Point c = poly[(i+2)%n];
    if(dcmp(Cross(a-b, c-b)) != 0) ans.push_back(b);
    }
    return ans;
}
int main()
{
    return 0;
}
