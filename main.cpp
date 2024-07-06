#include <cstdio>
#include <cmath>


int bd = 4;
double eps = 0.00001;
double delta_x = 0.01, delta_y = 0.01, CFL = 0.3, dtx, dty;
double c0 = 1.0, u0 = 3, r0 = 1.4, p0 = 1.0, dt = 0.0001, tm = 0, tmm = 0, r = 1.4;
double ma1, ma2, r1, f1p, f1m, g1p, g1m, phi;
int t = 0;
int i, j, k, rk;
double maxx, maxy;
int maxi, maxj;
double rp[310][310];
double u[310][110], v[310][110], rho[310][310], pres[310][110], cc[310][110];
double u_x[310][110], v_x[310][110], rho_x[310][110], pres_x[310][110], cc_x[310][110];
double u_y[310][110], v_y[310][110], rho_y[310][110], pres_y[310][110], cc_y[310][110];
int map[310][110];
double lmd[5], lmdp[5], lmdm[5], miu[5], miup[5], mium[5];
double qx0p[5], qx1p[5], qx2p[5], qx0m[5], qx1m[5], qx2m[5], ISx0p[5], ISx1p[5], ISx2p[5], ISx0m[5], ISx1m[5], ISx2m[5];
double qy0p[5], qy1p[5], qy2p[5], qy0m[5], qy1m[5], qy2m[5], ISy0p[5], ISy1p[5], ISy2p[5], ISy0m[5], ISy1m[5], ISy2m[5];
double fplus[310][110][5], fminus[310][110][5], gplus[310][110][5], gminus[310][110][5], un[310][110][5], utemp[310][110][5], fx[310][110][5], fy[310][110][5];
double ft[5], gt[5], t1[5], t2[5];
double om0 = 0.1, om1 = 0.6, om2 = 0.3, al0p, al1p, al2p, al0m, al1m, al2m, totalp, totalm, fom0p, fom1p, fom2p, fom0m, fom1m, fom2m;

void init(double u0, double r0, double p0);
void begin();
void loop(int t);
void Runge_Kutta(int i, int j, int rk);


int main()
{
    init(u0, r0, p0);
    begin();
}
void init(double u0, double r0, double p0)   //初始值设定 
{
    int i, j;
    for (i = 0 + bd; i <= 300 + bd; i++)
        for (j = 0 + bd; j <= 100 + bd; j++)
            map[i][j] = 1;
    for (i = 61 + bd; i <= 300 + bd; i++)
        for (j = 0 + bd; j <= 19 + bd; j++)
            map[i][j] = 0;
    for (i = -3 + bd; i <= 303 + bd; i++)
        for (j = -3 + bd; j <= 103 + bd; j++)
        {
            u[i][j] = u0; v[i][j] = 0; rho[i][j] = r0; pres[i][j] = p0; cc[i][j] = c0;
        }
}
double max(double a, double b, double c)
{
    if ((a >= b) && (a >= c)) return(a);
    if (b >= c)
    {
        maxi = i; maxj = j; return(b);
    }
    maxi = i; maxj = j;
    return(c);
}
void side(int t)    //处理边界条件 
{
    int i, j;
    for (i = 0 + bd; i <= 100 + bd; i++)
    {              //左边界
        u[-1 + bd][i] = u0; v[-1 + bd][i] = 0; pres[-1 + bd][i] = p0; rho[-1 + bd][i] = r0; cc[-1 + bd][i] = c0;
        u[-2 + bd][i] = u0; v[-2 + bd][i] = 0; pres[-2 + bd][i] = p0; rho[-2 + bd][i] = r0; cc[-2 + bd][i] = c0;
        u[-3 + bd][i] = u0; v[-3 + bd][i] = 0; pres[-3 + bd][i] = p0; rho[-3 + bd][i] = r0; cc[-3 + bd][i] = c0;
        u[-4 + bd][i] = u0; v[-4 + bd][i] = 0; pres[-4 + bd][i] = p0; rho[-4 + bd][i] = r0; cc[-4 + bd][i] = c0;
    }
    for (i = 0 + bd; i <= 100 + bd; i++)
    {             //右边界 
        u[301 + bd][i] = u[300 + bd][i]; v[301 + bd][i] = v[300 + bd][i]; pres[301 + bd][i] = pres[300 + bd][i]; rho[301 + bd][i] = rho[300 + bd][i]; cc[301 + bd][i] = cc[300 + bd][i];
        u[302 + bd][i] = u[300 + bd][i]; v[302 + bd][i] = v[300 + bd][i]; pres[302 + bd][i] = pres[300 + bd][i]; rho[302 + bd][i] = rho[300 + bd][i]; cc[302 + bd][i] = cc[300 + bd][i];
        u[303 + bd][i] = u[300 + bd][i]; v[303 + bd][i] = v[300 + bd][i]; pres[303 + bd][i] = pres[300 + bd][i]; rho[303 + bd][i] = rho[300 + bd][i]; cc[303 + bd][i] = cc[300 + bd][i];
        u[304 + bd][i] = u[300 + bd][i]; v[304 + bd][i] = v[300 + bd][i]; pres[304 + bd][i] = pres[300 + bd][i]; rho[304 + bd][i] = rho[300 + bd][i]; cc[304 + bd][i] = cc[300 + bd][i];
    }
    for (i = 0 + bd; i <= 300 + bd; i++)
    {              //上边界
        u[i][101 + bd] = u[i][99 + bd]; v[i][101 + bd] = -v[i][99 + bd]; pres[i][101 + bd] = pres[i][99 + bd]; rho[i][101 + bd] = rho[i][99 + bd]; cc[i][101 + bd] = cc[i][99 + bd];
        u[i][102 + bd] = u[i][98 + bd]; v[i][102 + bd] = -v[i][98 + bd]; pres[i][102 + bd] = pres[i][98 + bd]; rho[i][102 + bd] = rho[i][98 + bd]; cc[i][102 + bd] = cc[i][98 + bd];
        u[i][103 + bd] = u[i][97 + bd]; v[i][103 + bd] = -v[i][97 + bd]; pres[i][103 + bd] = pres[i][97 + bd]; rho[i][103 + bd] = rho[i][97 + bd]; cc[i][103 + bd] = cc[i][97 + bd];
        u[i][104 + bd] = u[i][96 + bd]; v[i][104 + bd] = -v[i][96 + bd]; pres[i][104 + bd] = pres[i][96 + bd]; rho[i][104 + bd] = rho[i][96 + bd]; cc[i][104 + bd] = cc[i][96 + bd];
    }
    for (i = 0 + bd; i <= 60 + bd; i++)
    {             //下边界1 
        u[i][-1 + bd] = u[i][1 + bd]; v[i][bd -1] = -v[i][1 + bd]; pres[i][-1 + bd] = pres[i][1 + bd]; rho[i][-1 + bd] = rho[i][1 + bd]; cc[i][-1 + bd] = cc[i][1 + bd];
        u[i][-2 + bd] = u[i][2 + bd]; v[i][bd -2] = -v[i][2 + bd]; pres[i][-2 + bd] = pres[i][2 + bd]; rho[i][-2 + bd] = rho[i][2 + bd]; cc[i][-2 + bd] = cc[i][2 + bd];
        u[i][-3 + bd] = u[i][3 + bd]; v[i][bd -3] = -v[i][3 + bd]; pres[i][-3 + bd] = pres[i][3 + bd]; rho[i][-3 + bd] = rho[i][3 + bd]; cc[i][-3 + bd] = cc[i][3 + bd];
        u[i][-4 + bd] = u[i][4 + bd]; v[i][bd -4] = -v[i][4 + bd]; pres[i][-4 + bd] = pres[i][4 + bd]; rho[i][-4 + bd] = rho[i][4 + bd]; cc[i][-4 + bd] = cc[i][4 + bd];
    }
    for (i = 61 + bd; i <= 300 + bd; i++)
    {            //下边界2 
        u[i][19 + bd] = u[i][21 + bd]; v[i][19 + bd] = -v[i][21 + bd]; pres[i][19 + bd] = pres[i][21 + bd]; rho[i][19 + bd] = rho[i][21 + bd]; cc[i][19 + bd] = cc[i][21 + bd];
        u[i][18 + bd] = u[i][22 + bd]; v[i][18 + bd] = -v[i][22 + bd]; pres[i][18 + bd] = pres[i][22 + bd]; rho[i][18 + bd] = rho[i][22 + bd]; cc[i][18 + bd] = cc[i][22 + bd];
        u[i][17 + bd] = u[i][23 + bd]; v[i][17 + bd] = -v[i][23 + bd]; pres[i][17 + bd] = pres[i][23 + bd]; rho[i][17 + bd] = rho[i][23 + bd]; cc[i][17 + bd] = cc[i][23 + bd];
        u[i][16 + bd] = u[i][24 + bd]; v[i][16 + bd] = -v[i][24 + bd]; pres[i][16 + bd] = pres[i][24 + bd]; rho[i][16 + bd] = rho[i][24 + bd]; cc[i][16 + bd] = cc[i][24 + bd];
    }
    for (i = 0 + bd; i <= 19 + bd; i++)
    {            //台阶边界 
        u[61 + bd][i] = -u[59 + bd][i]; v[61 + bd][i] = v[59 + bd][i]; pres[61 + bd][i] = pres[59 + bd][i]; rho[61 + bd][i] = rho[59 + bd][i]; cc[61 + bd][i] = cc[59 + bd][i];
        u[62 + bd][i] = -u[58 + bd][i]; v[62 + bd][i] = v[58 + bd][i]; pres[62 + bd][i] = pres[58 + bd][i]; rho[62 + bd][i] = rho[58 + bd][i]; cc[62 + bd][i] = cc[58 + bd][i];
        u[63 + bd][i] = -u[57 + bd][i]; v[63 + bd][i] = v[57 + bd][i]; pres[63 + bd][i] = pres[57 + bd][i]; rho[63 + bd][i] = rho[57 + bd][i]; cc[63 + bd][i] = cc[57 + bd][i];
        u[64 + bd][i] = -u[56 + bd][i]; v[64 + bd][i] = v[56 + bd][i]; pres[64 + bd][i] = pres[56 + bd][i]; rho[64 + bd][i] = rho[56 + bd][i]; cc[64 + bd][i] = cc[56 + bd][i];
    }
}
void Corner_x()    //x方向角点处理 
{
    int i, j;
    for (i = 16 + bd; i <= 19 + bd; i++)
    {            //台阶边界 
        u[61 + bd][i] = -u[59 + bd][i]; v[61 + bd][i] = v[59 + bd][i]; pres[61 + bd][i] = pres[59 + bd][i]; rho[61 + bd][i] = rho[59 + bd][i]; cc[61 + bd][i] = cc[59 + bd][i];
        u[62 + bd][i] = -u[58 + bd][i]; v[62 + bd][i] = v[58 + bd][i]; pres[62 + bd][i] = pres[58 + bd][i]; rho[62 + bd][i] = rho[58 + bd][i]; cc[62 + bd][i] = cc[58 + bd][i];
        u[63 + bd][i] = -u[57 + bd][i]; v[63 + bd][i] = v[57 + bd][i]; pres[63 + bd][i] = pres[57 + bd][i]; rho[63 + bd][i] = rho[57 + bd][i]; cc[63 + bd][i] = cc[57 + bd][i];
        u[64 + bd][i] = -u[56 + bd][i]; v[64 + bd][i] = v[56 + bd][i]; pres[64 + bd][i] = pres[56 + bd][i]; rho[64 + bd][i] = rho[56 + bd][i]; cc[64 + bd][i] = cc[56 + bd][i];
    }
}
void Corner_y()    //y方向角点处理 
{
    int i, j;
    for (i = 61 + bd; i <= 64 + bd; i++)
    {            //下边界2 
        u[i][19 + bd] = u[i][21 + bd]; v[i][19 + bd] = -v[i][21 + bd]; pres[i][19 + bd] = pres[i][21 + bd]; rho[i][19 + bd] = rho[i][21 + bd]; cc[i][19 + bd] = cc[i][21 + bd];
        u[i][18 + bd] = u[i][22 + bd]; v[i][18 + bd] = -v[i][22 + bd]; pres[i][18 + bd] = pres[i][22 + bd]; rho[i][18 + bd] = rho[i][22 + bd]; cc[i][18 + bd] = cc[i][22 + bd];
        u[i][17 + bd] = u[i][23 + bd]; v[i][17 + bd] = -v[i][23 + bd]; pres[i][17 + bd] = pres[i][23 + bd]; rho[i][17 + bd] = rho[i][23 + bd]; cc[i][17 + bd] = cc[i][23 + bd];
        u[i][16 + bd] = u[i][24 + bd]; v[i][16 + bd] = -v[i][24 + bd]; pres[i][16 + bd] = pres[i][24 + bd]; rho[i][16 + bd] = rho[i][24 + bd]; cc[i][16 + bd] = cc[i][24 + bd];
    }
}
void Steger_Warming_X(double rho, double u, double v, double c, int i, int j)   //x方向SW分裂 
{
    lmd[1] = u;
    lmd[2] = lmd[1];
    lmd[3] = lmd[1] - c;
    lmd[4] = lmd[1] + c;
    double rr = rho / (r * 2.0);
    maxx = max(maxx, abs(lmd[3]), abs(lmd[4]));
    int k;
    for (k = 1; k <= 4; k++)
    {
        lmdp[k] = 0.5 * (lmd[k] + sqrt(lmd[k] * lmd[k] + eps * eps));
        lmdm[k] = 0.5 * (lmd[k] - sqrt(lmd[k] * lmd[k] + eps * eps));
    }
    fplus[i][j][1] = rr * (2.0 * (r - 1.0) * lmdp[1] + lmdp[3] + lmdp[4]);
    fplus[i][j][2] = rr * (2.0 * (r - 1.0) * u * lmdp[1] + (u - c) * lmdp[3] + (u + c) * lmdp[4]);
    fplus[i][j][3] = rr * (2.0 * (r - 1.0) * v * lmdp[1] + v * lmdp[3] + v * lmdp[4]);
    fplus[i][j][4] = rr * ((r - 1.0) * (u * u + v * v) * lmdp[1] + 0.5 * ((u - c) * (u - c) + v * v) * lmdp[3] + 0.5 * ((u + c) * (u + c) + v * v) * lmdp[4] + ((3.0 - r) / (2.0 * (r - 1)) * c * c * (lmdp[3] + lmdp[4])));
    fminus[i][j][1] = rr * (2.0 * (r - 1.0) * lmdm[1] + lmdm[3] + lmdm[4]);
    fminus[i][j][2] = rr * (2.0 * (r - 1.0) * u * lmdm[1] + (u - c) * lmdm[3] + (u + c) * lmdm[4]);
    fminus[i][j][3] = rr * (2.0 * (r - 1.0) * v * lmdm[1] + v * lmdm[3] + v * lmdm[4]);
    fminus[i][j][4] = rr * ((r - 1.0) * (u * u + v * v) * lmdm[1] + 0.5 * ((u - c) * (u - c) + v * v) * lmdm[3] + 0.5 * ((u + c) * (u + c) + v * v) * lmdm[4] + ((3.0 - r) / (2.0 * (r - 1)) * c * c * (lmdm[3] + lmdm[4])));
}
void Steger_Warming_Y(double rho, double u, double v, double c, int i, int j)   //y方向SW分裂 
{
    miu[1] = v;
    miu[2] = miu[1];
    miu[3] = miu[1] - c;
    miu[4] = miu[1] + c;
    double rr = rho / (r * 2);
    maxy = max(maxy, abs(miu[3]), abs(miu[4]));
    int k;
    for (k = 1; k <= 4; k++)
    {
        miup[k] = 0.5 * (miu[k] + sqrt(miu[k] * miu[k] + eps * eps));
        mium[k] = 0.5 * (miu[k] - sqrt(miu[k] * miu[k] + eps * eps));
    }
    gplus[i][j][1] = rr * (2.0 * (r - 1.0) * miup[1] + miup[3] + miup[4]);
    gplus[i][j][2] = rr * (2.0 * (r - 1.0) * u * miup[1] + u * miup[3] + u * miup[4]);
    gplus[i][j][3] = rr * (2.0 * (r - 1.0) * v * miup[1] + (v - c) * miup[3] + (v + c) * miup[4]);
    gplus[i][j][4] = rr * ((r - 1.0) * (u * u + v * v) * miup[1] + 0.5 * (u * u + (v - c) * (v - c)) * miup[3] + 0.5 * (u * u + (v + c) * (v + c)) * miup[4] + ((3.0 - r) / (2.0 * (r - 1)) * c * c * (miup[3] + miup[4])));
    gminus[i][j][1] = rr * (2.0 * (r - 1.0) * mium[1] + mium[3] + mium[4]);
    gminus[i][j][2] = rr * (2.0 * (r - 1.0) * u * mium[1] + u * mium[3] + u * mium[4]);
    gminus[i][j][3] = rr * (2.0 * (r - 1.0) * v * mium[1] + (v - c) * mium[3] + (v + c) * mium[4]);
    gminus[i][j][4] = rr * ((r - 1.0) * (u * u + v * v) * mium[1] + 0.5 * (u * u + (v - c) * (v - c)) * mium[3] + 0.5 * (u * u + (v + c) * (v + c)) * mium[4] + ((3.0 - r) / (2.0 * (r - 1)) * c * c * (mium[3] + mium[4])));
}

double WENO_X(int i, int j, int k)      //x方向WENO差分 
{
    qx0p[k] = (1.0 / 3.0) * fplus[i - 2][j][k] - (7.0 / 6.0) * fplus[i - 1][j][k] + (11.0 / 6.0) * fplus[i][j][k];        //x方向
    qx1p[k] = -(1.0 / 6.0) * fplus[i - 1][j][k] + (5.0 / 6.0) * fplus[i][j][k] + (1.0 / 3.0) * fplus[i + 1][j][k];
    qx2p[k] = (1.0 / 3.0) * fplus[i][j][k] + (5.0 / 6.0) * fplus[i + 1][j][k] - (1.0 / 6.0) * fplus[i + 2][j][k];

    t1[k] = fplus[i - 2][j][k] - 2.0 * fplus[i - 1][j][k] + fplus[i][j][k]; t2[k] = fplus[i - 2][j][k] - 4.0 * fplus[i - 1][j][k] + 3.0 * fplus[i][j][k];
    ISx0p[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);
    t1[k] = fplus[i - 1][j][k] - 2.0 * fplus[i][j][k] + fplus[i + 1][j][k]; t2[k] = fplus[i - 1][j][k] - fplus[i + 1][j][k];
    ISx1p[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);
    t1[k] = fplus[i][j][k] - 2.0 * fplus[i + 1][j][k] + fplus[i + 2][j][k]; t2[k] = 3.0 * fplus[i][j][k] - 4.0 * fplus[i + 1][j][k] + fplus[i + 2][j][k];
    ISx2p[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);

    qx0m[k] = -(1.0 / 6) * fminus[i - 1][j][k] + (5.0 / 6) * fminus[i][j][k] + (1.0 / 3) * fminus[i + 1][j][k];
    qx1m[k] = (1.0 / 3) * fminus[i][j][k] + (5.0 / 6) * fminus[i + 1][j][k] - (1.0 / 6) * fminus[i + 2][j][k];
    qx2m[k] = (11.0 / 6) * fminus[i + 1][j][k] - (7.0 / 6) * fminus[i + 2][j][k] + (1.0 / 3) * fminus[i + 3][j][k];

    t1[k] = fminus[i - 1][j][k] - 2 * fminus[i][j][k] + fminus[i + 1][j][k]; t2[k] = fminus[i - 1][j][k] - 4 * fminus[i][j][k] + 3 * fminus[i + 1][j][k];
    ISx0m[k] = (13.0 / 12) * (t1[k] * t1[k]) + (1.0 / 4) * (t2[k] * t2[k]);
    t1[k] = fminus[i][j][k] - 2 * fminus[i + 1][j][k] + fminus[i + 2][j][k]; t2[k] = fminus[i + 2][j][k] - fminus[i][j][k];
    ISx1m[k] = (13.0 / 12) * (t1[k] * t1[k]) + (1.0 / 4) * (t2[k] * t2[k]);
    t1[k] = fminus[i + 1][j][k] - 2 * fminus[i + 2][j][k] + fminus[i + 3][j][k]; t2[k] = 3 * fminus[i + 1][j][k] - 4 * fminus[i + 2][j][k] + fminus[i + 3][j][k];
    ISx2m[k] = (13.0 / 12) * (t1[k] * t1[k]) + (1.0 / 4) * (t2[k] * t2[k]);

    al0p = om0 / ((ISx0p[k] + eps) * (ISx0p[k] + eps)); al1p = om1 / ((ISx1p[k] + eps) * (ISx1p[k] + eps)); al2p = om2 / ((ISx2p[k] + eps) * (ISx2p[k] + eps));
    totalp = al0p + al1p + al2p;
    fom0p = al0p / totalp; fom1p = al1p / totalp; fom2p = al2p / totalp;

    al0m = om2 / ((ISx0m[k] + eps) * (ISx0m[k] + eps)); al1m = om1 / ((ISx1m[k] + eps) * (ISx1m[k] + eps)); al2m = om1 / ((ISx2m[k] + eps) * (ISx2m[k] + eps));
    totalm = al0m + al1m + al2m;
    fom0m = al0m / totalm; fom1m = al1m / totalm; fom2m = al2m / totalm;

    return(qx0m[k] * fom0m + qx1m[k] * fom1m + qx2m[k] * fom2m + qx0p[k] * fom0p + qx1p[k] * fom1p + qx2p[k] * fom2p);
}
double WENO_Y(int i, int j, int k)     //y方向WENO差分 
{
    qy0p[k] = (1.0 / 3.0) * gplus[i][j - 2][k] - (7.0 / 6.0) * gplus[i][j - 1][k] + (11.0 / 6.0) * gplus[i][j][k];        //y方向
    qy1p[k] = -(1.0 / 6.0) * gplus[i][j - 1][k] + (5.0 / 6.0) * gplus[i][j][k] + (1.0 / 3.0) * gplus[i][j + 1][k];
    qy2p[k] = (1.0 / 3.0) * gplus[i][j][k] + (5.0 / 6.0) * gplus[i][j + 1][k] - (1.0 / 6.0) * gplus[i][j + 2][k];

    t1[k] = gplus[i][j - 2][k] - 2.0 * gplus[i][j - 1][k] + gplus[i][j][k]; t2[k] = gplus[i][j - 2][k] - 4.0 * gplus[i][j - 1][k] + 3.0 * gplus[i][j][k];
    ISy0p[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);
    t1[k] = gplus[i][j - 1][k] - 2.0 * gplus[i][j][k] + gplus[i][j + 1][k]; t2[k] = gplus[i][j - 1][k] - gplus[i][j + 1][k];
    ISy1p[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);
    t1[k] = gplus[i][j][k] - 2.0 * gplus[i][j + 1][k] + gplus[i][j + 2][k]; t2[k] = 3 * gplus[i][j][k] - 4.0 * gplus[i][j + 1][k] + gplus[i][j + 2][k];
    ISy2p[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);

    qy0m[k] = -(1.0 / 6.0) * gminus[i][j - 1][k] + (5.0 / 6.0) * gminus[i][j][k] + (1.0 / 3.0) * gminus[i][j + 1][k];
    qy1m[k] = (1.0 / 3.0) * gminus[i][j][k] + (5.0 / 6.0) * gminus[i][j + 1][k] - (1.0 / 6.0) * gminus[i][j + 2][k];
    qy2m[k] = (11.0 / 6.0) * gminus[i][j + 1][k] - (7.0 / 6.0) * gminus[i][j + 2][k] + (1.0 / 3.0) * gminus[i][j + 3][k];

    t1[k] = gminus[i][j - 1][k] - 2.0 * gminus[i][j][k] + gminus[i][j + 1][k]; t2[k] = gminus[i][j - 1][k] - 4.0 * gminus[i][j][k] + 3.0 * gminus[i][j + 1][k];
    ISy0m[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);
    t1[k] = gminus[i][j][k] - 2.0 * gminus[i][j + 1][k] + gminus[i][j + 2][k]; t2[k] = gminus[i][j + 2][k] - gminus[i][j][k];
    ISy1m[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);
    t1[k] = gminus[i][j + 1][k] - 2.0 * gminus[i][j + 2][k] + gminus[i][j + 3][k]; t2[k] = 3.0 * gminus[i][j + 1][k] - 4.0 * gminus[i][j + 2][k] + gminus[i][j + 3][k];
    ISy2m[k] = (13.0 / 12.0) * (t1[k] * t1[k]) + (1.0 / 4.0) * (t2[k] * t2[k]);

    al0p = om0 / ((ISy0p[k] + eps) * (ISy0p[k] + eps)); al1p = om1 / ((ISy1p[k] + eps) * (ISy1p[k] + eps)); al2p = om2 / ((ISy2p[k] + eps) * (ISy2p[k] + eps));
    totalp = al0p + al1p + al2p;
    fom0p = al0p / totalp; fom1p = al1p / totalp; fom2p = al2p / totalp;

    al0m = om2 / ((ISy0m[k] + eps) * (ISy0m[k] + eps)); al1m = om1 / ((ISy1m[k] + eps) * (ISy1m[k] + eps)); al2m = om0 / ((ISy2m[k] + eps) * (ISy2m[k] + eps));
    totalm = al0m + al1m + al2m;
    fom0m = al0m / totalm; fom1m = al1m / totalm; fom2m = al2m / totalm;

    return(qy0m[k] * fom0m + qy1m[k] * fom1m + qy2m[k] * fom2m + qy0p[k] * fom0p + qy1p[k] * fom1p + qy2p[k] * fom2p);
}
void begin()
{
    int i, j;
    int b1 = 1, b2 = 1, b3 = 1, b4 = 1, b5 = 1;
    while (t < 40000)
    {
        loop(t);
        t = t + 1;
        tm = tm + dt;
        if (t % 500 == 0)
        {
            printf("Loading... %d%%\n", t / 500);
        }
        if ((t == 4000) && (b1 = 1))
        {
            b1 = 0;
            for (i = 0 + bd; i <= 300; i++)     //输出到output1.txt，用python进行绘图 
                for (j = 0 + bd; j <= 300; j++)
                    rp[i][j] = rho[j][i];
            FILE* p1;
            p1 = fopen("output1.txt", "w");
            for (i = 0 + bd; i <= 100 + bd; i++) {
                for (j = 0 + bd; j <= 300 + bd; j++)
                    fprintf(p1, "%.1f ", rp[i][j]);
                fprintf(p1, "\n");
            }
        }
        if ((t == 8000) && (b2 = 1))
        {
            b2 = 0;
            for (i = 0 + bd; i <= 300; i++)     //输出到output2.txt，用python进行绘图 
                for (j = 0 + bd; j <= 300; j++)
                    rp[i][j] = rho[j][i];
            FILE* p2;
            p2 = fopen("output2.txt", "w");
            for (i = 0 + bd; i <= 100 + bd; i++) {
                for (j = 0 + bd; j <= 300 + bd; j++)
                    fprintf(p2, "%.1f ", rp[i][j]);
                fprintf(p2, "\n");
            }
        }
        if ((t == 12000) && (b3 = 1))
        {
            b3 = 0;
            for (i = 0 + bd; i <= 300; i++)     //输出到output3.txt，用python进行绘图 
                for (j = 0 + bd; j <= 300; j++)
                    rp[i][j] = rho[j][i];
            FILE* p3;
            p3 = fopen("output3.txt", "w");
            for (i = 0 + bd; i <= 100 + bd; i++) {
                for (j = 0 + bd; j <= 300 + bd; j++)
                    fprintf(p3, "%.1f ", rp[i][j]);
                fprintf(p3, "\n");
            }
        }
        if ((t == 24000) && (b4 = 1))
        {
            b4 = 0;
            for (i = 0 + bd; i <= 300; i++)     //输出到output4.txt，用python进行绘图 
                for (j = 0 + bd; j <= 300; j++)
                    rp[i][j] = rho[j][i];
            FILE* p4;
            p4 = fopen("output4.txt", "w");
            for (i = 0 + bd; i <= 100 + bd; i++) {
                for (j = 0 + bd; j <= 300 + bd; j++)
                    fprintf(p4, "%.1f ", rp[i][j]);
                fprintf(p4, "\n");
            }
        }
        if ((t == 40000) && (b5 = 1))
        {
            b5 = 0;
            for (i = 0 + bd; i <= 300; i++)     //输出到output5.txt，用python进行绘图 
                for (j = 0 + bd; j <= 300; j++)
                    rp[i][j] = rho[j][i];
            FILE* p5;
            p5 = fopen("output5.txt", "w");
            for (i = 0 + bd; i <= 100 + bd; i++) {
                for (j = 0 + bd; j <= 300 + bd; j++)
                    fprintf(p5, "%.1f ", rp[i][j]);
                fprintf(p5, "\n");
            }
        }
    }
    printf("Steps: %d\n", t);
    printf("Time: %f", tm);
}
void loop(int t)      //时间循环 
{
    int i, j, rk;
    for (rk = 1; rk <= 3; rk++)
    {
        maxx = 0; maxy = 0;
        side(t);
        //进行通量分裂 
        Corner_x();
        for (j = -4 + bd; j <= 104 + bd; j++)
            for (i = -4 + bd; i <= 304 + bd; i++)
                Steger_Warming_X(rho[i][j], u[i][j], v[i][j], cc[i][j], i, j);
        Corner_y();
        for (j = -4 + bd; j <= 104 + bd; j++)
            for (i = -4 + bd; i <= 304 + bd; i++)
                Steger_Warming_Y(rho[i][j], u[i][j], v[i][j], cc[i][j], i, j);
        //WENO差分格式 
        for (j = -1 + bd; j <= 100 + bd; j++)
            for (i = -1 + bd; i <= 300 + bd; i++)
                for (k = 1; k <= 4; k++)
                {
                    if ((i == 64) && (j == 9))
                        int deb = 0;
                    fx[i][j][k] = WENO_X(i, j, k);
                    fy[i][j][k] = WENO_Y(i, j, k);
                }
        //进行时间推进 
        for (j = 0 + bd; j <= 100 + bd; j++)
            for (i = 0 + bd; i <= 300 + bd; i++)
                if (map[i][j] == 1)
                {
                    Runge_Kutta(i, j, rk);
                }
    }
}
void Runge_Kutta(int i, int j, int rk)    //三步R-K法 
{
    double rt = 0, ut = 0, vt = 0, pt = 0;
    int k;
    if (rk == 1)
    {
        if (t == 0)
        {
            rt = rho[i][j]; ut = u[i][j]; vt = v[i][j]; pt = pres[i][j];
            un[i][j][1] = rt; un[i][j][2] = rt * ut; un[i][j][3] = rt * vt; un[i][j][4] = pt / (r - 1) + (ut * ut + vt * vt) * 0.5 * rt;
        }
        else
            for (k = 1; k <= 4; k++) un[i][j][k] = utemp[i][j][k];
        for (k = 1; k <= 4; k++)
            utemp[i][j][k] = un[i][j][k] - dt * ((fx[i][j][k] - fx[i - 1][j][k]) / delta_x + (fy[i][j][k] - fy[i][j - 1][k]) / delta_y);
        rho[i][j] = utemp[i][j][1];
        u[i][j] = utemp[i][j][2] / rho[i][j]; v[i][j] = utemp[i][j][3] / rho[i][j];
        pres[i][j] = (r - 1) * (utemp[i][j][4] - 0.5 * utemp[i][j][1] * ((u[i][j]) * (u[i][j]) + (v[i][j] * v[i][j])));
        if (pres[i][j] < eps) pres[i][j] = eps;
        cc[i][j] = sqrt(abs(r * pres[i][j] / rho[i][j]));
    }
    else if (rk == 2)
    {
        for (k = 1; k <= 4; k++)
            utemp[i][j][k] = (3.0 / 4) * un[i][j][k] + (1.0 / 4) * utemp[i][j][k] - (1.0 / 4) * dt * ((fx[i][j][k] - fx[i - 1][j][k]) / delta_x + (fy[i][j][k] - fy[i][j - 1][k]) / delta_y);
        rho[i][j] = utemp[i][j][1];
        u[i][j] = utemp[i][j][2] / rho[i][j]; v[i][j] = utemp[i][j][3] / rho[i][j];
        pres[i][j] = (r - 1) * (utemp[i][j][4] - 0.5 * utemp[i][j][1] * ((u[i][j]) * (u[i][j]) + (v[i][j] * v[i][j])));
        if (pres[i][j] < eps) pres[i][j] = eps;
        cc[i][j] = sqrt(abs(r * pres[i][j] / rho[i][j]));
    }
    else
    {
        for (k = 1; k <= 4; k++)
            utemp[i][j][k] = (1.0 / 3) * un[i][j][k] + (2.0 / 3) * utemp[i][j][k] - (2.0 / 3) * dt * ((fx[i][j][k] - fx[i - 1][j][k]) / delta_x + (fy[i][j][k] - fy[i][j - 1][k]) / delta_y);
        rho[i][j] = utemp[i][j][1];
        u[i][j] = utemp[i][j][2] / rho[i][j]; v[i][j] = utemp[i][j][3] / rho[i][j];
        pres[i][j] = (r - 1) * (utemp[i][j][4] - 0.5 * utemp[i][j][1] * ((u[i][j]) * (u[i][j]) + (v[i][j] * v[i][j])));
        if (pres[i][j] < eps) pres[i][j] = eps;
        cc[i][j] = sqrt(abs(r * pres[i][j] / rho[i][j]));
    }
}





