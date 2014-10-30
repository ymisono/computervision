#define NRANSI
#include "nrutil.h"

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
/*
 * Ax=b を解いてベクトル x を求める．
 * A はあらかじめ svdcmp で u[1..m][1..n], w[1..n], v[1..n][1..n] に特異値分解しておく．
 * m, n はそれぞれ A の行数，列数である．
 * 正方行列なら m, n は等しい．
 * b[1..m] は連立方程式の右辺(入力)，x[1..n] は解ベクトル(出力)である．
 * 入力したものはどれも破壊されないので，b を取り替えて何度も呼び出すことができる
 */
{
	int jj,j,i;
	double s,*tmp;

	tmp=dvector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_dvector(tmp,1,n);
}
#undef NRANSI
