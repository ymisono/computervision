/****************************************************************
 *
 *  Photopmetric Stereo Sample
 *
 *   filename	: photometricStereoSample.c
 *   writer	: Masahito Aoyama
 *   date	: Dec.17, 2007
 *   date	: Nov.19, 2009
 *   date	: Dec.1, 2010
 *   date	: Nov.28, 2011
 *   date	: Oct.9, 2012
 *   date	: Sep.30, 2014
 *   compile    : make PN=photometricStereoSample
 *   example    : ./photometricStereoSample 5 img/sphere{top,045,135,225,315}.pgm
 *              : ./photometricStereoSample 5 img/cone{top,045,135,225,315}.pgm
 *              : ./photometricStereoSample 5 img/torus{top,045,135,225,315}.pgm
 *              : ./photometricStereoSample 17 img/sphere{top,000,030,045,060,090,120,135,150,180,210,225,240,270,300,315,330}.pgm
 *              : ./photometricStereoSample 17 img/cone{top,000,030,045,060,090,120,135,150,180,210,225,240,270,300,315,330}.pgm
 *              : ./photometricStereoSample 17 img/torus{top,000,030,045,060,090,120,135,150,180,210,225,240,270,300,315,330}.pgm
 ****************************************************************/
#include <stdlib.h>
#include <math.h>
#include "header.h"

/**** External Variables ****/
Light light[LIGHTNUM]; // 光源

/**** Prototype Declaration  ****/
/**** Numerical Recipes      ****/
void svdcmp(double **a, int m, int n, double w[], double **v);
void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
/********************************/
/**** Others                 ****/
int main(int argc, char *argv[]);
int usage(char *command);  
int load_light(void);
uchar **readimg(char *iFname, char *magic, int *xsize, int *ysize, int *maxval);
int writeimg(uchar **optr, char *oFname, char *magic, int xsize, int ysize, int maxval); // -> 未使用
int calc_depth(double **V, uchar **iptr[], int inum, int xsize, int ysize);
int set_light(char *inp_fname, double *V);
int make_I_i(double **I, double *i, uchar **iptr[], int inum, int l0, int l1);
int calc_Ii(double *Ii, double **I, double *i, int inum);
int calc_IV(double **IV, double **I, double **V, int inum);
int safety_check(double **g, double **nx, double **ny, double **nz, double **err, double **p, double **q, int xsize, int ysize);
int make_depth(double **dep, double **p, double **q, int xsize, int ysize);
int calc_svd(double **IV, double *w, double **v, double *Ii, double *x, int inum);
int calc_vals(double *g, double *x, double *nx, double *ny, double *nz, double *p, double *q);
int out_all(double **g, double **nx, double **ny, double **nz, double **dep, double **err, int xsize, int ysize);
/****************************************************************/

/**** Routines ****/

int main(int argc, char *argv[])
{
  char    *iFname;  // 入力画像ファイル名
  uchar  **inp_img; // 画像一枚
  uchar ***iptr;    // 画像集合（照明方向が異なる入力画像一式）
  int      inum;    // 入力画像数

  char magic[IHLENGTH];     // 画像ヘッダ用
  int xsize, ysize, maxval; // 画像情報（横，縦，最大輝度値）

  int       i;

  double  **V; // 光源ベクトル

  // 引数のチェック
  if (argc == 1){
    usage(argv[0]);
    exit(1);
  }
  inum = atoi(argv[1]); // 画像数
  if ((inum == 0) || (argc != inum + 2)){
    usage(argv[0]);
    exit(1);
  }

  load_light();
  V  = (double **)dmatrix(1,inum,1,DIM); // 光源用
  iptr   = (uchar ***)malloc(sizeof(uchar **)*inum);

  // 入力画像の読み込み
  for (i = 0; i < atoi(argv[1]); i++){
    iFname = argv[i+2];
    if ((iptr[i] = readimg(iFname, magic, &xsize, &ysize, &maxval)) == NULL) return(-1);
    set_light(iFname, V[i+1]);
  }
	
  // 実際の処理
  if (calc_depth(V, iptr, inum, xsize, ysize) < 0){
    fprintf(stderr, "Error in BODY!!\n");
    exit(1);
  }

  // 領域の解放
  for (i = 0; i < atoi(argv[1]); i++){
    free_cmatrix(iptr[i], 0, xsize-1, 0, ysize-1);
  }
  free(iptr);
  return(0);
}

// 使用法の出力
int usage(char	*command)
{
    fprintf(stderr, "photometricStereo\n");
    fprintf(stderr, "Usage: %s ", command);
    fprintf(stderr, "<number of image> <input_file1> <input_file2> ... <input_file(number of image)>\n");

    return(0);
}

int load_light(void)
{
  // 用意している画像と光源ベクトル
  // 天頂
  light[ 0].angle = "top";
  light[ 0].x =          0.; light[ 0].y = 0.;
  // 0度
  light[ 1].angle = "000";
  light[ 1].x =          1.; light[ 1].y = 0.;
  // 30度
  light[ 2].angle = "030";
  light[ 2].x = sqrt(3.)/2.; light[ 2].y = 1./2.;
  // 45度
  light[ 3].angle = "045";
  light[ 3].x = sqrt(2.)/2.; light[ 3].y = sqrt(2.)/2.;
  // 60度
  light[ 4].angle = "060";
  light[ 4].x =       1./2.; light[ 4].y = sqrt(3.)/2.;
  // 90度
  light[ 5].angle = "090";
  light[ 5].x =          0.; light[ 5].y = 1.;
  // 120度
  light[ 6].angle = "120";
  light[ 6].x =-light[4].x;  light[ 6].y = light[4].y;
  // 135度
  light[ 7].angle = "135";
  light[ 7].x =-light[3].x;  light[ 7].y = light[3].y;
  // 150度
  light[ 8].angle = "150";
  light[ 8].x =-light[2].x;  light[ 8].y = light[2].y;
  // 180度
  light[ 9].angle = "180";
  light[ 9].x =-light[1].x;  light[ 9].y = light[1].y;
  // 210度
  light[10].angle = "210";
  light[10].x = light[8].x;  light[10].y =-light[8].y;
  // 225度
  light[11].angle = "225";
  light[11].x = light[7].x;  light[11].y =-light[7].y;
  // 240度
  light[12].angle = "240";
  light[12].x = light[6].x;  light[12].y =-light[6].y;
  // 270度
  light[13].angle = "270";
  light[13].x = light[5].x;  light[13].y =-light[5].y;
  // 300度
  light[14].angle = "300";
  light[14].x = light[4].x;  light[14].y =-light[4].y;
  // 315度
  light[15].angle = "315";
  light[15].x = light[3].x;  light[15].y =-light[3].y;
  // 330度
  light[16].angle = "330";
  light[16].x = light[2].x;  light[16].y =-light[2].y;

  return(0);
}

// 形状復元を行う関数本体
int calc_depth(double **V, uchar **iptr[], int inum, int xsize, int ysize)
{
  double *i;  // 各画像点における測定結果のベクトル
  double **I;
  double *Ii;
  double **IV;

  double **g;              // アルベド
  double **nx, **ny, **nz; // 法線
  double **p, **q, **err;  // p, q, err
  double **dep;            // 奥行き写像

  int l0, l1, l2; // 画像用ループ変数
  int flag;

  double *w, **v, *x; // SVD用変数

  // メモリ確保
  // ベクトル，行列
  i  = (double *)dvector(1,inum);
  I  = (double **)dmatrix(1,inum,1,inum);
  Ii = (double *)dvector(1,inum);
  IV = (double **)dmatrix(1,inum,1,DIM);
  // SVD用
  w = (double *)dvector(1,DIM);
  v = (double **)dmatrix(1,DIM,1,DIM);
  x = (double *)dvector(1,DIM);

  // アルベド，法線，p，q，誤差，奥行き写像
  g  = (double **)dmatrix(0,ysize-1,0,xsize-1);
  nx = (double **)dmatrix(0,ysize-1,0,xsize-1);
  ny = (double **)dmatrix(0,ysize-1,0,xsize-1);
  nz = (double **)dmatrix(0,ysize-1,0,xsize-1);
  p  = (double **)dmatrix(0,ysize-1,0,xsize-1);
  q  = (double **)dmatrix(0,ysize-1,0,xsize-1);
  err= (double **)dmatrix(0,ysize-1,0,xsize-1);
  dep= (double **)dmatrix(0,ysize-1,0,xsize-1);

  // 画像配列の各々の点について処理する
  for (l1 = 0; l1 < ysize; l1++){
    for (l0 = 0; l0 < xsize; l0++){

      // 画素(l0,l1)の輝度がすべての画像で0でなかったら計算する
      flag = FALSE;
      for (l2 = 0; l2 < inum; l2++){
				if (iptr[l2][l1][l0] != 0){
					flag = TRUE;
					break;
				}
      }

      if (flag == TRUE){
				// 各画像の画素値からベクトルiと対角行列Iを構築する
				make_I_i(I, i, iptr, inum, l0, l1);
						
				// IVg=Iiを解き，この点でのアルベドgを計算し，それを使って法線nx，ny，nzとp，qを求める
				calc_Ii(Ii, I, i, inum);
				calc_IV(IV, I, V, inum);
				// my:IV, Iiが入力, xが出力？
				calc_svd(IV, w, v, Ii, x, inum);
				// my:xが出力?
				calc_vals(&g[l1][l0], x, &nx[l1][l0], &ny[l1][l0], &nz[l1][l0], &p[l1][l0], &q[l1][l0]);
      }

    }
  }

  safety_check(g, nx, ny, nz, err, p, q, xsize, ysize); // 安全性のチェック

  make_depth(dep, p, q, xsize, ysize);                  // 奥行き写像の計算

  out_all(g, nx, ny, nz, dep, err, xsize, ysize);       // 結果の出力

  // メモリ解放
  free_dmatrix(V,1,inum,1,DIM);
  free_dvector(i,1,inum);
  free_dmatrix(I,1,inum,1,inum);
  free_dvector(Ii,1,inum);
  free_dmatrix(IV,1,inum,1,DIM);
  free_dvector(w,1,DIM);
  free_dmatrix(v,1,DIM,1,DIM);
  free_dvector(x,1,DIM);

  free_dmatrix(g,0,ysize-1,0,xsize-1);
  free_dmatrix(nx,0,ysize-1,0,xsize-1);
  free_dmatrix(ny,0,ysize-1,0,xsize-1);
  free_dmatrix(nz,0,ysize-1,0,xsize-1);
  free_dmatrix(p,0,ysize-1,0,xsize-1);
  free_dmatrix(q,0,ysize-1,0,xsize-1);
  free_dmatrix(err,0,ysize-1,0,xsize-1);
  free_dmatrix(dep,0,ysize-1,0,xsize-1);

  return(0);
}

int set_light(char *inp_fname, double *V)
{
  int i, j;
  double tmp; // 一時変数

  for (i = 0; i < LIGHTNUM; i++){
    if (strncmp(inp_fname+strlen(inp_fname)-7, light[i].angle, ANGLE) == SAME){ // 後ろから7文字目からANGLE文字比較
      V[1] = light[i].x;
      V[2] = light[i].y;
      V[3] = 1.; // z
      // 単位ベクトルにする
      tmp = 0.;
      for (j = 1; j <= DIM; j++){
	      tmp += sqr(V[j]);
      }
      for (j = 1; j <= DIM; j++){
      	V[j] /= sqrt(tmp);
      }
      return(0);
    }
  }
  fprintf(stderr, "please check image file %s\n", inp_fname);
  exit(1);
}

int make_I_i(double **I, double *i, uchar **iptr[], int inum, int l0, int l1)
{
  int l2, l3;

  // 初期化
  for (l2 = 1; l2 <= inum; l2++){
    for (l3 = 1; l3 <= inum; l3++){
      I[l2][l3] = 0.;
    }
  }

  for (l2 = 1; l2 <= inum; l2++){
    I[l2][l2] = i[l2] = (double)iptr[l2-1][l1][l0]/(double)IMGMAX;
  }

  return(0);
}

int calc_Ii(double *Ii, double **I, double *i, int inum)
{
  int l2, l3;

  // Iiの計算
  for (l2 = 1; l2 <= inum; l2++){
    Ii[l2] = 0.;
    for (l3 = 1; l3 <= inum; l3++){
      Ii[l2] += I[l2][l3]*i[l3];
    }
  }

  return(0);
}

int calc_IV(double **IV, double **I, double **V, int inum)
{
  int l2, l3, l4;

  // IVの計算
  for (l2 = 1; l2 <= inum; l2++){
    for (l3 = 1; l3 <= DIM; l3++){
      IV[l2][l3] = 0.;
      for (l4 = 1; l4 <= inum; l4++){
      	IV[l2][l3] += I[l2][l4]*V[l4][l3];
      }
    }
  }

  return(0);
}

int calc_svd(double **IV, double *w, double **v, double *Ii, double *x, int inum)
{
  double wmax, wmin;
  int l3;

  // SVD
  svdcmp(IV, inum, DIM, w, v);
  wmax = 0.;
  for (l3 = 1; l3 <= DIM; l3++){
    if (w[l3] > wmax){
      wmax = w[l3];
    }
  }
  wmin = wmax * 1.0e-6;
  for (l3 = 1; l3 <= DIM; l3++){
    if (w[l3] < wmin) w[l3] = 0.;
  }
  svbksb(IV, w, v, inum, DIM, Ii, x);

  return(0);
}

int calc_vals(double *g, double *x, double *nx, double *ny, double *nz, double *p, double *q)
{
  // アルベドは|g|，法線はg/|g|，pはN1/N3，qはN2/N3
  // この中身を作成
  // この引数のxはro(x,y)*N(x,y).これの絶対値を取るとgになる
  // xは3次元のベクトル
  //アルベドを求める（xのノルムを計算）
  *g = sqrt( x[1]*x[1] + x[2]*x[2] + x[3]*x[3] );
  //法線を計算
  *nx = x[1] / *g;
  *ny = x[2] / *g;
  *nz = x[3] / *g;

  *p = *nx / *nz;
  *q = *ny / *nz;

  return(0);
}

int safety_check(double **g, double **nx, double **ny, double **nz, double **err, double **p, double **q, int xsize, int ysize)
{
  int l0, l1; // 画像用ループ変数
  int cnt = 0;

  int py=0;
  int qx=0;

  // 安全性のチェック
  // この中身を作成
  // pとqの偏微分を取る。偏微分は隣り合う画素との差を考える。
  // 端は無視する。 
  for(l1=1;l1<ysize-1;l1++)
    for(l0=1;l0<xsize-1;l0++) {
      //pのy方向の偏微分
      py = (p[l1+1][l0] - p[l1-1][l0]) / 2;
      //qのx方向の偏微分
      qx = (q[l1][l0+1] - q[l1][l0-1]) / 2;

      int check = ( py - qx ) * ( py - qx );
      if ( check > THRESHOLD ) {
			  err[l1][l0] = check;
			  cnt++;
		  }
    }
  
	printf("Safety check count: %d\n", cnt);

  return(0);
}

int make_depth(double **dep, double **p, double **q, int xsize, int ysize)
{
  int l0, l1, l2; // 画像用ループ変数
  double **tmp[SNUM]; // 一時変数，tmp[0]: 左上，tmp[1]: 右上，tmp[2]: 右下，tmp[3]: 左下

  // 奥行き写像の計算
  // できれば画像の4角それぞれを起点として計算した奥行きの平均を用いる
  // （配布資料では左上を起点としたときの手順が書かれている）
  // この中身を作成

  for (l2=0;l2<SNUM;l2++) {
    tmp[l2] = (double **)dmatrix(0,ysize-1,0,xsize-1);
  }
  //初期化
  for (l1=0;l1<ysize;l1++){
    for(l0=0;l0<xsize;l0++){
      dep[l1][l0] = 0;
      for( l2=0; l2<SNUM;l2++){
        tmp[l2][l1][l0] = 0;
      }
    }
  }

  //左上起点
  for(l1=1;l1<ysize-1;l1++) {
    //左上から左下へのy方向の積分を1行求める
    if (l1==1) {
      tmp[0][l1][0] += q[l1][0];
    } else if (l1 > 1) {
      tmp[0][l1][0] += q[l1][0] + tmp[0][l1-1][0];
    }
    
    for(l0=1;l0<xsize-1;l0++) {
      //x方向の積分（上のy方向の値を使う）
      if (l0 > 1) {
        tmp[0][l1][l0] += p[l1][l0] + tmp[0][l1][l0-1];
      }
      dep[l1][l0] = tmp[0][l1][0] + tmp[0][l1][l0];
    }
  }

  /*
  //左下
  for(l0=0;l0<xsize-1;l0++){
    //x方向の積分(左下から右下へ)
    if (l1==1) {
      tmp[1][l0][0] += q[l0][0];
    } else if (l1 > 1) {
      tmp[1][l0][0] += q[l0][0] + tmp[0][l1-1][0];
    }
  }
  */
  

  return(0);
}

int out_all(double **g, double **nx, double **ny, double **nz, double **dep, double **err, int xsize, int ysize)
{
  FILE *albedo, *normalVec, *depth, *depthInv, *error, *colorMap;
  int l0, l1; // 画像用ループ変数

  if (((albedo = fopen("albedo","w"))==NULL) || ((normalVec = fopen("normalVec","w"))==NULL)
      || ((depth = fopen("depth","w"))==NULL) || ((depthInv = fopen("depthInv","w"))==NULL) 
      || ((error = fopen("error","w"))==NULL) || ((colorMap = fopen("colorMap.ppm","w"))==NULL)){
    fprintf(stderr,"cannot open file\n");
    exit(1);
  }

  for (l1 = 0; l1 < ysize; l1=l1+ysize/32){
    for (l0 = 0; l0 < xsize; l0=l0+xsize/32){

//  for (l1 = 0; l1 < ysize; l1=l1++){
//    for (l0 = 0; l0 < xsize; l0=l0++){

      fprintf(albedo, "%d %d %g\n", l0, l1, g[l1][l0]);
      fprintf(normalVec, "%d %d %d %g %g %g\n", l0, l1, 0, nx[l1][l0], ny[l1][l0], nz[l1][l0]);
      fprintf(depth, "%d %d %g\n", l0, l1, dep[l1][l0]);
      fprintf(depthInv, "%d %d %g\n", l0, l1, -dep[l1][l0]);
      fprintf(error, "%d %d %g\n", l0, l1, err[l1][l0]);

    }
  }
// カラーマップ
  fprintf(colorMap, "P3\n%d %d\n%d\n", xsize, ysize, IMGMAX);
  for (l1 = 0; l1 < ysize; l1=l1++){
    for (l0 = 0; l0 < xsize; l0=l0++){
      fprintf(colorMap, "%d %d %d\n",
              (int)(((nx[l1][l0]+1.0)/2.0)*(double)(IMGMAX)),
              (int)(((ny[l1][l0]+1.0)/2.0)*(double)(IMGMAX)),
              (int)(((nz[l1][l0]+1.0)/2.0)*(double)(IMGMAX)));
    }
  }

  fclose(albedo);
  fclose(normalVec);
  fclose(depth);
  fclose(depthInv);
  fclose(error);
  fclose(colorMap);

  return(0);
}
