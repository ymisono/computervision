#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "nrutil.h"

/**** Typedefs ****/
typedef unsigned char uchar;

/**** Defines ****/
#define IHLENGTH     512 // 画像ヘッダの一行のサイズ
#define	IMGMAX	     255 // 画像の最高輝度値
#define	DIM	       3 // 次元数
#define THRESHOLD 1.0e-3 // 安全性チェックのしきい値
#define LIGHTNUM      17 // 光源の数
#define ANGLE          3 // 角度表現の桁数
#define FALSE         (0)
#define TRUE          (!FALSE)
#define SNUM           4 // 画像の辺数
#define SAME          (0)

/**** Macros ****/
#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

typedef struct _Light
{
  char *angle; // 光源の名前(角度)
  double x;
  double y;
  // double z; // 常に1.
} Light;
