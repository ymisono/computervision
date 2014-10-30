#include "header.h"

uchar **readimg(char *iFname, char *magic, int *xsize, int *ysize, int *maxval)
{
  FILE *fp;
  uchar **iptr;
  int i, j;
  char buf[IHLENGTH];

  // 入力画像の読み込み
  if ((fp = fopen(iFname, "rb")) == NULL) exit(1);

  fgets( magic, 256, fp );
  fgets( buf, 256, fp );
  while ((int)buf[0] == 35){ // #の行捨てる
    fgets( buf, 256, fp );
  }
  sscanf(buf, "%d %d\n", xsize, ysize );
  fgets( buf, 256, fp );
  sscanf(buf, "%d\n", maxval );

  // 入力画像領域
  iptr = cmatrix(0, *xsize-1, 0, *ysize-1);

  // 入力画像読み込み
  for (j = 0; j < *ysize; j++){
    fread(iptr[j], sizeof(uchar), *xsize, fp);
  } 
  fclose(fp);

  return( iptr );
}

int writeimg(uchar **optr, char *oFname, char *magic, int xsize, int ysize, int maxval)
{
  FILE *fp;
  int j;
  char buf[IHLENGTH];

  // 出力画像の吐き出し
  if ((fp = fopen(oFname, "wb")) == NULL) return(-1);

  fputs( magic, fp );
  sprintf( buf, "%d %d\n", xsize, ysize);
  fputs( buf, fp );
  sprintf( buf, "%d\n", maxval);
  fputs( buf, fp );

  for (j = 0; j < ysize; j++){
    fwrite(optr[j], sizeof(uchar), xsize, fp);
  }

  fclose(fp);

  return(0);
}
