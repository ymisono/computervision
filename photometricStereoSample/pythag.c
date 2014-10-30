#include <math.h>
#define NRANSI
#include "nrutil.h"

double pythag(double a, double b)
/*
 * $(a^2 + b^2)^\frac{1}{2}$を計算する．
 * 悪いアンダーフローやオーバーフローをしない．
 */
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
#undef NRANSI
