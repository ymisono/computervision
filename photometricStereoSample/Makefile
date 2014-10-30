PROGRAM =	${PN}
OBJS	=	${PN}.o\
                nrutil.o\
		imgio.o\
		svdcmp.o\
		pythag.o\
		svbksb.o
CC	=	gcc
LIBS	=	-lm

${PROGRAM}:	${OBJS}
	${CC} -O ${OBJS} -o $@ ${LIBS}

.c.o:	
	${CC} ${CFLAGS} -c $<

${OBJS}:	header.h

clean:
	- rm -f ${PROGRAM} ${OBJS} *~ albedo normalVec depth depthInv error colorMap.ppm
