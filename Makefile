#CC = gcc
CFLAGS = -O4
#CFLAGS = -g 
#CFLAGS = -O4

LIB = -lm

OBJS=	\
		main.o \
		fast.o \
		kmeres.o

SRCS= $(OBJS:.o=.c)

INCS=	\
		fasta.h \
		ssaha.h

PROGRAM = kmeres 

$(PROGRAM): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	rm -f $(PROGRAM) *.o
