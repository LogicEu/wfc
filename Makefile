# wfc makefile

STD=-std=c99
WFLAGS=-Wall -Wextra -pedantic
OPT=-O2
IDIR=-I. -Iutopia -Ispxe
LIBS=imgtool
CC=gcc
NAME=wfc
SRC=main.c

LDIR=lib
IDIR += $(patsubst %,-I%/,$(LIBS))
LSTATIC=$(patsubst %,lib%.a,$(LIBS))
LPATHS=$(patsubst %,$(LDIR)/%,$(LSTATIC))
LFLAGS=$(patsubst %,-L%,$(LDIR))
LFLAGS += $(patsubst %,-l%,$(LIBS))
LFLAGS += -lz -lpng -ljpeg
GL=-lglfw

SCRIPT=build.sh

OS=$(shell uname -s)
ifeq ($(OS),Darwin)
    #OSFLAGS=-mmacos-version-min=10.10
    GL+=-framework OpenGL
else 
	OSFLAGS=-lm -lpthread -D_POSIX_C_SOURCE=199309L
    GL+=-lGL -lGLEW
endif

CFLAGS=$(STD) $(WFLAGS) $(OPT) $(IDIR)

$(NAME): $(LPATHS) $(SRC)
	$(CC) -o $@ $(SRC) $(CFLAGS) $(LFLAGS) $(GL) $(OSFLAGS)

$(LDIR)/$(LDIR)%.a: $(LDIR)%.a $(LDIR)
	mv $< $(LDIR)/

$(LDIR): 
	@[ -d $@ ] || mkdir $@ && echo "mkdir $@"

$(LDIR)%.a: %
	cd $^ && make && mv $@ ../

exe:
	$(CC) -o $(NAME) $(SRC) $(CFLAGS) $(LFLAGS) $(GL) $(OSFLAGS)

clean: $(SCRIPT)
	./$^ $@
