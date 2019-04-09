### MENDIETA OPTIONS ####
#EXTRAS += -DDEBUG       # it will check for anomalities
#EXTRAS += -DCOMPUTE_EP
#EXTRAS += -DGETPOSITIONS
#EXTRAS += -DGETINDEX
#EXTRAS += -DSUBBOXES


### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1.0    #Positions in Kpc/h
#EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DENERGIES
#EXTRAS += -DREASIGNA

### IDENFOF options #############
#EXTRAS += -DREADIDENFOF
#EXTRAS += -DNHALO=1     # set to the number of halo to analyse
#EXTRAS += -Dlimpiamelo  #
#EXTRAS += -DCOMPUTE_FOF_PROPERTIES

## IDENSUB options ##############
EXTRAS += -DIDENSUB
EXTRAS += -DASSIGN_CLOSEST_GROUP

## INCLUDES #####################
#INC += -I/usr/local/include/gsl/

#Intel CC
#OMPP:=kinst-ompp-papi
#DOMPP:=-DompP
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=4
CFLAGS := -Wall -O3 -fopenmp -g
GSLL   := -lgsl -lgslcblas
LIBS   := -lm
LIBS   += $(GSLL)

#CFLAGS += -fdump-tree-vect-blocks=mendieta.dump
#CFLAGS += -fdump-tree-pre=stderr
#CFLAGS += -fopt-info
#CFLAGS += -fopt-info-missed

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := octree.o leesnap.o grid.o deltas.o compute_prop.o \
				variables.o limpieza.o io.o iden.o

HEADERS := $(patsubst %.o,%.h,$(OBJS))

EXEC := mendieta.x

todo: $(EXEC)

%.o: %.c $(MAKEFILE) $(HEADERS)
	$(CC) $(INC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

mendieta.x: mendieta.o $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

prop.x: prop.o propiedades.old.o octree.o leesnap.o grid.o deltas.o variables.o io.c
	$(CC) $(CFLAGS) $(EXTRAS) $^ -o $@ $(LIBS)

limpiador.x: $(OBJS) limpiador.o limpieza.o
	$(CC) $(CFLAGS) $(EXTRAS) $^ -o $@ $(GSLL)

clean:
	rm -rf $(OBJS)

cleanall: clean
	rm -rf $(EXEC)
