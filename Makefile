### MENDIETA OPTIONS ####
#EXTRAS += -DDEBUG       # it will check for anomalities
#EXTRAS += -DCOMPUTE_EP
#EXTRAS += -DGETPOSITIONS
#EXTRAS += -DSUBBOXES


### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0    #Positions in Kpc/h
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DENERGIES

### IDENFOF options #############
#EXTRAS += -DREADIDENFOF
#EXTRAS += -DNHALO=1     # set to the number of halo to analyse
EXTRAS += -Dlimpiamelo  #
#EXTRAS += -DCOMPUTE_FOF_PROPERTIES

## IDENSUB options ##############
EXTRAS += -DIDENSUB
#EXTRAS += -DASSIGN_CLOSEST_GROUP
#EXTRAS += -DPROP
#EXTRAS += -I/usr/local/include/gsl/

#Intel CC
#OMPP:=kinst-ompp-papi
#DOMPP:=-DompP
#CC     := $(OMPP) icc $(DOMPP)
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=6
CFLAGS := -Wall -O3 -fopenmp -g
#GSLL   := -L/usr/local/lib/ -lgsl -lgslcblas
GSLL   := -lgsl -lgslcblas
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := octree.o leesnap.o grid.o deltas.o compute_prop.o peano.o \
				variables.o limpieza.o io.o iden.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := mendieta.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

mendieta.x: mendieta.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $^  -o $@ $(GSLL) $(LIBS)

prop.x: $(OBJS) prop.o propiedades.o
	$(CC) $(CFLAGS) $(EXTRAS) $^ -o $@ $(GSLL)

clean:
	rm -rf $(OBJS)
	rm -rf mendieta.o

cleanall: clean
	rm -rf $(EXEC)
