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
#EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DENERGIES
#EXTRAS += -DREASIGNA

### IDENFOF options #############
#EXTRAS += -DREADIDENFOF
#EXTRAS += -DNHALO=1     # set to the number of halo to analyse
#EXTRAS += -Dlimpiamelo  #
#EXTRAS += -DCOMPUTE_PROPERTIES

## IDENSUB options ##############
EXTRAS += -DIDENSUB
EXTRAS += -DASSIGN_CLOSEST_GROUP

## INCLUDES #####################
INC += -I/usr/local/include/gsl/

#Intel CC
#OMPP:=kinst-ompp-papi
#DOMPP:=-DompP
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=8
CFLAGS := -Wall -O3 -fopenmp
GSLL   := -L/usr/local/lib/ -lgsl -lgslcblas -lm
LIBS   := -lm

#CFLAGS += -fdump-tree-vect-blocks=mendieta.dump
#CFLAGS += -fdump-tree-pre=stderr
#CFLAGS += -fopt-info
#CFLAGS += -fopt-info-missed

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := octree.o leesnap.o grid.o deltas.o \
				variables.o io.o iden.o limpieza.o # compute_prop.o

HEADERS := $(patsubst %.o,%.h,$(OBJS))

EXEC := mendieta.x

todo: $(EXEC)

%.o: %.c $(MAKEFILE) $(HEADERS)
	$(CC) $(INC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

mendieta.x: mendieta.o $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(GSLL) $(LIBS)

prop.x: prop.o propiedades.o $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(GSLL)

clean:
	rm -rf $(OBJS)
	rm -rf $(EXEC)

cleanall: clean
	rm -rf $(EXEC)
