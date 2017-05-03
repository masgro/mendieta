#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "peano.h"

void peano_hilbert()
{
  int i, Tid;
  type_real DomainFac = 1.0 / (cp.lbox*1.001) * (((peanokey) 1) << (BITS_PER_DIMENSION));
  peanokey *Min, MinTot;

  peanokey *Key;            /*!< a table used for storing Peano-Hilbert keys for particles */

  Key = (peanokey *) malloc(sizeof(peanokey)*cp.npart);
  
  Min = (peanokey *) malloc(sizeof(peanokey)*NTHREADS);
  MinTot = (((peanokey)1)<<(3*BITS_PER_DIMENSION));

  for(i=0;i<NTHREADS;i++) 
    Min[i] = MinTot;  

  #pragma omp parallel for num_threads(NTHREADS) schedule(static) \
  default(none) private(Tid) shared(Key,Min,P,cp,DomainFac)
  for(i=0;i<cp.npart;i++) 
  {
    Tid = omp_get_thread_num();
    Key[i] = peano_hilbert_key(P[i].Pos[0]*DomainFac,P[i].Pos[1]*DomainFac,P[i].Pos[2]*DomainFac,BITS_PER_DIMENSION);  
    Min[Tid] = Key[i] < Min[Tid] ? Key[i] : Min[Tid] ;
  }
  
  for(i=0;i<NTHREADS;i++) 
    MinTot = Min[i] < MinTot ? Min[i] : MinTot;

  #pragma omp parallel for num_threads(NTHREADS) schedule(static) \
  default(none) shared(Key,MinTot,cp)
  for(i=0;i<cp.npart;i++) 
    Key[i] -= MinTot;
 
  free(Min);
  
  peano_hilbert_order(Key);

  free(Key);

  return;
}



/*! This function puts the particles into Peano-Hilbert order by sorting them
 *  according to their keys. The latter half already been computed in the
 *  domain decomposition. Since gas particles need to stay at the beginning of
 *  the particle list, they are sorted as a separate block.
 */
void peano_hilbert_order(peanokey *Key)
{
  int i;
  int *Id;
  struct peano_hilbert_data *mp;

  printf("Comienza Peano-Hilbert...\n");

  mp = (struct peano_hilbert_data *) malloc(sizeof(struct peano_hilbert_data) * cp.npart);
  Id = (int *) malloc(sizeof(int) * cp.npart);

  #pragma omp parallel for num_threads(NTHREADS) default(none) \
  shared(Key,mp,cp)
  for(i=0;i<cp.npart;i++)
  {
    mp[i].index = i;
    mp[i].key = Key[i];
  }

  qsort(mp, cp.npart, sizeof(struct peano_hilbert_data), compare_key);

  #pragma omp parallel for num_threads(NTHREADS) default(none) \
  shared(Id,mp,cp)
  for(i=0;i<cp.npart;i++)
    Id[mp[i].index] = i;

  reorder_particles(Id);

  free(Id);
  free(mp);

  return;
}

/*! This function is a comparison kernel for sorting the Peano-Hilbert keys. */
int compare_key(const void *a, const void *b)
{
  if(((struct peano_hilbert_data *) a)->key < (((struct peano_hilbert_data *) b)->key))
    return -1;

  if(((struct peano_hilbert_data *) a)->key > (((struct peano_hilbert_data *) b)->key))
    return +1;

  return 0;
}

/*! This function brings the collisionless particles into the same order as
 *  the sorted keys. (The sort is first done only on the keys themselves and
 *  done directly on the particles in order to reduce the amount of data that
 *  needs to be moved in memory. Only once the order is established, the
 *  particles are rearranged, such that each particle has to be moved at most
 *  once.)
 */
void reorder_particles(int *Id)
{
  int i;
  struct particle_data Psave, Psource;
  int idsource, idsave, dest;

  for(i=0;i<cp.npart;i++)
  {
      if(Id[i] != i)
	    {

    	  Psource = P[i];
	      idsource = Id[i];
    	  dest = Id[i];

	      do
	      {
	         Psave = P[dest];
	         idsave = Id[dest];

	         P[dest] = Psource;
	         Id[dest] = idsource;

	         if(dest == i)  break;

	         Psource = Psave;
	         idsource = idsave;

	         dest = idsource;

        }while(1);

	    }
   }

   return;
}

static int quadrants[24][2][2][2] = {
  /* rotx=0, roty=0-3 */
  {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
  {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
  {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
  {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
  /* rotx=1, roty=0-3 */
  {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
  {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
  {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
  {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
  /* rotx=2, roty=0-3 */
  {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
  {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
  {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
  {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
  /* rotx=3, roty=0-3 */
  {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
  {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
  {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
  {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
  /* rotx=4, roty=0-3 */
  {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
  {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
  {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
  {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
  /* rotx=5, roty=0-3 */
  {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
  {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
  {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
  {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};


static int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
  11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static int rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
static int roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };

static int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };

/*! This function computes a Peano-Hilbert key for an integer triplet (x,y,z),
 *  with x,y,z in the range between 0 and 2^bits-1.
 */
peanokey peano_hilbert_key(int x, int y, int z, int bits)
{
  int i, quad, bitx, bity, bitz;
  int mask, rotation, rotx, roty, sense;
  peanokey key;


  mask = 1 << (bits - 1);
  key = 0;
  rotation = 0;
  sense = 1;


  for(i = 0; i < bits; i++, mask >>= 1)
  {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;

      quad = quadrants[rotation][bitx][bity][bitz];

      key <<= 3;
      key += (sense == 1) ? (quad) : (7 - quad);

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];

      while(rotx > 0)
	    {
        rotation = rotxmap_table[rotation];
    	  rotx--;
	    }

      while(roty > 0)
    	{
	      rotation = rotymap_table[rotation];
    	  roty--;
    	}
  }

  return key;
}



