#ifndef TOOLS3B_H
#define TOOLS3B_H

#define LOOP3B( NPARTICLES,  BODY )	\
{ \
  int kj=0;  \
  int ki2=0; \
  int ki = 0; \
  for(int k=0;k<NPARTICLES;k++) \
    { \
      int ji=0;\
      for (int j=0;j<k;j++ & kj++) \
	{ \
	 ki=ki2; \
	  for(int i=0;i<j;i++ & ji++ & ki++) \
	    { \
              BODY \
	    }  \
	}\
      ki2+=k; \
      }	      \
} \

#define LOOP3B_DIS( N1, N2 , N3,  BODY )		\
{ \
  int ij=0;  \
  int ik2=0; \
  int ik = 0;		\
  for(int i=0;i<N1;i++) \
    { \
      int jk=0;\
      for (int j=0;j<N2;j++ & ij++) \
	{ \
	 ik=ik2; \
	  for(int k=0;k<N3;k++ & jk++ & ik++) \
	    { \
              BODY \
	    }  \
	}\
      ik2+=N3; \
      }	      \
} \


#endif
