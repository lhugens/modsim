/*
   Revision of
   Carlos Vega & Santiago Lago
   Computers Chem. 18, 55-59, 1994

   Subrutine to evaluate the shortest distance between two rods of
   different length

   The original code did not give the symmetry property of the distance for almost parallel rods.
   The coordinates of the centers of the rods should be given in a periodic system

   r1,r2: centers of rods
   w1,w2: unit orientation vectors of rods
   lh1,lh2: halves of the length of rods
   Lv.x,Lv.y,Lv.z the edges of the periodic simulation cell
   */

#include <math.h>

//----------------- VECTOR operations: -----------------------------------------------------

#define VECT_COMMA ,
#define VECT_PAR (
#define VECT_PSEQ(_,SEP) (_ x)) SEP (_ y)) SEP (_ z))

#define VECT_COMP(x) .x
#define VECT_OP(A,COMP,OP,x) A COMP(x) OP
#define VECT_A_OP_B(A,OP,B,x) VECT_OP(A,VECT_COMP,OP,x) VECT_OP(B,VECT_COMP,,x)

#define VECT_OSEQ_(A,OP,B,SEP,_) \
    VECT_PSEQ(VECT_A_OP_B VECT_PAR A VECT_COMMA OP VECT_COMMA B VECT_COMMA,SEP##_)

#define VECT_OSEQ(A,OP,B,SEP) VECT_OSEQ_(A,OP,B,SEP,)
#define VECT_PROD(A,B) VECT_OSEQ(A,*,B,+)  /* product of A and B */
#define VECT_NORM2(A) VECT_PROD(A,A)  /* square of the norm of A */

#define VECT_OLIST(A,OP,B) VECT_OSEQ_(A,OP,B,VECT_COMMA,) /* (A.x OP B.x), ... */

#define VECT_SEQ(V,SEP) V(x) SEP V(y) SEP V(z)  /* because of the single macro expansion */
#define VECT_LIST(V) VECT_SEQ(V,VECT_COMMA)  /* V(x), ... */

typedef struct { double VECT_LIST(); } coo_t;

//---------------------------------------------------------------------------------------

extern coo_t Lv;


// Minimum distance in the periodic system:

#define MIN_RIJ(x) \
    ( FX= fabs(rij.x),(FX<Lv.x-FX)?rij.x:(rij.x-((rij.x >0)?Lv.x:-Lv.x) ) )


#define PW2(x) (x*x)

static inline double sign(double a,double b) { return a= fabs(a),(b<0)?-a:a; }

//---------------- Distance of two rods: -------------------------------------

double dist2_rods(coo_t r1,coo_t r2,coo_t w1,coo_t w2,double lh1,double lh2)
{
    coo_t rij= { VECT_OLIST(r2,-,r1) };
    double FX;
    coo_t min_rij= { VECT_LIST(MIN_RIJ) };
    double
        xla,xmu,
        rr= VECT_NORM2(min_rij),
        rw1= VECT_PROD(min_rij,w1),
        rw2= VECT_PROD(min_rij,w2),
        w1w2= VECT_PROD(w1,w2),
        cc= 1-PW2(w1w2);

    // Checking whether the rods are or not parallel:
    // The original code is modified to have symmetry:

    if(cc<1e-6) {
        if(rw1 && rw2) {
            xla= rw1/2;
            xmu= -rw2/2;
        }
        else return rr;
    }

    else {

        // Step 1

        xla= (rw1-w1w2*rw2)/cc;
        xmu= (-rw2+w1w2*rw1)/cc;
    }

    // Step 2

    if( fabs(xla)>lh1 || fabs(xmu)>lh2 ) {

        // Step 3 - 7

        if(fabs(xla)-lh1>fabs(xmu)-lh2) {
            xla= sign(lh1,xla);
            xmu= xla*w1w2-rw2;
            if( fabs(xmu)>lh2 ) xmu= sign(lh2,xmu);
        }
        else {
            xmu= sign(lh2,xmu);
            xla= xmu*w1w2+rw1;
            if( fabs(xla)>lh1 ) xla= sign(lh1,xla);
        }
    }

    // Step 8

    return rr+PW2(xla)+PW2(xmu) + 2*(xmu*rw2 -xla*(rw1+xmu*w1w2));
}

int main(){
}
