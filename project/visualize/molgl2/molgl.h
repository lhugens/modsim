#define MAXN 25000 /* maximum number of molecules to plot */
#define NUMBW 256
#define MAXAT 16 /* maximum number of aroms per molecule allowed (equal to the number of greys) */
#define MGL_USELIST
#define PI 2.0*acos(0.0)
#define TWOPI (4.0*acos(0.0))
#define PID2 acos(0.0)
#define STACKS 20
#define SLIDES 20
#define MGL_DISK_STACKS 15
#define MGL_DISK_SLIDES 15
#define MGL_NO_FADE  1
#define MGL_FADE_LIN 2
#define MGL_FADE_QUAD 3
#define MGL_MAX_FRAMES 4
#define MGL_REDBEFQUIT 5
#define SIGN(X) ((X>0)?1.0:(-1.0)) 
#define SIGN2(X)  ((X>0)?1.0:((X==0)?0.0:(-1.0)))
typedef struct {
	double x,y,z;
} XYZ;

#define Sqr(x) (x)*(x)
struct colStruct 
{
  float rgba[4];
  char name[32];
} *mgl_col;
enum atom_types {MGL_ATOM_SPHERE, MGL_ATOM_DISK, MGL_ATOM_CYLINDER, MGL_ATOM_SUPELLIPS, MGL_ATOM_SUPQUADRICS, MGL_ATOM_SPHERE_SPOT, MGL_ATOM_SPHERE_MSPOT};
typedef enum atom_types atom_types_e;

struct atom_common {
  atom_types_e type;
  double rx;
  double ry;
  double rz;
  int atcol;
  float atred;
  float atgreen;
  float atblue;
  int greyLvl;
  double transp;
};
struct atom_sphere
{
  struct atom_common common;
  double radius;
};
struct atom_supellips
{
  struct atom_common common;
  double R[3][3]; 
  double a; /* semi-axes */
  double b;
  double c;
  int n1;  /* integer for super-ellipsoid */
  int n2; 
  double tbeg;
  double tend;
};
struct atom_supquadrics
{
  struct atom_common common;
  double R[3][3]; 
  double a; /* semi-axes */
  double b;
  double c;
  double n1;  /* integer for super-ellipsoid */
  double n2; 
  double n3;
  double tbeg;
  double tend;
};
struct atom_sphere_spot
{
  struct atom_common common;
  double R[3][3]; 
  double a; /* semi-axes */
  double b;
  double c;
  int n1;  /* integer for super-ellipsoid */
  int n2; 
  double tbeg;
  double tend;
  int spotcol;
  float red; /* 0.0-1.0 */
  float green;
  float blue;
  double nx;
  double ny;
  double nz;
};
struct spotlst
{
  double n[3];
  int spotcol;
  float red;
  float green;
  float blue;
  double spotangle;
  struct spotlst *next;
};
struct atom_sphere_mspot
{
  struct atom_common common;
  double R[3][3]; 
  double a; /* semi-axes */
  double b;
  double c;
  int n1;  /* integer for super-ellipsoid */
  int n2; 
  double tbeg;
  double tend;
  int nspot;
  struct spotlst* sl;
};
struct atom_disk
{
  struct atom_common common;
  double nx; 
  double ny; 
  double nz;
  double toprad;
  double radius;
  double height;
};
struct atom_cylinder
{
  struct atom_common common;
  double nx; 
  double ny; 
  double nz;
  double toprad;
  double botrad;
  double height;
};
union atom 
{
  struct atom_common common; 
  struct atom_disk disk;
  struct atom_sphere sphere;
  struct atom_cylinder cylinder;
  struct atom_supellips supellips;
  struct atom_supquadrics supquadrics;
  struct atom_sphere_spot sphere_spot;
  struct atom_sphere_mspot sphere_mspot; 
};
enum bond_types {NONE, MGL_BOND_WIRE, MGL_BOND_CYLINDER};
typedef enum bond_types bond_types_e;
typedef union atom atom_u;
struct bond 
{
  /* indici relativi ai due atomi tra cui ci deve essere un bond */
  int from;
  int to;
  /* tipo del bond: ossia come deve essere rappresentato graficamente */
  bond_types_e type;
  double thickness; /* spessore del bond */
  int color; /* colore del bond */ 
  float red;
  float green;
  float blue;
  double transp;
};
typedef union atom atom_s; 
typedef struct bond bond_s;
struct molecule 
{
  atom_s *atom;
  bond_s *bond; 
  int numbonds;
  int numat;
};
struct global_settings
{
  int saved_counter;
  /*int saveimage=0, savedimg=0;*/
#ifdef MGL_USELIST
  int mgl_uselist;
#endif
  int saveandquit;
  char *savefile;
  int drawcube;
  double *sig;
  double* a;
  double* b;
  double* c;
  /*double *height;*/
  int numAt;
  int setdiameter;
  int setsemiax;
  double sa;
  double sb;
  double sc;
  int setheight;
  int *NumMols;
  int stacks;
  int slides;
  int Width;
  int Height;
  int infos;
  int frameNo; 
  int fadeMode; /* default = linear fading */
  int NA; 
  int axon;/* perspective by default */
  int bw;
  double L; /* lenght of an edge of the cubic box */
  float degx;
  float degy;
  float degz;
  float deginc;
  int *colIdxCol;
  int *colIdxBw;
  double viewangle;
  double near;
  double far;
  int setvp;
  double diameter;
  double height;
  double  dist;
  double ivpx;
  double ivpy;
  double ivpz;
  int default_bw;
  int default_col;
  double defbondthick;
  int defbondcol;
  double defbondtransp;
  double deftransp;
  int nrefresh;
  int twolights;
  int depthmask;
  int exitDelay;
  GLfloat light_pos0[4];
  GLfloat light_pos1[4];
};

