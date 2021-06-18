#ifdef MGL_MACOS
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "molgl.h"
#include <zlib.h>
#define SQ_CALC_NORM
#define SQ_REND_SYM
/* NOTA 27/04/2010: DELSQ*TWOPI/n1 e DELSQ*TWOPI/n2 sono i passi in radianti per stimare numericamente 
   il gradiente delle superquadriche */
#define DELSQ 1E-5

const int NUMCOLS = 746;
float mgl_bw[NUMBW][4];
char *mglrgb[]={
#include "mglrgb.h"
};
struct molecule **mols = NULL;
struct global_settings globset;

/* array con i valori di default */
char inputFile[512];
GLuint *atomsList = NULL;

/* mgl_*[colIdx*[j]] is the color of the j-th atom */ 
void setLight(void)
{
  GLfloat light_ambient0[] = { 1.0, 1.0, 1.0, 0.15 };
  GLfloat light_diffuse0[] = { 1.0, 1.0, 1.0, 0.15 };
  GLfloat light_specular0[] = { 1.0, 1.0, 1.0, 0.15 };
  GLfloat light_ambient1[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light_diffuse1[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light_specular1[] = { 1.0, 1.0, 1.0, 1.0 }; 
 
  /*	light_position is NOT default value	*/
  //GLfloat light_position0[] = { 10.0, 10.0, 10.0, 0.0 };
  //GLfloat light_position1[] = { -10.0, -10.0, -10.0, 0.0};
  GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 0.15 };
  GLfloat local_view[] = { 0.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);
  
 
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

  glLightfv (GL_LIGHT0, GL_AMBIENT, light_ambient0);
  glLightfv (GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
  glLightfv (GL_LIGHT0, GL_SPECULAR, light_specular0);
  glLightfv (GL_LIGHT0, GL_POSITION, globset.light_pos0);
  glLightfv (GL_LIGHT1, GL_AMBIENT, light_ambient1);
  glLightfv (GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
  glLightfv (GL_LIGHT1, GL_SPECULAR, light_specular1);
  glLightfv (GL_LIGHT1, GL_POSITION, globset.light_pos1);
  // tuning light spot size
  //glLightf (GL_LIGHT0, GL_SPOT_CUTOFF, 30.f);
  //glLightf (GL_LIGHT1, GL_SPOT_CUTOFF, 30.f);
}

/* =============================== >>> setBW <<< ===========================*/
void setBW(void)
{
  int nc;
  for (nc = 0; nc < NUMBW; ++nc)
    {
      mgl_bw[nc][0] = ((float) nc) / 255.0;
      mgl_bw[nc][1] = ((float) nc) / 255.0;
      mgl_bw[nc][2] = ((float) nc) / 255.0;
      mgl_bw[nc][3] = 1.0; 
    }  
  /* 255 levels of gray */
}

/*  Initialize material property and light source. */
void myinit (void)
{
    setLight();
    //glFrontFace (GL_CW);
    glEnable (GL_LIGHTING);
    glEnable (GL_LIGHT0);
    if (globset.twolights)
      glEnable(GL_LIGHT1);
#if 1
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glEnable(GL_CULL_FACE);
#else
    glEnable(GL_POLYGON_STIPPLE);
#endif
    //glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
#if 0
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
#endif
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glClearColor(1.0, 1.0, 1.0, 0.0);
}

/* ========================== >>> setColor <<< =============================*/
void setColor(float col[4], double ff)
{
  /* col is the specular color, diffusie and the ambient are calculated
     scaling this by df and af respectively */
  float mat[4]; /* atoms color */
  // OLD VALURE PRE 29/01/2010: float df = 0.95, af = 0.2, sf=1.0;
  float df = 0.85, af = 0.4, sf=0.85;
#if 0
  int i;
#endif
#if 0
    mat[0] = ff*af*col[0]; mat[1] = ff*af*col[1]; 
    mat[2] = ff*af*col[2]; mat[3] = ff*col[3];
    glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
    mat[0] = ff*df*col[0]; mat[1] = ff*df*col[1]; mat[2] = ff*df*col[2];	
    glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
    mat[0] = ff*col[0]; mat[1] = ff*col[1]; mat[2] = ff*col[2];
    glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
    glMaterialf (GL_FRONT, GL_SHININESS, 0.9*128.0);
#endif
  mat[0] = af*col[0]; mat[1] = af*col[1]; 
  mat[2] = af*col[2]; mat[3] = ff*col[3];
  glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
  mat[0] = df*col[0]; mat[1] = df*col[1]; mat[2] = df*col[2];	
  glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
  mat[0] = sf*col[0]; mat[1] = sf*col[1]; mat[2] = sf*col[2];
  glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
  //OLD// glMaterialf (GL_FRONT, GL_SHININESS, 0.9*128.0);
  glMaterialf (GL_FRONT, GL_SHININESS, 0.9*128.0);
}
/* ========================== >>> setColor <<< =============================*/
void setColorRGB(int ncol, double ff, float red, float green, float blue)
{
  float col[4];
  int i;
  if (ncol==-2)
    {
      col[0] = red;
      col[1] = green;
      col[2] = blue;
      col[3] = 1.0;
    }
  else
    {
      for (i=0; i < 4; i++)
	col[i] = mgl_col[ncol].rgba[i];
    }
  //printf("col=%f %f %f %f\n", col[0], col[1], col[2], col[3]);
  setColor(col, ff);
}
/* ======================== >>> calcFadeFact <<< ===========================*/
double  calcFadeFact(int mode, int nf)
{
  double ff;
  if (mode == MGL_FADE_LIN)
    {
      ff =  1.0 - ((double)nf) / ((double)globset.frameNo) ;
    }
  else if(mode == MGL_FADE_QUAD)
    {
      ff = 1.0 / Sqr(((double)nf) / ((double)globset.frameNo));
    }
  else ff = 1;

  return ff;
}
void vectProd(double r1x, double r1y, double r1z, 
	 double r2x, double r2y, double r2z, 
	 double* r3x, double* r3y, double* r3z)
{
  /* DESCRIPTIOM:
     r3 = [ r1, r2 ] where [ , ] the vectorial product */
  *r3x = r1y * r2z - r1z * r2y; 
  *r3y = r1z * r2x - r1x * r2z;
  *r3z = r1x * r2y - r1y * r2x;
}
/*
   Create a superellipse
   "method" is 0 for quads, 1 for triangles
      (quads look nicer in wireframe mode)/
   This is a "unit" ellipsoid (-1 to 1) at the origin,
      glTranslate() and glScale() as required.
*/
void EvalSuperEllipse(double t1,double t2,double p1,double p2,
		      double a, double b, double c, XYZ *p);
void Normalise(XYZ *p)
{
   double length;

   length = p->x * p->x + p->y * p->y + p->z * p->z;
   if (length > 0) 
     {
       length = sqrt(length);
       p->x /= length;
       p->y /= length;
       p->z /= length;
     }
   else
     {
       p->x = 0;
       p->y = 0;
       p->z = 0;
     }	
}

XYZ CalcNormal(XYZ p, XYZ p1, XYZ p2)
{
  XYZ n,pa,pb;

   pa.x = p1.x - p.x;
   pa.y = p1.y - p.y;
   pa.z = p1.z - p.z;
   pb.x = p2.x - p.x;
   pb.y = p2.y - p.y;
   pb.z = p2.z - p.z;
   n.x = pa.y * pb.z - pa.z * pb.y;
   n.y = pa.z * pb.x - pa.x * pb.z;
   n.z = pa.x * pb.y - pa.y * pb.x;
   n.x = - n.x;
   n.y = - n.y;
   n.z = - n.z;
   Normalise(&n);

   return(n);
}
void EvalSuperQuadrics(double t1,double t2,double p1,double p2,double p3, double a, double b, double c, XYZ *p);
void EvalSuperQuadricsNorm(double t1,double t2,double p1,double p2,double p3, double a, double b, double c, XYZ *p, XYZ *en);

void CreatePartialSuperQuadrics(double power1,double power2, double power3, 
				double a, double b, double c,
			int n1, int n2, int method, double thetaBeg, double thetaEnd)
{
   int i,j;
   double theta1,theta2,theta3;
   XYZ p,p1,p2,en;
   int n1beg, n1end;
   double delta1, delta2;
   /* n1 = stacks
    * n2 = slides */
   /* Shall we just draw a point? */
#if 0
   if (n1 < 4 && n2 < 4) {
      glBegin(GL_POINTS);
      glVertex3f(0.0,0.0,0.0);
      glEnd();
      return;
   }

   /* Shall we just draw a plus */
   if (power1 > 10 && power2 > 10) {
      glBegin(GL_LINES);
      glVertex3f(-1.0, 0.0, 0.0);
      glVertex3f( 1.0, 0.0, 0.0);
      glVertex3f( 0.0,-1.0, 0.0);
      glVertex3f( 0.0, 1.0, 0.0);
      glVertex3f( 0.0, 0.0,-1.0);
      glVertex3f( 0.0, 0.0, 1.0);
      glEnd();
      return;
   }
#endif
   delta1 = DELSQ*TWOPI / (double)n1;
   delta2 =  DELSQ*TWOPI / (double)n2;
   //printf("boh=%.15G n1/2=%d thetaBeg =%.15G TWOPI=%.15G thetaBeg/TWOPI=%.15G\n", ((double)(n1/2))*thetaBeg/TWOPI, n1/2, thetaBeg, TWOPI, thetaBeg/TWOPI);
   if (thetaBeg > 0)
     n1beg = (int) rint(n1*thetaBeg/TWOPI);
   else 
     n1beg = 0;
   if (thetaEnd > 0)
     n1end = (int) rint(n1*thetaEnd/TWOPI);
   else
     n1end = n1/2;
   //printf("thetaBeg: %.15G thetaEnd: %.15G n1beg=%d n1end=%d n1=%d n2=%d\n", thetaBeg, thetaEnd, n1beg, n1end, n1, n2);
   for (j=0;j<n1/2;j++) {
      if (!(j >= n1beg && j < n1end))
	continue;
      theta1 = j * TWOPI / (double)n1 - PID2;
      theta2 = (j+1) * TWOPI / (double)n1 - PID2;
      if (method==2)
	glBegin(GL_TRIANGLE_FAN);
      else if (method == 0)
         glBegin(GL_QUAD_STRIP);
      else
         glBegin(GL_TRIANGLE_STRIP);
      for (i=0;i<=n2;i++) {
         if (i == 0 || i == n2)
            theta3 = 0;
         else
            theta3 = i * TWOPI / n2;
   
#ifndef SQ_CALC_NORM
         EvalSuperQuadrics(theta2,theta3,power1,power2,power3,a,b,c,&p);
         EvalSuperQuadrics(theta2+delta1,theta3,power1,power2,power3,a,b,c,&p1);
         EvalSuperQuadrics(theta2,theta3+delta2,power1,power2,power3,a,b,c,&p2);
         en = CalcNormal(p,p1,p2);
#else
	 EvalSuperQuadricsNorm(theta2,theta3,power1,power2,power3,a,b,c,&p,&en);
#endif
	 glNormal3f(en.x,en.y,en.z);
         //glTexCoord2f(i/(double)n,2*(j+1)/(double)n);
	 //glColor4f(1,1,1,0.1);
         glVertex3f(p.x,p.y,p.z);

#ifndef SQ_CALC_NORM
         EvalSuperQuadrics(theta1,theta3,power1,power2,power3,a,b,c,&p);
         EvalSuperQuadrics(theta1+delta1,theta3,power1,power2,power3,a,b,c,&p1);
         EvalSuperQuadrics(theta1,theta3+delta2,power1,power2,power3,a,b,c,&p2);
         en = CalcNormal(p,p1,p2);
#else
	 EvalSuperQuadricsNorm(theta1,theta3,power1,power2,power3,a,b,c,&p,&en);
#endif
	 glNormal3f(en.x,en.y,en.z);
         //glTexCoord2f(i/(double)n,2*j/(double)n);
	 //glColor4f(1,1,1,0.1);
         glVertex3f(p.x,p.y,p.z);
      }
      glEnd();
   }
}

void CreateSuperQuadrics(double power1,double power2,double power3,double a, double b, double c,
			int n1, int n2, int method)
{
   int i,j;
   double theta1,theta2,theta3;
   XYZ p,p1,p2,en;
   double delta1, delta2;
   /* n1 = stacks
    * n2 = slides */
#if 0
   /* Shall we just draw a point? */
   if (n1 < 4 && n2 < 4) {
      glBegin(GL_POINTS);
      glVertex3f(0.0,0.0,0.0);
      glEnd();
      return;
   }

   /* Shall we just draw a plus */
   if (power1 > 10 && power2 > 10) {
      glBegin(GL_LINES);
      glVertex3f(-1.0, 0.0, 0.0);
      glVertex3f( 1.0, 0.0, 0.0);
      glVertex3f( 0.0,-1.0, 0.0);
      glVertex3f( 0.0, 1.0, 0.0);
      glVertex3f( 0.0, 0.0,-1.0);
      glVertex3f( 0.0, 0.0, 1.0);
      glEnd();
      return;
   }
#endif
   delta1 = DELSQ*TWOPI / (double)n1;
   delta2 =  DELSQ*TWOPI / (double)n2;
   //printf("n1=%d n2=%d\n", n1, n2);
   for (j=0;j<(n1/2);j++) {
      theta1 = j * TWOPI / (double)n1 - PID2;
      theta2 = (j+1) * TWOPI / (double)n1 - PID2;
      if (method==2)
	glBegin(GL_TRIANGLE_FAN);
      else if (method == 0)
      	glBegin(GL_QUAD_STRIP);
      else
	glBegin(GL_TRIANGLE_STRIP);

      for (i=0;i<=n2;i++) {
         if (i == 0 || i == n2)
            theta3 = 0;
         else
            theta3 = i * TWOPI / n2;// - TWOPI/2.0;
   
#ifndef SQ_CALC_NORM
         EvalSuperQuadrics(theta2,theta3,power1,power2,power3,a,b,c,&p);
         EvalSuperQuadrics(theta2+delta1,theta3,power1,power2,power3,a,b,c,&p1);
         EvalSuperQuadrics(theta2,theta3+delta2,power1,power2,power3,a,b,c,&p2);
     	 en = CalcNormal(p,p1,p2);
#else
         EvalSuperQuadricsNorm(theta2,theta3,power1,power2,power3,a,b,c,&p,&en);
#endif
         glNormal3f(en.x,en.y,en.z);
         //glTexCoord2f(i/(double)n,2*(j+1)/(double)n);
	 //glColor4f(1,1,1,0.1);
         glVertex3f(p.x,p.y,p.z);

#ifndef  SQ_CALC_NORM
         EvalSuperQuadrics(theta1,theta3,power1,power2,power3,a,b,c,&p);
         EvalSuperQuadrics(theta1+delta1,theta3,power1,power2,power3,a,b,c,&p1);
         EvalSuperQuadrics(theta1,theta3+delta2,power1,power2,power3,a,b,c,&p2);
     	 en = CalcNormal(p,p1,p2);
#else
         EvalSuperQuadricsNorm(theta1,theta3,power1,power2,power3,a,b,c,&p,&en);
#endif
     	 glNormal3f(en.x,en.y,en.z);
         //glTexCoord2f(i/(double)n,2*j/(double)n);
	 //glColor4f(1,1,1,0.1);
         glVertex3f(p.x,p.y,p.z);
      }
      glEnd();
   }
}
#if 1
void CreatePartialSuperEllipse(double power1,double power2, double a, double b, double c,
			int n1, int n2, int method, double thetaBeg, double thetaEnd)
{
   int i,j;
   double theta1,theta2,theta3;
   XYZ p,p1,p2,en;
   int n1beg, n1end;
   double delta1, delta2;
   /* n1 = stacks
    * n2 = slides */
   /* Shall we just draw a point? */
   if (n1 < 4 && n2 < 4) {
      glBegin(GL_POINTS);
      glVertex3f(0.0,0.0,0.0);
      glEnd();
      return;
   }

   /* Shall we just draw a plus */
   if (power1 > 10 && power2 > 10) {
      glBegin(GL_LINES);
      glVertex3f(-1.0, 0.0, 0.0);
      glVertex3f( 1.0, 0.0, 0.0);
      glVertex3f( 0.0,-1.0, 0.0);
      glVertex3f( 0.0, 1.0, 0.0);
      glVertex3f( 0.0, 0.0,-1.0);
      glVertex3f( 0.0, 0.0, 1.0);
      glEnd();
      return;
   }
   delta1 = 1E-6*TWOPI / (double)n1;
   delta2 =  1E-6*TWOPI / (double)n2;
   //printf("boh=%.15G n1/2=%d thetaBeg =%.15G TWOPI=%.15G thetaBeg/TWOPI=%.15G\n", ((double)(n1/2))*thetaBeg/TWOPI, n1/2, thetaBeg, TWOPI, thetaBeg/TWOPI);
   if (thetaBeg > 0)
     n1beg = (int) rint(n1*thetaBeg/TWOPI);
   else 
     n1beg = 0;
   if (thetaEnd > 0)
     n1end = (int) rint(n1*thetaEnd/TWOPI);
   else
     n1end = n1/2;
   //printf("thetaBeg: %.15G thetaEnd: %.15G n1beg=%d n1end=%d n1=%d n2=%d\n", thetaBeg, thetaEnd, n1beg, n1end, n1, n2);
   for (j=0;j<n1/2;j++) {
      if (!(j >= n1beg && j < n1end))
	continue;
      theta1 = (j+1) * TWOPI / (double)n1 - PID2;
      theta2 = j * TWOPI / (double)n1 - PID2;
      if (method==2)
	glBegin(GL_TRIANGLE_FAN);
      else if (method == 0)
         glBegin(GL_QUAD_STRIP);
      else
         glBegin(GL_TRIANGLE_STRIP);
      for (i=0;i<=n2;i++) {
         if (i == 0 || i == n2)
            theta3 = 0;
         else
            theta3 = i * TWOPI / n2;
   
         EvalSuperEllipse(theta2,theta3,power1,power2,a,b,c,&p);
         EvalSuperEllipse(theta2+delta1,theta3,power1,power2,a,b,c,&p1);
         EvalSuperEllipse(theta2,theta3+delta2,power1,power2,a,b,c,&p2);
         en = CalcNormal(p,p1,p2);
         glNormal3f(en.x,en.y,en.z);
         //glTexCoord2f(i/(double)n,2*(j+1)/(double)n);
	 //glColor4f(1,1,1,0.1);
         glVertex3f(p.x,p.y,p.z);

         EvalSuperEllipse(theta1,theta3,power1,power2,a,b,c,&p);
         EvalSuperEllipse(theta1+delta1,theta3,power1,power2,a,b,c,&p1);
         EvalSuperEllipse(theta1,theta3+delta2,power1,power2,a,b,c,&p2);
         en = CalcNormal(p,p1,p2);
         glNormal3f(en.x,en.y,en.z);
         //glTexCoord2f(i/(double)n,2*j/(double)n);
	 //glColor4f(1,1,1,0.1);
         glVertex3f(p.x,p.y,p.z);
      }
      glEnd();
   }
}
#endif

void CreateSuperEllipse(double power1,double power2, double a, double b, double c,
			int n1, int n2, int method)
{
   int i,j;
   double theta1,theta2,theta3;
   XYZ p,p1,p2,en;
   double delta1, delta2;
   /* n1 = stacks
    * n2 = slides */
   /* Shall we just draw a point? */
   if (n1 < 4 && n2 < 4) {
      glBegin(GL_POINTS);
      glVertex3f(0.0,0.0,0.0);
      glEnd();
      return;
   }

   /* Shall we just draw a plus */
   if (power1 > 10 && power2 > 10) {
      glBegin(GL_LINES);
      glVertex3f(-1.0, 0.0, 0.0);
      glVertex3f( 1.0, 0.0, 0.0);
      glVertex3f( 0.0,-1.0, 0.0);
      glVertex3f( 0.0, 1.0, 0.0);
      glVertex3f( 0.0, 0.0,-1.0);
      glVertex3f( 0.0, 0.0, 1.0);
      glEnd();
      return;
   }
   delta1 = 1E-6*TWOPI / (double)n1;
   delta2 =  1E-6*TWOPI / (double)n2;
   for (j=0;j<n1/2;j++) {
      theta1 = (j+1) * TWOPI / (double)n1 - PID2;
      theta2 = j * TWOPI / (double)n1 - PID2;
      if (method==2)
	glBegin(GL_TRIANGLE_FAN);
      else if (method == 0)
         glBegin(GL_QUAD_STRIP);
      else
         glBegin(GL_TRIANGLE_STRIP);
      for (i=0;i<=n2;i++) {
         if (i == 0 || i == n2)
            theta3 = 0;
         else
            theta3 = i * TWOPI / n2;
   
         EvalSuperEllipse(theta2,theta3,power1,power2,a,b,c,&p);
         EvalSuperEllipse(theta2+delta1,theta3,power1,power2,a,b,c,&p1);
         EvalSuperEllipse(theta2,theta3+delta2,power1,power2,a,b,c,&p2);
         en = CalcNormal(p,p1,p2);
         glNormal3f(en.x,en.y,en.z);
         //glTexCoord2f(i/(double)n,2*(j+1)/(double)n);
	 //glColor4f(1,1,1,0.1);
         glVertex3f(p.x,p.y,p.z);

         EvalSuperEllipse(theta1,theta3,power1,power2,a,b,c,&p);
         EvalSuperEllipse(theta1+delta1,theta3,power1,power2,a,b,c,&p1);
         EvalSuperEllipse(theta1,theta3+delta2,power1,power2,a,b,c,&p2);
         en = CalcNormal(p,p1,p2);
         glNormal3f(en.x,en.y,en.z);
         //glTexCoord2f(i/(double)n,2*j/(double)n);
	 //glColor4f(1,1,1,0.1);
         glVertex3f(p.x,p.y,p.z);
      }
      glEnd();
   }
}
void EvalSuperEllipse(double t1,double t2,double p1,double p2,
		      double a, double b, double c, XYZ *p)
{
   double tmp;
   double ct1,ct2,st1,st2;

   ct1 = cos(t1);
   ct2 = cos(t2);
   st1 = sin(t1);
   st2 = sin(t2);

   tmp  = SIGN(ct1) * pow(fabs(ct1),p1);
   p->x = a * tmp * SIGN(ct2) * pow(fabs(ct2),p2);
   p->y = b * SIGN(st1) * pow(fabs(st1),p1);
   p->z = c * tmp * SIGN(st2) * pow(fabs(st2),p2);
}
void swap_axes(double *a, double *b, double *c, double *p1, double *p2, double *p3);

void EvalSuperQuadricsNorm(double t1,double t2,double p1,double p2,double p3,
		      double a, double b, double c, XYZ *p, XYZ *en)
{
  double ct1,ct2,st1,st2;
  
  ct1 = cos(t1);
  ct2 = cos(t2);
  st1 = sin(t1);
  st2 = sin(t2);
  //p1=p2=p3=2.0;
#ifdef SQ_REND_SYM
  swap_axes(&a,&b,&c,&p1,&p2,&p3);
#endif
  p->x = a*SIGN2(ct1) * pow(fabs(ct1),2.0/p1) * SIGN2(ct2) * pow(fabs(ct2),2.0/p1);
  //printf("x actual=%.15G ellips=%.15G\n", p->x, a*SIGN2(ct1) * pow(fabs(ct1),1.0/2) * SIGN2(ct2) * pow(fabs(ct2),1.0/2));
  p->y = b*SIGN2(ct1) * pow(fabs(ct1),2.0/p2) * SIGN2(st2) * pow(fabs(st2),2.0/p2);
  //printf("y actual=%.15G ellips=%.15G\n", p->y, b*SIGN2(ct1) * pow(fabs(ct1),1.0/2) * SIGN2(st2) * pow(fabs(st2),1.0/2));
  p->z = c*SIGN2(st1) * pow(fabs(st1),2.0/p3);
  //printf("z actual=%.15G ellips=%.15G\n", c*SIGN2(st1) * pow(fabs(st1),1.0/p3), c*SIGN2(st1) * pow(fabs(st1),1.0/2.0));
  // evaluate SQ normal vector here 
  en->x = SIGN2(p->x)*p1*pow(fabs(p->x),p1-1.0)/pow(a,p1);
  en->y = SIGN2(p->y)*p2*pow(fabs(p->y),p2-1.0)/pow(b,p2);
  en->z = SIGN2(p->z)*p3*pow(fabs(p->z),p3-1.0)/pow(c,p3); 

  Normalise(en);
  
#if 0
  tmp  = SIGN(ct1) * pow(fabs(ct1),p1);
   p->x = a * tmp * SIGN(ct2) * pow(fabs(ct2),p2);
   p->y = b * SIGN(st1) * pow(fabs(st1),p1);
   p->z = c * tmp * SIGN(st2) * pow(fabs(st2),p2);
#endif
}
void EvalSuperQuadrics(double t1,double t2,double p1,double p2,double p3,
		      double a, double b, double c, XYZ *p)
{
  double ct1,ct2,st1,st2, tmp;

  ct1 = cos(t1);
  ct2 = cos(t2);
  st1 = sin(t1);
  st2 = sin(t2);
  //p1=p2=p3=2.0;

  p->x = a*SIGN2(ct1) * pow(fabs(ct1),2.0/p1) * SIGN2(ct2) * pow(fabs(ct2),2.0/p1);
  //printf("x actual=%.15G ellips=%.15G\n", p->x, a*SIGN2(ct1) * pow(fabs(ct1),1.0/2) * SIGN2(ct2) * pow(fabs(ct2),1.0/2));
  p->y = b*SIGN2(ct1) * pow(fabs(ct1),2.0/p2) * SIGN2(st2) * pow(fabs(st2),2.0/p2);
  //printf("y actual=%.15G ellips=%.15G\n", p->y, b*SIGN2(ct1) * pow(fabs(ct1),1.0/2) * SIGN2(st2) * pow(fabs(st2),1.0/2));
  p->z = c*SIGN2(st1) * pow(fabs(st1),2.0/p3);
  //printf("z actual=%.15G ellips=%.15G\n", c*SIGN2(st1) * pow(fabs(st1),1.0/p3), c*SIGN2(st1) * pow(fabs(st1),1.0/2.0));
#if 0
  tmp  = SIGN(ct1) * pow(fabs(ct1),p1);
  p->x = a * tmp * SIGN(ct2) * pow(fabs(ct2),p2);
  p->y = b * SIGN(st1) * pow(fabs(st1),p1);
  p->z = c * tmp * SIGN(st2) * pow(fabs(st2),p2);
#endif
}
void render_one_spot(double nx, double ny, double nz, double spotradius, 
		     int spotcol, double spotangle, float fadeFact, int red, int green, int blue)
{
  double rax, ray, raz, rotangle, normra, normn;
  double Pi;
  Pi = 2.0*acos(0);
  float col[4];
 
  glPushMatrix();
  vectProd(0,1,0, nx, ny, nz, &rax, &ray, &raz);
  if (rax==0 && ray==0 && raz==0)
    {
      if (ny < 0)
	{
	  rax = 0;
	  ray = 0;
	  raz = 1;
	  rotangle=180;
	}
      else
	{
	  rax = 0;
	  ray = 0;
	  raz = 1;
	  rotangle = 0;
	}
    }
  else
    {
      normra = sqrt(Sqr(rax)+Sqr(ray)+Sqr(raz));
      normn = sqrt(Sqr(nx)+Sqr(ny)+Sqr(nz));
      rotangle = 180.0*acos(ny/normn)/Pi;
      //printf("n=%f %f %f r=%f %f %f\n", nx, ny, nz, rax, ray, raz);
      //printf("rotangle=%f\n", rotangle);
    }
  glRotatef(rotangle, rax, ray, raz);
  glRotatef(180, 1, 0, 0);/* up-down flip */
  /* render the spot here! */
  setColorRGB(spotcol, fadeFact, red, green, blue);
  CreatePartialSuperEllipse(1, 1, spotradius, spotradius, spotradius, globset.stacks, 
			    globset.slides, 1, 0.0, spotangle);
  glPopMatrix(); 
}
#ifdef SQ_REND_SYM 
void swap_axes(double *a, double *b, double *c, double *p1, double *p2, double *p3)
{
  double aL, bL, cL, p1L, p2L, p3L;
  aL = *a;
  bL = *b;
  cL = *c;
  p1L = *p1;
  p2L = *p2;
  p3L = *p3;
  /* scambia i parametri in modo che l'asse di simmetria giaccia sempre lungo l'asse z */
  if ((*p2==*p3) && (*b==*c)) 
    {
      *a = bL;
      *b = cL;
      *c = aL;
      *p1 = p2L;
      *p2 = p3L;
      *p3 = p1L;
    }
  else if ((*p1==*p3) && (*a==*c))
    {
      *a = aL;
      *c = bL;
      *b = cL; 
      *p1 = p1L;
      *p3 = p2L;
      *p2 = p3L;
    }
}
int create_rotate_to_z(atom_s *atom, float rotm[16])
{
  /* it returns 1 if rotation is not indentity */
  double a,b,c,p1,p2,p3;
  int k1, k2;
  /* ruota l'asse di simmetria z sull'asse x o y dopo che 
     swap_axes() ha portato l'asse di simmetria sull'asse z */
  p1 = atom->supquadrics.n1;
  p2 = atom->supquadrics.n2;
  p3 = atom->supquadrics.n3;
  a = atom->supquadrics.a;
  b = atom->supquadrics.b;
  c = atom->supquadrics.c;
  for (k1 = 0; k1 < 4; k1++)
    for (k2 = 0; k2 < 4; k2++)
      {
	if (k1==k2)
	  rotm[k1*4+k2] = 1.0;
	else
	  rotm[k1*4+k2]=0.0;
      }
  if (p2==p3 && b==c)
    {
      rotm[0] = 0.0;
      rotm[10] = 0.0;
      rotm[2] = -1.0;
      rotm[8] = 1.0;
      return 1;
    }
  else if (p1==p3 && a==c)
    {
      rotm[6] = -1.0;
      rotm[5] = 0.0;
      rotm[9] = 1.0;
      rotm[10] = 0.0;
      return 1;
    }
  else
    {
      return 0;
    }
}
#endif
/* ========================== >>> displayMol <<< ===========================*/
void displayAtom(int nf, int nm, int na)
{
  float fadeFact;
  float rotm[16], rotz[16];
  struct spotlst *sl;
  GLUquadricObj *ss, *ss2, *ss3;
  atom_s *atom;
  int k1, k2, do_rotz;
  double rax, ray, raz, rotangle, normra, normn, Pi, redrad, ang, dh;
  int ms;
  glPushMatrix();
  Pi = 2.0*acos(0);
  atom = &mols[nf][nm].atom[na];
  glTranslatef(atom->common.rx,atom->common.ry,atom->common.rz);/* 1st atom */ 
  
  fadeFact = atom->common.transp*calcFadeFact(globset.fadeMode, nf);
  if (atom->common.greyLvl)
    {
      setColor(mgl_bw[atom->common.greyLvl], fadeFact);
    }
  else 
    {
      if (globset.bw)
	{
	  if (globset.NA)
	    setColor(mgl_bw[globset.colIdxBw[na]], fadeFact);
	  else
	    setColor(mgl_bw[globset.default_bw], fadeFact);
	}
      else
	{
	  if (atom->common.atcol>=0 && atom->common.atcol<NUMCOLS)
	    {
	      setColor(mgl_col[atom->common.atcol].rgba, fadeFact);
	    }
	  else if (atom->common.atcol==-2) /* -2 = RGB format */
	    {
	      setColorRGB(atom->common.atcol, fadeFact, atom->common.atred, atom->common.atgreen, atom->common.atblue);
	    }	    
	  else
	    {
	      if (globset.NA)
		setColor(mgl_col[globset.colIdxCol[na]].rgba, fadeFact);
	      else
		setColor(mgl_col[globset.default_col].rgba, fadeFact);
	    }
	}
    }
  /* 
  glEnable (GL_BLEND);
  if (atom->common.transp < 1.0)
    {
      //glDepthMask (GL_FALSE);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
  else
    {
      //glDepthMask (GL_TRUE);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    */
  if (atom->common.type==MGL_ATOM_SPHERE)
    {
      glutSolidSphere (atom->sphere.radius, globset.stacks, globset.slides);
      
    }
  else if (atom->common.type==MGL_ATOM_SPHERE_MSPOT)
    {
      /* draw the sphere first with reduced radius to arrange spot later on it */
      //printf("quiiiii\n");
      ms = (globset.stacks < globset.slides)?globset.stacks:globset.slides;
      ang = TWOPI/ms;
      dh = atom->sphere_mspot.a*(1.0-cos(ang*0.5));
      //printf("ms=%d rad=%.15G dh=%.15G\n", ms, atom->sphere_mspot.a, dh);
      redrad = atom->sphere_mspot.a-dh*2.0;
      glutSolidSphere (redrad, globset.stacks, globset.slides);
      sl = atom->sphere_mspot.sl;
      //printf("sl=%p\n", sl);
      while (sl)
	{
	  //printf("sl=%p col=%d angle=%.15G\n", sl, sl->spotcol, sl->spotangle);
	  render_one_spot(sl->n[0], sl->n[1], sl->n[2], atom->sphere_mspot.a, 
			  sl->spotcol, sl->spotangle, fadeFact, sl->red, sl->green, sl->blue);
	  sl = sl->next;
	}
    }
  else if (atom->common.type==MGL_ATOM_SPHERE_SPOT)
    {
      /* qui si deve orientare il superellissoide */
#if 0
      for (k1 = 0; k1 < 4; k1++)
	for (k2 = 0; k2 < 4; k2++)
	  {
	    if (k1 < 3 && k2 < 3)
	      {
		rotm[k1*4+k2]=atom->sphere_spot.R[k2][k1];
	      }
	    else if (k1==3 && k2 ==3)
	      rotm[15] = 1.0;
	    else
	      rotm[k1*4+k2] = 0.0;
	    //printf("rotm[%d]:%f\n", k1*4+k2, rotm[k1*4+k2]);
	  }
#endif
      /* notare che x' = R x quindi:
       * x = Inversa(R) x' = Trasposta(R) x'*/
      //glMultTransposeMatrixf(rotm);
      vectProd(0,1,0,atom->sphere_spot.nx, atom->sphere_spot.ny, atom->sphere_spot.nz,
		   &rax, &ray, &raz);

      if (rax==0 && ray==0 && raz==0)
	{
	  if (atom->sphere_spot.ny < 0)
	    {
	      rax = 0;
	      ray = 0;
	      raz = 1;
	      rotangle=180;
	    }
	  else
	    {
	      rax = 0;
	      ray = 0;
	      raz = 1;
	      rotangle = 0;
	    }
	}
      else
       	{
	  normra = sqrt(Sqr(rax)+Sqr(ray)+Sqr(raz));
	  normn = sqrt(Sqr(atom->sphere_spot.nx)+Sqr(atom->sphere_spot.ny)+Sqr(atom->sphere_spot.nz));
	  rotangle = 180.0*acos(atom->sphere_spot.ny/normn)/Pi;
	}
      glRotatef(rotangle, rax, ray, raz);
      glRotatef(180, 1, 0, 0);/* up-down flip */
      //printf("rotangle=%.15G rax=%f %f %f\n", rotangle, rax, ray, raz);
      CreatePartialSuperEllipse(atom->sphere_spot.n1, 
				atom->sphere_spot.n2, atom->sphere_spot.a, 
				atom->sphere_spot.b, atom->sphere_spot.c, globset.stacks, 
				globset.slides, //1, 0.0, atom->sphere_spot.tbeg);
				1, atom->sphere_spot.tbeg, TWOPI/2.0);
      
      setColorRGB(atom->sphere_spot.spotcol, fadeFact, atom->sphere_spot.red, atom->sphere_spot.green,
		  atom->sphere_spot.blue);
      CreatePartialSuperEllipse(atom->sphere_spot.n1, 
				atom->sphere_spot.n2, atom->sphere_spot.a, 
				atom->sphere_spot.b, atom->sphere_spot.c, globset.stacks, 
				globset.slides, //1, atom->sphere_spot.tbeg, TWOPI/2.0);
				1, 0.0, atom->sphere_spot.tbeg);
    }
  else if (atom->common.type==MGL_ATOM_DISK)
    {
#if 0
      printf("qui radius=%f height=%f\n",  atom->disk.radius,atom->disk.height );
#endif
      ss = gluNewQuadric();
      ss2 =  gluNewQuadric();
      ss3 =  gluNewQuadric();
      vectProd(0,0,1,atom->disk.nx, atom->disk.ny, atom->disk.nz,
	       &rax, &ray, &raz);
      if (rax==0 && ray==0 && raz==0)
	{
	  if (atom->sphere_spot.nz < 0)
	    {
	      rax = 0;
	      ray = 1;
	      raz = 0;
	      rotangle = 180;
	    }
	  else
	    {
	      rax = 0;
	      ray = 1;
	      raz = 0;
	      rotangle = 0;
	    }
	}
      else
	{
	  normra = sqrt(Sqr(rax)+Sqr(ray)+Sqr(raz));
	  normn = sqrt(Sqr(atom->disk.nx)+Sqr(atom->disk.ny)+Sqr(atom->disk.nz));
	  rotangle = 180.0*acos(atom->disk.nz/normn)/Pi; 	       
	}
#if 0
      printf("Rotation Angle: %f around (%f,%f,%f) n(%f,%f,%f)\n", 
	     rotangle, rax, ray, raz,atom->disk.nx, atom->disk.ny,acos(atom->disk.nz/normn)/Pi );
#endif
      glRotatef(rotangle, rax, ray, raz);
      glPushMatrix();
      glTranslatef(0, 0, atom->disk.height*0.5);
      gluDisk(ss2, 0, atom->disk.radius, globset.stacks, globset.slides);
      glPopMatrix();
      /*
      if (atom->common.transp < 1.0)
	glDepthMask (GL_FALSE);*/
      glTranslatef(0, 0, -atom->disk.height*0.5);
      gluCylinder(ss, atom->disk.radius, 
		      atom->disk.radius, 
		      atom->disk.height, 
		      globset.stacks, globset.slides);
      /* if (atom->common.transp < 1.0)
	glDepthMask (GL_TRUE); */
      glPushMatrix();
      gluQuadricOrientation(ss3, GLU_INSIDE);
      gluDisk(ss3, 0, atom->disk.radius, globset.stacks, globset.slides);
      glPopMatrix();
    }
  else if (atom->common.type==MGL_ATOM_CYLINDER)
    {
      /* orienta il cilindro */
      ss = gluNewQuadric();
      vectProd(0,0,1,atom->disk.nx, atom->disk.ny, atom->disk.nz,
	       &rax, &ray, &raz);
      normra = sqrt(Sqr(rax)+Sqr(ray)+Sqr(raz));
      normn = sqrt(Sqr(atom->disk.nx)+Sqr(atom->disk.ny)+Sqr(atom->disk.nz));
      rotangle = 180.0*acos(atom->disk.nz/normn)/Pi; 	
      glRotatef(rotangle, rax, ray, raz);
      glTranslatef(0, 0, -atom->disk.height*0.5);
      /*
	 if (atom->common.transp < 1.0)
	 glDepthMask (GL_FALSE);*/
      gluCylinder(ss, atom->cylinder.toprad, 
		      atom->cylinder.botrad, 
		      atom->cylinder.height, 
		      globset.stacks, globset.slides);
      /*if (atom->common.transp < 1.0)
	glDepthMask (GL_TRUE);*/
    }
  else if (atom->common.type==MGL_ATOM_SUPELLIPS)
    {
      /* qui si deve orientare il superellissoide */
      for (k1 = 0; k1 < 4; k1++)
	for (k2 = 0; k2 < 4; k2++)
	  {
	    if (k1 < 3 && k2 < 3)
	      {
		rotm[k1*4+k2]=atom->supellips.R[k2][k1];
	      }
	    else if (k1==3 && k2 ==3)
	      rotm[15] = 1.0;
	    else
	      rotm[k1*4+k2] = 0.0;
	    //printf("rotm[%d]:%f\n", k1*4+k2, rotm[k1*4+k2]);
	  }
      /* notare che x' = R x quindi:
       * x = Inversa(R) x' = Trasposta(R) x'*/
      glMultTransposeMatrixf(rotm);
      //glEnable (GL_BLEND);
      /*if (atom->common.transp < 1.0)
	glDepthMask (GL_FALSE);*/
      //glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      /* for now disabled */
      if (atom->supellips.n1==1 && atom->supellips.n2==1)
	{
	  //printf("qui st=%d sl=%d\n", globset.stacks, globset.slides);
	  for (k1 = 0; k1 < 4; k1++)
	    for (k2 = 0; k2 < 4; k2++)
	      {
		if (k1 < 3 && k2 < 3 && k1==k2)
		  {
		    switch(k1)
		      {
		      case 0:
			rotm[k1*4+k2]=atom->supellips.a;
		      break;
		      case 1:
			rotm[k1*4+k2]=atom->supellips.b;
		      break;
		      case 2:
			rotm[k1*4+k2]=atom->supellips.c;
		      break;
		      }
		  }
		else if (k1==3 && k2 ==3)
		  rotm[15] = 1.0;
	    else
	      rotm[k1*4+k2] = 0.0;
	    //printf("rotm[%d]:%f\n", k1*4+k2, rotm[k1*4+k2]);
	  }
	  glEnable(GL_NORMALIZE);
	  glMultMatrixf(rotm);
       	  glutSolidSphere (1, globset.stacks, globset.slides);
	  glDisable(GL_NORMALIZE);
	}	  
      else
	{
	  if (atom->supellips.tbeg > 0.0 || atom->supellips.tend > 0.0)
	    {
	      CreatePartialSuperEllipse(atom->supellips.n1, 
	  				atom->supellips.n2, atom->supellips.a, 
					atom->supellips.b, atom->supellips.c, globset.stacks, 
					globset.slides, 1, atom->supellips.tbeg, atom->supellips.tend);
	    }
	  else
	    CreateSuperEllipse(atom->supellips.n1, 
		  	       atom->supellips.n2, atom->supellips.a, 
		  	       atom->supellips.b, atom->supellips.c, globset.stacks, 
		  	       globset.slides, 1);
	}
      /*if (atom->common.transp < 1.0)
	glDepthMask (GL_TRUE);*/
      //glDisable (GL_BLEND); 
    }
  else if (atom->common.type==MGL_ATOM_SUPQUADRICS)
    {

      /* qui si deve orientare il superellissoide */
      for (k1 = 0; k1 < 4; k1++)
	for (k2 = 0; k2 < 4; k2++)
	  {
	    if (k1 < 3 && k2 < 3)
	      {
		rotm[k1*4+k2]=atom->supquadrics.R[k2][k1];
	      }
	    else if (k1==3 && k2 ==3)
	      rotm[15] = 1.0;
	    else
	      rotm[k1*4+k2] = 0.0;
	    //printf("rotm[%d]:%f\n", k1*4+k2, rotm[k1*4+k2]);
	  }
      /* notare che x' = R x quindi:
       * x = Inversa(R) x' = Trasposta(R) x'*/
      glMultTransposeMatrixf(rotm);
      /* porta l'asse di simmetria a coincidere con l'asse z*/
      //glEnable (GL_BLEND);
      /*if (atom->common.transp < 1.0)
	glDepthMask (GL_FALSE);*/
      //glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      /* for now disabled */
#ifdef SQ_REND_SYM
      do_rotz = create_rotate_to_z(atom, rotz);
      if (do_rotz)
	glMultMatrixf(rotz);
#endif
      if (atom->supquadrics.n1==2 && atom->supquadrics.n2==2 &&
     	  atom->supquadrics.n3==2)
	{
	  //printf("qui st=%d sl=%d\n", globset.stacks, globset.slides);
	  for (k1 = 0; k1 < 4; k1++)
	    for (k2 = 0; k2 < 4; k2++)
	      {
		if (k1 < 3 && k2 < 3 && k1==k2)
		  {
		    switch(k1)
		      {
		      case 0:
			rotm[k1*4+k2]=atom->supquadrics.a;
		      break;
		      case 1:
			rotm[k1*4+k2]=atom->supquadrics.b;
		      break;
		      case 2:
			rotm[k1*4+k2]=atom->supquadrics.c;
		      break;
		      }
		  }
		else if (k1==3 && k2 ==3)
		  rotm[15] = 1.0;
	    else
	      rotm[k1*4+k2] = 0.0;
	    //printf("rotm[%d]:%f\n", k1*4+k2, rotm[k1*4+k2]);
	  }
	  glEnable(GL_NORMALIZE);
	  glMultMatrixf(rotm);
       	  glutSolidSphere (1, globset.stacks, globset.slides);
	  glDisable(GL_NORMALIZE);
	}	  
      else
	{
	  if (atom->supquadrics.tbeg > 0.0 || atom->supquadrics.tend > 0.0)
	    {
	      CreatePartialSuperQuadrics(atom->supquadrics.n1, 
	  				atom->supquadrics.n2, 
					atom->supquadrics.n3,
				       	atom->supquadrics.a, 
					atom->supquadrics.b, atom->supquadrics.c, globset.stacks, 
					globset.slides, 1, atom->supquadrics.tbeg, 
					atom->supquadrics.tend);
	    }
	  else
	    {
	      //printf("boh n1=%G n2=%G n3=%G\n", atom->supquadrics.n1, atom->supquadrics.n2,
//atom->supquadrics.n3);
#if 0 
	      CreateSuperEllipse(atom->supquadrics.n1, 
	   			  atom->supquadrics.n2,
	   			  atom->supquadrics.a, 
	   			  atom->supquadrics.b, atom->supquadrics.c, globset.stacks, 
	   			  globset.slides, 1);
#else
	      CreateSuperQuadrics(atom->supquadrics.n1, 
	   			  atom->supquadrics.n2, atom->supquadrics.n3,
	   			  atom->supquadrics.a, 
	   			  atom->supquadrics.b, atom->supquadrics.c, globset.stacks, 
	   			  globset.slides, 1);
#endif
	      //printf("qui\n");
	    }
	}
      /*if (atom->common.transp < 1.0)
	glDepthMask (GL_TRUE);*/
      //glDisable (GL_BLEND); 
    }
  /*if (atom->common.transp < 1.0)
    glDepthMask (GL_TRUE); 
    glDisable (GL_BLEND); */ 
/*
    gluSphere ss = gluNewQuadric();(ss, sig[na], 12, 12);
  */
  glPopMatrix();

}
void renderSolidCylinder(double rx, double ry, double rz, 
			 double radius, double height)
{
  GLUquadricObj *ss, *ss2, *ss3;
  glPushMatrix();
  glTranslatef(0, 0, height);
  ss = gluNewQuadric();
  ss2 =  gluNewQuadric();
  ss3 =  gluNewQuadric();
  gluDisk(ss2, 0, radius, STACKS, SLIDES);
  glPopMatrix();
  gluCylinder(ss, radius, 
	      radius, 
	      height, 
	      STACKS, SLIDES);
  glPushMatrix();
  gluQuadricOrientation(ss3, GLU_INSIDE);
  gluDisk(ss3, 0, radius, STACKS, SLIDES);
  glPopMatrix();
}
void displayBonds(int nf, int i)
{
  int nb, from, to;
  double rcmx, rcmy, rcmz, rax, ray, raz, normra, normn;
  double nx, ny, nz, Pi, rotangle;
  float fadeFact;
  
  fadeFact = calcFadeFact(globset.fadeMode, nf);
  
  Pi = 2.0 * acos(0);
  for (nb = 0; nb < mols[nf][i].numbonds; nb++)
    {
      glPushMatrix();
      from = mols[nf][i].bond[nb].from;
      to   = mols[nf][i].bond[nb].to;
#if 0
      printf("nb: %d from:%d to:%d\n", nb, from, to);
#endif
      if (from < 0 || from >= mols[nf][i].numat ||
	  to < 0 || to >= mols[nf][i].numat)
	{
	  fprintf(stderr,
		  "WARNING: Bond %d-%d of Molecule N. %i in frame %d refers to non-existant atom\n", 
		  from, to, i, nf);
	  continue;
	}
      rcmx = mols[nf][i].atom[from].common.rx;
      rcmy = mols[nf][i].atom[from].common.ry;
      rcmz = mols[nf][i].atom[from].common.rz;
      nx = mols[nf][i].atom[to].common.rx - mols[nf][i].atom[from].common.rx;
      ny = mols[nf][i].atom[to].common.ry - mols[nf][i].atom[from].common.ry;
      nz = mols[nf][i].atom[to].common.rz - mols[nf][i].atom[from].common.rz;
#if 0
      printf("rx[%d]:%f rx[%d]:%f\n", from,  mols[nf][i].atom[from].common.rx, to,
	     mols[nf][i].atom[to].common.rx);
#endif
      glTranslatef(rcmx, rcmy, rcmz);
      vectProd(0,0,1, nx, ny, nz, &rax, &ray, &raz);
      normra = sqrt(Sqr(rax)+Sqr(ray)+Sqr(raz));
      normn = sqrt(Sqr(nx)+Sqr(ny)+Sqr(nz));
      //printf("normn:%f\n", normn);
      if (normn==0)
	continue;
      rotangle = 180.0*acos(nz/normn)/Pi; 	       
#if 0
      printf("Rotation Angle: %f around (%f,%f,%f) n(%f,%f,%f)\n", 
	     rotangle, rax, ray, raz,atom->disk.nx, atom->disk.ny,acos(atom->disk.nz/normn)/Pi );
#endif
      glRotatef(rotangle, rax, ray, raz);
      if (mols[nf][i].bond[nb].color>=0 && mols[nf][i].bond[nb].color<NUMCOLS)
	{
	  setColor(mgl_col[mols[nf][i].bond[nb].color].rgba, fadeFact);
	}
      else if (mols[nf][i].bond[nb].color==-2)
    	setColorRGB(mols[nf][i].bond[nb].color, fadeFact, mols[nf][i].bond[nb].red,
		    mols[nf][i].bond[nb].green, mols[nf][i].bond[nb].blue);
      else
	{
	  setColor(mgl_col[globset.default_col].rgba, fadeFact);
	}

      /*if (mols[nf][i].bond[nb].transp < 1.0)
	glDepthMask (GL_FALSE);*/
      renderSolidCylinder(rcmx, rcmy, rcmz, mols[nf][i].bond[nb].thickness, normn);
      /*if (mols[nf][i].bond[nb].transp < 1.0)
	glDepthMask (GL_TRUE);*/
      glPopMatrix();
    }
}
/* ========================= >>> buildAtomsList <<< ======================= */
void buildAtomsList()
{
  int i, j, nf;
  /*float col[4];*/

  /* NOTE: frameno is the number of frames read from postions file */
  for (nf = 0; nf < globset.frameNo; ++nf )
    {
      atomsList[nf] = glGenLists(1);
      glNewList(atomsList[nf], GL_COMPILE);
      /* printf("numols[%d]:%d\n", nf, globset.NumMols[nf]);*/
      for(i = 0; i < globset.NumMols[nf]; ++i)
	{
	  /*printf("mols[%d][%d].numat:%d\n", nf, i, mols[nf][i].numat);*/
	  for (j=0; j < mols[nf][i].numat; j++)
	    displayAtom(nf, i, j);
	  displayBonds(nf, i);
	}
      glEndList();
    }
}

void myReshape(int w, int h);

/* ======================== >>> printStr <<< =============================== */
void printStr(int x, int y, const char* text)
{
  /* print a string at the position (x,y) */
  int i;
  
  glRasterPos2i(x, y);
  for (i = 0; i < strlen(text); ++i)
    {
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
    }
}
/* ======================== >>> onScreenInfo <<< ========================== */
void onScreenInfo()
{
  char text[255];
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, globset.Width, 0, globset.Height);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glColor3f(0.8, 0.0, 0.0);
  sprintf(text, "degx: %.1f degy: %.1f degz: %.1f", globset.degx, globset.degy, globset.degz);
  printStr(10, globset.Height - 20, text);
  sprintf(text, "L: %.1f deginc: %.1f", globset.L, globset.deginc);
  printStr(10, globset.Height - 45, text);
  myReshape(globset.Width, globset.Height);
  glPopAttrib();
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
}
#ifdef MGL_MACOS
#include<png.h>
#else
#include<png.h>
#endif
/*extern PNG_EXPORT(void,png_set_zbuf_size);*/
void save_image(void)
{
   FILE *fp;
   char *fn;
   char numstr[300];
   png_structp png_ptr;
   png_infop info_ptr;
   png_bytep *row_pointer;
   unsigned char *image;
   int width,height;
   unsigned int NbBytes;
   int i, j ;
   int pixel_size=3;
   /* Allocate our buffer for the image */
   width= glutGet(GLUT_WINDOW_WIDTH);
   height=glutGet(GLUT_WINDOW_HEIGHT) ;
   NbBytes= 3*width*height*sizeof(char);
   if ((image = malloc(NbBytes)) == NULL) {
      fprintf(stderr,"Failed to allocate memory for image\n");
      return;
   }

   glPixelStorei(GL_PACK_ALIGNMENT,1);

   /* Open the file */
   fn = malloc(sizeof(char)*(strlen(globset.savefile)+257));
   strcpy(fn,globset.savefile);
   if (globset.saved_counter > 0)
     {
       snprintf(numstr, 256, "%d",globset.saved_counter);
       strcat(fn, numstr);
     }
   if ((fp = fopen(fn,"wb")) == NULL) {
      fprintf(stderr,"Failed to open file for window dump\n");
      free(fn);
      free(image);
      return;
   }

   png_ptr = png_create_write_struct
       (PNG_LIBPNG_VER_STRING, NULL,
	NULL, NULL);
   if (!png_ptr)
     {
       free(image);
       free(fn);
       return;
     }
   info_ptr = png_create_info_struct(png_ptr);
   if (!info_ptr)
     {
       png_destroy_write_struct(&png_ptr,
			 	(png_infopp)NULL);
       free(image);
       free(fn);
       return;
     }
   if (setjmp(png_jmpbuf(png_ptr)))
    {
       png_destroy_write_struct(&png_ptr, &info_ptr);
       free(image);
       free(fn);
       fclose(fp);
       return;
    }
   png_init_io(png_ptr, fp);
    /* turn on or off filtering, and/or choose
       specific filters.  You can use either a single
       PNG_FILTER_VALUE_NAME or the logical OR of one
       or more PNG_FILTER_NAME masks. */
   /*
   png_set_filter(png_ptr, 0,
       PNG_FILTER_NONE  | PNG_FILTER_VALUE_NONE |
       PNG_FILTER_SUB   | PNG_FILTER_VALUE_SUB  |
       PNG_FILTER_UP    | PNG_FILTER_VALUE_UP   |
       PNG_FILTER_AVE   | PNG_FILTER_VALUE_AVE  |
       PNG_FILTER_PAETH | PNG_FILTER_VALUE_PAETH|
       PNG_ALL_FILTERS);*/
   /* set the zlib compression level */
   png_set_compression_level(png_ptr,
		      	     Z_BEST_COMPRESSION);
   
   /* set other zlib parameters */
   png_set_compression_mem_level(png_ptr, 8);
   png_set_compression_strategy(png_ptr,
				Z_DEFAULT_STRATEGY);
   png_set_compression_window_bits(png_ptr, 15);
   png_set_compression_method(png_ptr, 8);
   png_set_compression_buffer_size(png_ptr, 8192);
   png_set_IHDR(png_ptr, info_ptr, width, height,
       8,  PNG_COLOR_TYPE_RGB,  PNG_INTERLACE_NONE,
       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT); 

   /* Copy the image into our buffer */
   //glReadBuffer(GL_BACK);
   glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);
   row_pointer=(png_bytep*)malloc(sizeof(png_bytep)*height);
   row_pointer = png_malloc(png_ptr,
      height*sizeof(png_bytep));
   for (i=0; i<height; i++)
      row_pointer[i]=png_malloc(png_ptr,
         width*pixel_size);
   /* Write the raw file */
   /* fprintf(fptr,"P6\n%d %d\n255\n",width,height); for ppm */
   for (j=height-1;j>=0;j--) {
      for (i=0;i<width;i++) {
         row_pointer[height-j-1][3*i+0]=image[3*j*width+3*i+0];
         row_pointer[height-j-1][3*i+1]=image[3*j*width+3*i+1];
         row_pointer[height-j-1][3*i+2]=image[3*j*width+3*i+2];
      }
   }
   png_set_rows(png_ptr, info_ptr, row_pointer);
   png_write_png(png_ptr, info_ptr,  PNG_TRANSFORM_IDENTITY, NULL);
   
   fclose(fp);
   /* Clean up */
   free(fn);
   free(image);
}
void display (void);
void timerCB(int val)
{
  globset.nrefresh=1;
  display();
}
/* ======================== >>> display <<< ===============================*/
void display (void)
{
  static int count=0;
  int nf;
  int i, j;
  /*
   * double fadeFact;
   */
#if 0
  if (globset.saveandquit)
    {
      glutHideWindow();
    }
#endif
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (globset.setvp)
    {
      glLoadIdentity();
      gluLookAt(globset.ivpx,globset.ivpy,globset.ivpz,0, 0, 0, 0, 1, 0);
      glPushMatrix ();
    }
  else
    {
      glLoadIdentity();
      
      if (!globset.axon) 
	glTranslatef(0.0, 0.0, -(globset.near+globset.L/2.0));
      else
	glTranslatef(0.0, 0.0,-globset.L/2.0);

      glTranslatef(0.0, 0.0, -globset.dist);
      glPushMatrix ();
      /* NOTE:
	 First arg is the degrees of the rotation, the others are the component
	 of the vector around which we perfomr the rotation */
      glRotatef(-globset.degz, 0.0, 0.0, 1.0);
      glRotatef(-globset.degy, 0.0, 1.0, 0.0);
      glRotatef(-globset.degx, 1.0, 0.0, 0.0);
    }
  setColor(mgl_bw[0], 1.0);
  if (globset.drawcube)
    glutWireCube(globset.L);
 
  for (nf = 0; nf < globset.frameNo; ++nf)
    {
#if 1
      if (nf == 0)
	{
	  if (globset.depthmask)
	    glDepthMask(GL_TRUE);
	  else
	    glDepthMask(GL_FALSE);
	}
      else 
	{
	  glDepthMask(GL_FALSE);
	}
#endif
#ifdef MGL_USELIST
      if (globset.mgl_uselist)
	{
	  glCallList(atomsList[nf]);
	}
      else
	{
	  for(i = 0; i < globset.NumMols[nf]; ++i)
	    {
	      for (j=0; j < mols[nf][i].numat; j++)
		displayAtom(nf, i, j);
	      displayBonds(nf, i);
	    }
	}
#else
      for(i = 0; i < globset.NumMols[nf]; ++i)
	{
	  for (j=0; j < mols[nf][i].numat; j++)
	    displayAtom(nf, i, j);
	  displayBonds(nf, i);
	}
#endif
    }
  glPopMatrix ();
  if (!globset.depthmask)
    glDepthMask(GL_TRUE); 
  if (globset.infos) onScreenInfo();
  if (globset.saveandquit==1)
    count++;
  if (globset.saveandquit==1)
    glFinish(); 
  glutSwapBuffers();
  if (globset.nrefresh > 1 && count==1 && globset.exitDelay==0)
    globset.exitDelay=1000; /* dopo 1 sec in ogni caso esce */
  //printf("count=%d exitDelay=%d nrefresh=%d QUI\n", count, globset.exitDelay, globset.nrefresh);
  if (globset.exitDelay > 0 && count==1)
    {
      //glutPostRedisplay();
      glutTimerFunc(globset.exitDelay, timerCB, 1);
      return;
    }
  if (globset.saveandquit==1)
    {
      if (globset.nrefresh==1 || (globset.nrefresh > 1 && count >= globset.nrefresh))
	{
	  //printf("nrefresh=%d count=%d exitDelay=%d\n", globset.nrefresh, count, globset.exitDelay);
	  save_image();
	  exit(0);
	}
    }
}


/* =========================== >>> setproj <<< ============================*/
void setproj(void)
{
  double L;
  L = globset.L;
  /*int w = Width, h = Height;*/
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  if (globset.axon)
    {
      /*printf("orto\n");*/
      glOrtho (-L, L, -L, 
      	       L, -2.0*L, 2.0*L);
      /*if (w <= h) 
	glOrtho (-L/2.0, L/2.0, -L*(GLfloat)h/(GLfloat)w/2.0, 
		 L*(GLfloat)h/(GLfloat)w/2.0, -L/2.0, L/2.0);
      else 
	glOrtho (-L*(GLfloat)w/(GLfloat)h/2.0, 
	L*(GLfloat)w/(GLfloat)h/2.0, -L/2.0, L/2.0, -L/2.0, L/2.0);*/
    }
  else 
    {
      globset.near = L / tan(PI*globset.viewangle/360.0) / 2.0;
      globset.far = 100*globset.near;
      gluPerspective(globset.viewangle, (GLdouble) globset.Width / (GLdouble) globset.Height, 
		     globset.near, globset.far);
    }

  glMatrixMode (GL_MODELVIEW);
}

/* ========================== >>> myReshape <<< ========================*/
void myReshape(int w, int h)
{
  globset.Width = w;
  globset.Height = h;
  glViewport (0, 0, w, h);
  setproj();
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity();
}
void print_usage(void)
{
  printf ("USAGE:\n");
  printf("molgl [-h/--help | --uselist/-ul | --saveandquit/-sq | --pngfile/-f <filename> \n");
  printf("| --viewpoint/-vp (x,y,z) | --diameter/-d <atoms_diameter> | --noinfos/-ni\n");
  printf("| --nobox|-nb | --semiax/-sa (a,b,c) | --stacks/-st <stacks>\n");
  printf("| --slides/-sl <slides> | --degreesx/-dgx <rotation_angle>\n");
  printf("| --degreesy/-dgy <rotation_angle> | --degreesz/dgz <rotatioan_angle>\n");
  printf("| --slides/-sl <slides> | --bondtransp/-br <transparency> | --transp/-r <transparency> \n");
  printf("| --boxsize/-L <box_size> | --distance/-di <distance_offset> \n");
  printf("| --twolights/-tl | --light_pos0/-lp0 (x,y,z) | --light_pos1/-lp1 (x,y,z)\n"); 
  printf("| --exitDelay/-ed | --nrefresh/-nr | --height <scr_height> | --width <scr_width> ] <input_file> \n"); 
}
/* ============================= >>> args <<< ============================= */
void args(int argc, char* argv[])
{
  int i=1;
  //printf("argc=%d %s,%s\n",argc, argv[0], argv[1]);
  if (argc == 1)
    {
      print_usage();
      exit(-1);
    }
  else if (argc == 2 && (!strcmp(argv[i], "-h") || !strcmp(argv[i],"--help")))
    {
      print_usage();
      exit(-1);
    }
  else if (argc > 2)
    {
      while (i < argc-1)
	{
	  if (!strcmp(argv[i], "-h") || !strcmp(argv[i],"--help"))
	    {
	      print_usage();
	      exit(-1);
	    }  
#ifdef MGL_USELIST
	  else if (!strcmp(argv[i],"--uselist") || !strcmp(argv[i],"-ul"))
	    {
	      globset.mgl_uselist = 1;
	    }
#endif
	  else if (!strcmp(argv[i],"--saveandquit") || !strcmp(argv[i],"-sq"))
	    {
	      globset.saveandquit = 1;
	    }
	  else if (!strcmp(argv[i],"--twolights") || !strcmp(argv[i],"-tl"))
	    {
	      globset.twolights = 1;
	    }
	  else if (!strcmp(argv[i],"--nodepthmask") || !strcmp(argv[i],"-ndm"))
	    {
	      globset.depthmask = 0;
	    }
	  else if (!strcmp(argv[i],"--pngfile") || !strcmp(argv[i],"-f"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply a file name!\n");
		  exit(-1);
		}
	      globset.savefile = malloc(sizeof(char)*(strlen(argv[i])+1));
	      strcpy(globset.savefile,argv[i]);
	    }
	  else if (!strcmp(argv[i],"--viewpoint")||!strcmp(argv[i],"-vp"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the viewpoint (x,y,z)!\n");
		  exit(-1);
		}
	      sscanf(argv[i], "(%lf,%lf,%lf)", &globset.ivpx, &globset.ivpy, &globset.ivpz);
	      globset.setvp = 1;
	    }
	  else if (!strcmp(argv[i],"--lightpos0")||!strcmp(argv[i],"-lp0"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the light N.0 position (x,y,z)!\n");
		  exit(-1);
		}
	      sscanf(argv[i], "(%f,%f,%f)", &globset.light_pos0[0], &globset.light_pos0[1], &globset.light_pos0[2]);
	    }
	  else if (!strcmp(argv[i],"--lightpos1")||!strcmp(argv[i],"-lp1"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the light N.1 position (x,y,z)!\n");
		  exit(-1);
		}
	      sscanf(argv[i], "(%f,%f,%f)", &globset.light_pos1[0], &globset.light_pos1[1], &globset.light_pos1[2]);
	    }
	  else if (!strcmp(argv[i],"--nrefresh")||!strcmp(argv[i],"-nr"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the number of refresh done before quitting!\n");
		  exit(-1);
		}
	      globset.nrefresh = atoi(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--exitDelay")||!strcmp(argv[i],"-ed"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the delay in msec before quitting!\n");
		  exit(-1);
		}
	      globset.exitDelay = atoi(argv[i]);
	    }
	   else if (!strcmp(argv[i],"--stacks")|| !strcmp(argv[i],"-st"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the number of stacks!\n");
		  exit(-1);
		}
	      globset.stacks = atoi(argv[i]);
	      if (globset.stacks > 1000 || globset.stacks < 1)
		{
		  printf("stacks is better between 1 and 100!\n");
		  exit(-1);
		}
	    }
	  else if (!strcmp(argv[i],"--slides")|| !strcmp(argv[i],"-sl"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the number of slides!\n");
		  exit(-1);
		}
	      globset.slides = atoi(argv[i]);
	      if (globset.slides > 1000 || globset.slides < 1)
		{
		  printf("stacks is better between 1 and 100!\n");
		  exit(-1);
		}
	    }
	  else if (!strcmp(argv[i],"--diameter")|| !strcmp(argv[i],"-d"))
	    {
	      i++;
	      globset.setdiameter = 1;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the diameter!\n");
		  exit(-1);
		}
	      globset.diameter = atof(argv[i]);
	      globset.NA = 1;
	    }
	  else if (!strcmp(argv[i],"--width"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the screen width!\n");
		  exit(-1);
		}
	      globset.Width = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--height"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the screen height!\n");
		  exit(-1);
		}
	      globset.Height = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--semiax")|| !strcmp(argv[i],"-sa"))
	    {
      	      i++;
	      globset.setsemiax = 1;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply three semi-axes (a,b,c)!\n");
		  exit(-1);
		}
	      sscanf(argv[i], "(%lf,%lf,%lf)", &globset.sa, &globset.sb, &globset.sc);
	      globset.NA = 1;
	    }
	  else if (!strcmp(argv[i],"--bondthickness")|| !strcmp(argv[i],"-bt"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the bond thickness!\n");
		  exit(-1);
		}
	      globset.defbondthick = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--bondcolor")|| !strcmp(argv[i],"-bc"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the bond color!\n");
		  exit(-1);
		}
	      globset.defbondcol = atoi(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--distance")|| !strcmp(argv[i],"-di"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the distance!\n");
		  exit(-1);
		}
	      globset.dist = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--boxsize")|| !strcmp(argv[i],"-L"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the box size!\n");
		  exit(-1);
		}
	      globset.L = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--degreesx") || !strcmp(argv[i],"-dgx"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the rotation angle around x-axis!\n");
		  exit(-1);
		}
	      globset.degx = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--degreesy")|| !strcmp(argv[i],"-dgy"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the rotation angle around y-axis!\n");
		  exit(-1);
		}
	      globset.degy = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--degreesz") || !strcmp(argv[i],"-dgz"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the rotation angle around z-axis!\n");
		  exit(-1);
		}
	      globset.degz = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--bondtransp")|| !strcmp(argv[i],"-br"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the default bond transparency!\n");
		  exit(-1);
		}
	      globset.defbondtransp = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--transp")|| !strcmp(argv[i],"-r"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the default transparency!\n");
		  exit(-1);
		}
	      globset.deftransp = atof(argv[i]);
	    }
	 else if (!strcmp(argv[i],"--height")|| !strcmp(argv[i],"-ht"))
	    {
	      i++;
	      globset.setheight = 1;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the viewpoint (x,y,z)!\n");
		  exit(-1);
		}
	      globset.height = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--nobox")|| !strcmp(argv[i],"-nb"))
	    {
	      globset.drawcube = 0;
	    }
	  else if (!strcmp(argv[i],"--noinfos")|| !strcmp(argv[i],"-ni"))
	    {
	      globset.infos = 0;
	    }
	  else
	    {
	      fprintf(stderr, "ERROR: Invalid argumet!\n");
	      exit(-1);
	    }
      	  i++;
	}
    }

  if (i == argc)
    {
      fprintf(stderr, "ERROR: You must supply an input file!\n");
      exit(-1);
    }
  
      
  if (globset.savefile == NULL)
    globset.savefile = "molglimg.png";

  strcpy(inputFile, argv[i]);
}
void dropSpaces(char *S);
int getColByName(const char* name);

int parsecol(char *str, double *transp, float *red, float *green, float *blue)
{
  int colNum;
  char* eptr;
  char cols[128];
  *red = *green = *blue = -1.0;
  /* guess if there is transparency */
  //printf("str:%s\n", str);
  if (sscanf(str, "%[^/]/%lf",cols,transp)==2)
    str=cols;
  else
    *transp = globset.deftransp;
  //printf("deftransp:%.15G\n", globset.deftransp);
  if (sscanf(str,"%f,%f,%f", red, green, blue)==3)
    {
      //printf("red=%f green=%f blue=%f\n", *red, *green, *blue);
      /* RGB format */
      return -2; /* -2 = rgb */
    }
  colNum = (int) strtod(str, &eptr);
  if (eptr == str) /* not a number */
    {
      dropSpaces(str);
      return getColByName(str);
      /* Find the number of the color named 's1'*/
      /*printf("col:%s:, %d\n", S, colIdxCol[j]);
      */
    }
  else
    {
      /*printf("Color Number: %d\n", colNum);*/
      return colNum;
    }
}
/* PARSE parse  parsing  PARSING */
/* ========================== >>> assignAtom <<< ===========================*/
char shl[4096];
void assignAtom(int nf, int i, int a, const char* L)
{
  char s1[128], s2[128], s3[128], s4[128], s5[128], s6[128], s7[128], s8[128], s9[128];
  char s10[128], s11[128], s12[128], s13[128], s14[128], s15[128], s16[128], s17[128], s18[128];
  char s19[128];
  char *ss;
  atom_s *at;
  struct spotlst** sl, *newsp;
  double t;
  at = &mols[nf][i].atom[a];
  //printf("read: %s\n", L);
  if (sscanf(L,"%s %s %s %s %s %s %s %s %s %s %s %s @ %s %s %s C[%[^]]] Q %s %s %s ", s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19) == 19)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
  	printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	atoi(s5),i, j);*/
      //printf("qui\n");
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SUPQUADRICS;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->supquadrics.R[0][0] = atof(s4);
      at->supquadrics.R[0][1] = atof(s5);
      at->supquadrics.R[0][2] = atof(s6);
      at->supquadrics.R[1][0] = atof(s7);
      at->supquadrics.R[1][1] = atof(s8);
      at->supquadrics.R[1][2] = atof(s9);
      at->supquadrics.R[2][0] = atof(s10);
      at->supquadrics.R[2][1] = atof(s11);
      at->supquadrics.R[2][2] = atof(s12);
      at->supquadrics.a = atof(s13);
      at->supquadrics.b = atof(s14);
      at->supquadrics.c = atof(s15);
      at->supquadrics.n1 = atof(s17);
      at->supquadrics.n2 = atof(s18);
      at->supquadrics.n3 = atof(s19);
      at->supquadrics.tbeg = -1.0;
      at->supquadrics.tend = -1.0;
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = parsecol(s16,&t, &(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;

    }
  else if (sscanf(L,"%s %s %s %s %s %s %s %s %s %s %s %s @ %s %s %s C[%[^]]] P %s %s ", s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18) == 18)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
  	printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	atoi(s5),i, j);*/
      //printf("qui\n");
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SUPELLIPS;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->supellips.R[0][0] = atof(s4);
      at->supellips.R[0][1] = atof(s5);
      at->supellips.R[0][2] = atof(s6);
      at->supellips.R[1][0] = atof(s7);
      at->supellips.R[1][1] = atof(s8);
      at->supellips.R[1][2] = atof(s9);
      at->supellips.R[2][0] = atof(s10);
      at->supellips.R[2][1] = atof(s11);
      at->supellips.R[2][2] = atof(s12);
      at->supellips.a = atof(s13);
      at->supellips.b = atof(s14);
      at->supellips.c = atof(s15);
      at->supellips.n1 = 1.0;
      at->supellips.n2 = 1.0;
      at->supellips.tbeg = atof(s17);
      at->supellips.tend = atof(s18);
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = parsecol(s16,&t,&(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;

    }
  else if (sscanf(L,"%s %s %s %s %s %s %s %s %s %s %s %s @ %s %s %s C[%[^]]", s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16) == 16)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
  	printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	atoi(s5),i, j);*/
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SUPELLIPS;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->supellips.R[0][0] = atof(s4);
      at->supellips.R[0][1] = atof(s5);
      at->supellips.R[0][2] = atof(s6);
      at->supellips.R[1][0] = atof(s7);
      at->supellips.R[1][1] = atof(s8);
      at->supellips.R[1][2] = atof(s9);
      at->supellips.R[2][0] = atof(s10);
      at->supellips.R[2][1] = atof(s11);
      at->supellips.R[2][2] = atof(s12);
      at->supellips.a = atof(s13);
      at->supellips.b = atof(s14);
      at->supellips.c = atof(s15);
      at->supellips.n1 = 1.0;
      at->supellips.n2 = 1.0;
      at->supellips.tbeg = -1.0;
      at->supellips.tend = -1.0;
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = parsecol(s16,&t,&(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;
    }
  else if (sscanf(L,"%s %s %s %s %s %s %s %s %s %s %s %s @ %s %s %s", s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15) == 15)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SUPELLIPS;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->supellips.R[0][0] = atof(s4);
      at->supellips.R[0][1] = atof(s5);
      at->supellips.R[0][2] = atof(s6);
      at->supellips.R[1][0] = atof(s7);
      at->supellips.R[1][1] = atof(s8);
      at->supellips.R[1][2] = atof(s9);
      at->supellips.R[2][0] = atof(s10);
      at->supellips.R[2][1] = atof(s11);
      at->supellips.R[2][2] = atof(s12);
      at->supellips.a = atof(s13);
      at->supellips.b = atof(s14);
      at->supellips.c = atof(s15);
      at->supellips.n1 = 1.0;
      at->supellips.n2 = 1.0;
      at->supellips.tbeg = -1.0;
      at->supellips.tend = -1.0;
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = -1;
      at->common.transp = globset.deftransp;
    }
  else if (sscanf(L,"%s %s %s %s %s %s %s %s %s %s %s %s C[%[^]]]", s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13) == 13)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SUPELLIPS;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->supellips.R[0][0] = atof(s4);
      at->supellips.R[0][1] = atof(s5);
      at->supellips.R[0][2] = atof(s6);
      at->supellips.R[1][0] = atof(s7);
      at->supellips.R[1][1] = atof(s8);
      at->supellips.R[1][2] = atof(s9);
      at->supellips.R[2][0] = atof(s10);
      at->supellips.R[2][1] = atof(s11);
      at->supellips.R[2][2] = atof(s12);
      if (globset.NA)
	{
    	  at->supellips.a = globset.a[a];
	  at->supellips.b = globset.b[a];
	  at->supellips.c = globset.c[a];
	}
      else
	{
	  at->supellips.a = 1.0;
	  at->supellips.b = 1.0;
	  at->supellips.c = 0.5;
	}
      at->supellips.n1 = 1.0;
      at->supellips.n2 = 1.0;
      at->supellips.tbeg = -1.0;
      at->supellips.tend = -1.0;
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = parsecol(s13, &t,&(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;
    }
  else if (sscanf(L,"%s %s %s %s %s %s @ %s %s C[%[^]]", s1, s2, s3, s4, s5, s6, s7, s8, s9) == 9)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_DISK;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->disk.nx = atof(s4);
      at->disk.ny = atof(s5);
      at->disk.nz = atof(s6);
      at->disk.radius = atof(s7);
      at->disk.height = atof(s8);
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = parsecol(s9, &t,&(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;
    }
  else if (sscanf(L,"%s %s %s %s %s %s @ %s %s", s1, s2, s3, s4, s5, s6, s7, s8) == 8)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_DISK;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->disk.nx = atof(s4);
      at->disk.ny = atof(s5);
      at->disk.nz = atof(s6);
      at->disk.radius = atof(s7);
      at->disk.height = atof(s8);
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = -1;
      at->common.transp = globset.deftransp;
    }
  else if (sscanf(L,"%s %s %s @ %s $ %s ", s1, s2, s3, s4, s5) == 5 )
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz= atof(s3);
      at->common.type = MGL_ATOM_SPHERE;
      /*greylLvl[j][i] = colIdxBW[j];*/
      /* default value of grey level*/
      at->sphere.radius = atof(s4);
      at->common.greyLvl = atoi(s5);
      at->common.atcol  = -1;
      at->common.transp = globset.deftransp;
    }
  else if (sscanf(L,"%s %s %s @ %s C[%[^]]] M %[^\n]\n", s1, s2, s3, s4, s5, shl) == 6)
    {
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SPHERE_MSPOT;
      //greylLvl[j][i] = colIdxBW[j];// default value of grey level
      at->sphere_mspot.R[0][0] = 1.0;
      at->sphere_mspot.R[0][1] = 0.0;
      at->sphere_mspot.R[0][2] = 0.0;
      at->sphere_mspot.R[1][0] = 0.0;
      at->sphere_mspot.R[1][1] = 1.0;
      at->sphere_mspot.R[1][2] = 0.0;
      at->sphere_mspot.R[2][0] = 0.0;
      at->sphere_mspot.R[2][1] = 0.0;
      at->sphere_mspot.R[2][2] = 1.0;
      at->sphere_mspot.n1 = 1.0;
      at->sphere_mspot.n2 = 1.0;
      at->sphere_mspot.a = atof(s4);
      at->sphere_mspot.b = atof(s4);
      at->sphere_mspot.c = atof(s4);
      /* loop over spots here */ 
      sl = &at->sphere_mspot.sl;
      *sl=NULL;
      //printf("shl=%s\n", shl);
      for (ss = strtok(shl, ":"); ss; ss = strtok(NULL, ":"))
	{ 
	  newsp = malloc(sizeof(struct spotlst));
	  newsp->next = *sl;
	  sscanf(ss, "%lf %lf %lf %lf C[%[^]]] ", &newsp->n[0], &newsp->n[1], &newsp->n[2], &newsp->spotangle, s10); 
	  //printf("scanning ss=%s\n", ss);
	  newsp->spotcol = parsecol(s10, &t,&(newsp->red), &(newsp->green), &(newsp->blue));
	  *sl = newsp;
	}
      /* ------------------- */
      at->common.transp = t;
      at->common.greyLvl = 0;
      at->common.atcol = parsecol(s5, &t,&(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;
    }
  else if (sscanf(L,"%s %s %s @ %s C[%[^]]] P %s %s %s %s C[%[^]]", s1, s2, s3, s4, s5, s6, s7, s8, s9, s10) == 10)
    {
      /* disegna una sfera di colore s5 con una calotta (che parte dall'angolo s6) di colore s7*/
      /* printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      */
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SPHERE_SPOT;
      //greylLvl[j][i] = colIdxBW[j];// default value of grey level
      at->sphere_spot.nx = atof(s6);
      at->sphere_spot.ny = atof(s7);
      at->sphere_spot.nz = atof(s8);
      at->sphere_spot.R[0][0] = 1.0;
      at->sphere_spot.R[0][1] = 0.0;
      at->sphere_spot.R[0][2] = 0.0;
      at->sphere_spot.R[1][0] = 0.0;
      at->sphere_spot.R[1][1] = 1.0;
      at->sphere_spot.R[1][2] = 0.0;
      at->sphere_spot.R[2][0] = 0.0;
      at->sphere_spot.R[2][1] = 0.0;
      at->sphere_spot.R[2][2] = 1.0;
      at->sphere_spot.n1 = 1.0;
      at->sphere_spot.n2 = 1.0;
      at->sphere_spot.a = atof(s4);
      at->sphere_spot.b = atof(s4);
      at->sphere_spot.c = atof(s4);
      at->sphere_spot.tbeg = atof(s9);
      at->sphere_spot.tend = atof(s9);
      at->sphere_spot.spotcol = parsecol(s10, &t,&(at->sphere_spot.red), &(at->sphere_spot.green), &(at->sphere_spot.blue));
      at->common.transp = t;
      at->common.greyLvl = 0;
      at->common.atcol = parsecol(s5, &t,&(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;
    }
  else if (sscanf(L,"%s %s %s @ %s C[%[^]]", s1, s2, s3, s4, s5) == 5)
    {
      /* printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      */
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SPHERE;
      //greylLvl[j][i] = colIdxBW[j];// default value of grey level
      at->sphere.radius = atof(s4);
      at->common.greyLvl = 0;
      at->common.atcol = parsecol(s5, &t,&(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;
    }
   else if (sscanf(L,"%s %s %s C[%[^]]", s1, s2, s3, s4) == 4 )
    {
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SPHERE;
      if (globset.NA)
	at->sphere.radius = globset.sig[a];
      else
	at->sphere.radius = globset.diameter/2.0;
      at->common.greyLvl = 0;
      at->common.atcol = parsecol(s4, &t,&(at->common.atred), &(at->common.atgreen), &(at->common.atblue));
      at->common.transp = t;
    }
  else if (sscanf(L,"%s %s %s @ %s ", s1, s2, s3, s4) == 4 )
    {
      /* printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      */
      at->common.type = MGL_ATOM_SPHERE;
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      //greylLvl[j][i] = colIdxBW[j];// default value of grey level
      at->sphere.radius = atof(s4);
      at->common.greyLvl = 0;
      at->common.atcol  = -1;
      at->common.transp = globset.deftransp;
    }
  else if (sscanf(L,"%s %s %s %s %s %s %s %s %s %s %s %s ", s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12) == 12)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_SUPELLIPS;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->supellips.R[0][0] = atof(s4);
      at->supellips.R[0][1] = atof(s5);
      at->supellips.R[0][2] = atof(s6);
      at->supellips.R[1][0] = atof(s7);
      at->supellips.R[1][1] = atof(s8);
      at->supellips.R[1][2] = atof(s9);
      at->supellips.R[2][0] = atof(s10);
      at->supellips.R[2][1] = atof(s11);
      at->supellips.R[2][2] = atof(s12);
      if (globset.NA)
	{
    	  at->supellips.a = globset.a[a];
	  at->supellips.b = globset.b[a];
	  at->supellips.c = globset.c[a];
	}
      else
	{
	  at->supellips.a = 1.0;
	  at->supellips.b = 1.0;
	  at->supellips.c = 0.5;
	}
      at->supellips.n1 = 1.0;
      at->supellips.n2 = 1.0;
      at->supellips.tbeg = -1.0;
      at->supellips.tend = -1.0;
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = -1;
      at->common.transp = globset.deftransp;
    }
  else if (sscanf(L,"%s %s %s %s %s %s ", s1, s2, s3, s4, s5, s6) == 6)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      printf("qui\n");
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      at->common.type = MGL_ATOM_DISK;
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      at->disk.nx = atof(s4);
      at->disk.ny = atof(s5);
      at->disk.nz = atof(s6);
      if (globset.NA)
	{
	  at->disk.radius = globset.sig[a];
	  at->disk.height = globset.height;
	}
      else
	{
	  at->disk.radius = globset.diameter/2.0;
	  at->disk.height = globset.height;
	}
      at->common.greyLvl = 0; /*colIdxBW[j];// default value of grey level */
      at->common.atcol  = -1;
      at->common.transp = globset.deftransp;
    }
  else if (sscanf(L,"%s %s %s ", s1, s2, s3) == 3 )
    {
      at->common.rx = atof(s1);
      at->common.ry = atof(s2);
      at->common.rz = atof(s3);
      //greylLvl[j][i] = colIdxBW[j];// default value of grey level
      at->common.type = MGL_ATOM_SPHERE;
      if (globset.NA)
	at->sphere.radius = globset.sig[a];
      else
	at->sphere.radius = globset.diameter/2.0;
      at->common.greyLvl = 0;
      at->common.atcol  = -1;
      at->common.transp = globset.deftransp;
    }
  else
    {
      printf("Line: %s\n", L);
      printf("ERROR: Bad format in input file! Exiting...\n");
      exit(-1);
    }
} 

/* ========================== >>> getColByName <<< ========================= */
int getColByName(const char* name)
{
  /* NOTE: Not optimized search !!!! */
  int nc;
  for (nc = 0; nc < NUMCOLS; ++nc)
    {
      if (!strcmp(mgl_col[nc].name, name))
	{
	  //printf("col name=%s %f %f %f\n", mgl_col[nc].name, mgl_col[nc].rgba[0], mgl_col[nc].rgba[1],
	  //	 mgl_col[nc].rgba[2]);
	  return nc;
	}
    }
  printf("ERROR: Unrecognized color %s!\n", name);
  exit(-1);
}

void add_mol(int cf, int curmol);
void add_frame(int curframes);

void add_atom(int cf, int curmol, int curat)
{
  if (mols==NULL)
    add_frame(0);
  if (mols[cf]==NULL)
    {
      add_mol(cf, 0);
    }
#if 0
  printf("adding atom (%d,%d, %d)\n", cf, curmol, curat);
#endif
  mols[cf][curmol].atom = realloc(mols[cf][curmol].atom, sizeof(atom_u)*(curat+1));
  /*mols[cf][curmol].atom[curat-1]*/
  mols[cf][curmol].numat++;
}

void add_mol(int cf, int curmol)
{ 
  if (mols == NULL)
    add_frame(0);
  mols[cf] = realloc(mols[cf], sizeof(struct molecule)*(curmol+1));
  mols[cf][curmol].atom = NULL;
  mols[cf][curmol].bond = NULL;
  mols[cf][curmol].numbonds = 0;
  mols[cf][curmol].numat = 0;
}

void add_frame(int curframes)
{
  int nf;
  nf = curframes + 1;
  globset.NumMols = realloc(globset.NumMols, sizeof(int)*nf);
  globset.NumMols[nf] = 0;
  atomsList = realloc(atomsList, sizeof(GLuint)*nf);
  mols = realloc(mols, sizeof(struct molecule*)*nf);
  mols[curframes] = NULL;
}

/* =========================== >>> readLine <<< =============================*/
void readLine(FILE* stream, char* L)
{
  if (fscanf(stream, "%[^\n]\n", L) < 1)
    {
      fprintf(stderr,"ERROR: Void line...aborting\n");
      exit(-1);
    }
}

/* ========================== >>> dropSpaces <<< ==========================*/
void dropSpaces(char *S)
{
  char s1[512];
  int i;

  sscanf(S, " %[^\n]", s1);/* drop initial spaces */
//  printf("----->>>>>%s<<<\n", s1);
  /* and now final spaces */
  for (i = strlen(s1) - 1; (i >= 0) && (s1[i] == ' '); --i)
    {
      s1[i] = '\0';  
    }
  strcpy(S, s1);
}

/* ========================== >>> assignCol <<< ===========================*/
void assignCol(char* S, int j)
{
  int colNum;
  char* eptr;

  colNum = (int) strtod(S, &eptr);
  if (eptr == S) /* not a number */
    {
      dropSpaces(S);
      globset.colIdxCol[j] = getColByName(S);
      /* Find the number of the color named 's1'*/
      /*printf("col:%s:, %d\n", S, colIdxCol[j]);
      */
    }
  else
    {
      //printf("Color Number: %d\n", colNum);
      globset.colIdxCol[j] = colNum;
    }
}

/* ========================== >>> pareseLine <<< =========================== */
int parseLine(const char* Line, int* nf, int* i, int *at, int alloc)
{
  /*
    return true if this is a parameter, 0 otherwise
  */
  char parName[1024], parVal[1024], s1[1024], s2[1024], s3[1024], s4[1024], *ns;
  int lett, j, a, nb;
  double defbondthick, t, defbondtransp;
  int defbondcolor; 
  float defbondcolor_red, defbondcolor_green, defbondcolor_blue;
  /* Syntax:
     <parname> : <value>
     where <parname> is of the form ".<name>"
   */
  
  lett =  sscanf(Line, " %[^:#\n] : %[^#\n]", parName, parVal);
  //printf("lett: %d pn:%s|pv: %s\n", lett, parName, parVal);
  if (lett == 0) return 1; /* 2 = comment only */
  if (parName[0] != '.') return 0; // no a parameter line!
  
  /* new frames */
  if (!strcmp(parName, ".newmol"))
    {
      if (!globset.NA)
	{
	  ++(*i);
	  *at = 0;
	  if (alloc)
	    add_mol(*nf,*i);
	}
      return 1; 
    }
  if (!strcmp(parName, ".newframe"))
    {
      if (!((*nf == 0) && (*i == 0)))
	{
	  globset.NumMols[*nf] = *i;
	  ++(*nf);
	  *i = 0;
	  if (alloc)
      	    /* allocate a new frame here */
	    add_frame(*nf);
	}
      return 1;
    }
      
  /* Atoms radi */
  if (!strcmp(parName, ".atomRad"))
    {
      /* Build a string of this type: "%f , %f , ..." with NA '%f' */
      //strcpy(s1, strtok(parVal, ","));
      //printf("radius[0]: %s\n", s1);
      ns = strtok(parVal, ",");
      a = 0;
      if (globset.sig)
	free(globset.sig);
      while(ns)
	{
	  //strcpy(s1, strtok(NULL, ","));
	  globset.sig = realloc(globset.sig,sizeof(double)*(a+1));
	  globset.sig[a] = atof(ns);
	  a++;
	  ns = strtok(NULL, ",");
	}
    
      //printf("---->%f %f\n", sig[0], sig[1]);
      return 1;
    }
  if (!strcmp(parName, ".semiAxes"))
    {
      /* Build a string of this type: "%f , %f , ..." with NA '%f' */
      //strcpy(s1, strtok(parVal, ","));
      //printf("radius[0]: %s\n", s1);
      //printf("parVal: %s\n", parVal);
      a = 0;
      if (globset.a != NULL)
	{
	  free(globset.a);
	  free(globset.b);
	  free(globset.c);
	  globset.a = globset.b = globset.c = NULL;
	}
   
      
      ns = strtok(parVal, ",");
      while(ns)
	{
	  //strcpy(s1, strtok(NULL, ","));
	  globset.a = realloc(globset.a,(a+1)*sizeof(double));
	  globset.b = realloc(globset.b,(a+1)*sizeof(double));
	  globset.c = realloc(globset.c,(a+1)*sizeof(double));
	  //globset.a[a] = atof(ns);
	  //printf("a: %f\n", globset.a[a]);
	  if (sscanf(ns, "%lf %lf %lf ", &(globset.a[a]), &(globset.b[a]), &(globset.c[a])) < 3)
	    {
	      printf("You must supply three axes!\n");
	      exit(-1);
	    }
	  ns = strtok(NULL, ",");
	  a++;
	}

      globset.NA = a;
      //printf("NA: %d\n", globset.NA);
      //printf("---->%f %f\n", sig[0], sig[1]);
      return 1;
    }
  if (!strcmp(parName, ".Vol"))
    {
      globset.L = atof(parVal);
      globset.L = cbrt(globset.L);
      return 1;
    }
  if (!strcmp(parName, ".Bonds"))
    {
      ns = strtok(parVal, ",");
      nb = 0;
      defbondthick = globset.defbondthick;
      defbondcolor = globset.defbondcol;
      defbondtransp = globset.defbondtransp;
      while(ns)
	{
	  nb++;
	  mols[*nf][*i].numbonds = nb;
	  if (alloc)
	    mols[*nf][*i].bond = realloc(mols[*nf][*i].bond, sizeof(bond_s)*nb);
	  
	  if (sscanf(ns, "%[^-]-%[^[][%[^:]:%[^]]]", s1, s2, s3, s4)==4) 
	    {
	      /* atomo-atomo[spessore,colore] */
	      mols[*nf][*i].bond[nb-1].from = atoi(s1);
	      mols[*nf][*i].bond[nb-1].to   = atoi(s2);
	      mols[*nf][*i].bond[nb-1].thickness = atof(s3);
	      mols[*nf][*i].bond[nb-1].color     = parsecol(s4,&t,&(mols[*nf][*i].bond[nb-1].red), &(mols[*nf][*i].bond[nb-1].green), &(mols[*nf][*i].bond[nb-1].blue));
	      mols[*nf][*i].bond[nb-1].transp    = t;
	      defbondthick = atoi(s3);
	      defbondcolor = parsecol(s4, &t,&(defbondcolor_red), &(defbondcolor_green), &(defbondcolor_blue));
	      defbondtransp = t;
	    }
	  else if (sscanf(ns, "[%[^,],%[^]]",s1,s2)==2)
	    {
	      /* [spessore,colore] */
	      defbondthick = atof(s1);
	      defbondcolor = parsecol(s2,&t,&(defbondcolor_red), &(defbondcolor_green), &(defbondcolor_blue));
	      defbondtransp = t;
	    }
	  else if (sscanf(ns, "%[^-]-%s", s1, s2)==2)
	    {
	      /* atomo-atomo */ 
	      mols[*nf][*i].bond[nb-1].from = atoi(s1);
	      mols[*nf][*i].bond[nb-1].to   = atoi(s2);
	      mols[*nf][*i].bond[nb-1].thickness = defbondthick;
	      mols[*nf][*i].bond[nb-1].color     = defbondcolor;
	      mols[*nf][*i].bond[nb-1].red = defbondcolor_red;
	      mols[*nf][*i].bond[nb-1].green = defbondcolor_green;
	      mols[*nf][*i].bond[nb-1].blue = defbondcolor_blue;
	      mols[*nf][*i].bond[nb-1].transp = defbondtransp;
#if 0
	      printf("qui [%s,%s] bondthick:%f\n", s1,s2, mols[*nf][*i].bond[nb-1].thickness );
#endif
	    }
 	  ns = strtok(NULL, ",");
	}
#if 0 
      printf("ns: %s nb:%d nf:%d i:%d mols:%d\n", ns, nb, *nf, *i, mols[*nf][*i].bond[0].to);
#endif
      return 1;
    }
  if (!strcmp(parName, ".fadeMode"))
    {
      globset.fadeMode = atoi(parVal);
      return 1;
    }
  
  if (!strcmp(parName, ".atomCol"))
    {
      ns = strtok(parVal, ",");
      j = 0; 
      if (globset.colIdxCol)
	free(globset.colIdxCol);
      while(ns)
	{
	  strcpy(s1, ns);
	  globset.colIdxCol = realloc(globset.colIdxCol, (j+1)*sizeof(int));
	  assignCol(s1, j);
	  ns = strtok(NULL, ",");
	  j++;
	}
    
      return 1;
    }
  if (!strcmp(parName, ".atomBw"))
    {
      /* Build a string of this type: "%f , %f , ..." with NA '%f' */
      ns = strtok(parVal, ",");
      if (globset.colIdxBw)
	free(globset.colIdxBw);
      a = 0;
      while (ns)
	{
	  globset.colIdxBw = realloc(globset.colIdxBw, (a+1)*sizeof(int));
	  globset.colIdxBw[a] = atoi(ns);
	  ns = strtok(NULL, ",");
	  a++;
	}
      return 1;
    }

  if (!strcmp(parName, ".numAt"))
    {
      globset.NA = atoi(parVal);
      return 1;
    }
  printf("ERROR: Invalid parameter %s!\n", parName);
  exit(-1);
}
void setdefaults_after_fakeread(void)
{
  int a;
  if (globset.NA)
    {
      globset.colIdxCol = malloc(sizeof(int)*globset.NA);
      globset.colIdxBw  = malloc(sizeof(int)*globset.NA);
      for (a = 0; a < globset.NA; ++a)
	{
	  /*  sig[i] = 0.5;*/
	  globset.colIdxCol[a] = a;
	  /* first (and only those ones) 25 grey are quite different in this way */
	  if (a < 25)
	    globset.colIdxBw[a] = a*10;
	  else 
	    globset.colIdxBw[a] = a;
	}
      globset.colIdxCol[0] = 466; /* red1 */
      if (globset.NA > 1)
	globset.colIdxCol[1] = 126; /* green */
      globset.colIdxBw[0] = 120;
      if (globset.NA > 1 )
	globset.colIdxBw[1] = 220;
    }
  if (globset.setdiameter && globset.NA)
    {
      globset.sig = malloc(sizeof(double)*globset.NA);
      for (a = 0; a < globset.NA; ++a)
	{
	  globset.sig[a] = globset.diameter;
	} 
    }

  if (globset.setsemiax && globset.NA)
    {
      globset.a = malloc(sizeof(double)*globset.NA);
      globset.b = malloc(sizeof(double)*globset.NA);
      globset.c = malloc(sizeof(double)*globset.NA);
      for (a = 0; a < globset.NA; ++a)
	{
	  globset.a[a] = globset.sa;
	  globset.b[a] = globset.sb;
	  globset.c[a] = globset.sc;
	} 
    }
}
/* ========================== >>> loadAtomPos <<< ===========================*/
void loadAtomPos(void)
{
  /* DESCRIPTION:
     The file should be an ascii file of the form:
     x0 y0 z0 
     x1 y1 z1 
     .
     .
     .
     where 0 refers to atom 0 and 1 refers to atom 1*/
  FILE* ifs;
  int i, a, nf;
  char line[1024];
  /*printf("Loading file %s\n", inputFile);*/
  i = 0;
  nf = 0;
  a = 0;
 
  ifs = fopen(inputFile, "r"); 
  if (!ifs)
    {
      fprintf(stderr, "ERROR: invalid input file or wrong arguments!\n");
      exit(-1);
    }
  /* fake reading to get number of particles etc. */
  while(!feof(ifs))
    {
#if 0
      fpostmp=ftell(ifs);
#endif
      readLine(ifs, line);

      if (parseLine(line, &nf, &i, &a, 1)) 
	continue;
      #if 0  
        if (first)
	{
	  first=0;
	  fpos = fpostmp;  
	}
#endif
      /*assignAtom(nf, i, 0, line);*/
      
      /*readLine(ifs, line);*/
      /*assignAtom(nf, i, j, line);*/
      
      add_atom(nf, i, a);
      ++a;
      if (globset.NA && a >= globset.NA) 
	{
	  a = 0;
	  i++;
	  add_mol(nf, i);
	}
    }

  /*printf("NumMols:%d", NumMols);*/
  
  if (nf == 0 && globset.NumMols[0] == 0 && mols == NULL)
    {
      printf("ERROR: Presumbly the input file you supplied is void!\n");
      exit(-1);
    }
  
  setdefaults_after_fakeread();
  globset.NumMols[nf] = i+1;
  globset.frameNo = ++nf;
#if 0
   printf("Read %i molecule\n", globset.NumMols[nf]);
   printf("Number of frames: %d\n", globset.frameNo);
#endif
  a = nf = i = 0;
#if 0
  if (fseek(ifs, fpos, SEEK_SET))
    {
      perror("Error seeking the positions file:");
      exit(-1);
    }
#else
  rewind(ifs);
#endif
  while(!feof(ifs))
    {
      readLine(ifs, line);

      if (parseLine(line, &nf, &i, &a, 0)) 
	continue;
      assignAtom(nf, i, a, line);
      ++a;
      if (globset.NA && a >= globset.NA)
	{
	  a = 0;
	  i++;
	}
    }
 fclose(ifs);
  /*
     printf("File with atoms positions successfully read.\n");*/
}


/* ======================== >>> readRGB <<< ==============================*/
void readRGB(void)
{
  int i, rc, gc, bc;
  /* determine the number of color to redad */
  mgl_col = (struct colStruct*) malloc(NUMCOLS * sizeof(struct colStruct));
  for (i=0; i < NUMCOLS; i++)
    {
      sscanf(mglrgb[i], "%d %d %d %[^\n]\n", &rc, &gc, &bc, 
	     mgl_col[i].name);
      //printf("rgb(%d, %d, %d)\n", rc, gc, bc);
      mgl_col[i].rgba[3] = 1.0;
      /* float normalized to [0.0, 1.0] */
      mgl_col[i].rgba[0] = ((float) rc) / 255.0; 
      mgl_col[i].rgba[1] = ((float) gc) / 255.0; 
      mgl_col[i].rgba[2] = ((float) bc) / 255.0; 
    }
}

/* ======================== >>> default_pars <<< ==========================*/
void default_pars(void)
{
  globset.L = 14.0;
  globset.saved_counter = 0;
  globset.saveandquit = 0;
#ifdef MGL_USELIST 
  globset.mgl_uselist = 0;
#endif
  globset.savefile = NULL;
  globset.drawcube = 1;
  globset.sig = NULL;
  /*globset.height = NULL;*/
  globset.a = NULL;
  globset.b = NULL;
  globset.c = NULL;
  globset.setdiameter = 1;
  globset.diameter = 1;
  globset.setsemiax = 0;
  globset.setheight = 0;
  globset.numAt = 0; /* atomi per molecola 0=illimitati a meno che non si usi .newmol*/
  globset.setheight = 0;
  globset.NumMols = NULL;
  globset.Width = 500;
  globset.Height = 500;
  globset.stacks = STACKS;
  globset.slides = SLIDES;
  globset.infos = 1;
  globset.frameNo = 1;
  globset.fadeMode = 2;
  globset.NA = 0;
  globset.axon = 0;
  globset.bw = 0;
  globset.degx = 0.0;
  globset.degy = 0.0;
  globset.degz = 0.0;
  globset.deginc = 5.0;
  globset.viewangle=45.0;
  globset.setvp = 0;
  globset.diameter = 1.0;
  globset.height = 0.2;
  globset.dist = 0.0;
  globset.default_bw=120;
  globset.default_col=466;
  globset.defbondthick = 1.0;
  globset.defbondcol = 465;
  globset.defbondtransp = 1.0;
  globset.deftransp = 1.0;
  globset.nrefresh=1;
  globset.exitDelay=0;/* ritardo in msec prima di fare save and quit */
  globset.light_pos0[0]=globset.light_pos0[1]=globset.light_pos0[2]=150.0;
  globset.light_pos0[3]=0.0;
  globset.light_pos1[0]=globset.light_pos1[1]=globset.light_pos1[2]=150.0;
  globset.light_pos1[3]=0.0;
  globset.twolights=0;
  globset.depthmask=1;
  readRGB();

  setBW();

}


/* =========================== >>> rotatex <<< =========================== */
void special(int k, int x, int y)
{
  switch (k) {
  case GLUT_KEY_UP:
    globset.degx += globset.deginc;
    //printf("degx: %f\n", degx);
    break;
  case GLUT_KEY_DOWN:
    globset.degx -= globset.deginc;
    //printf("degx: %f\n", degx);
    break;
  case GLUT_KEY_LEFT:
    globset.degy -= globset.deginc;
    //printf("degy: %f\n", degy);
    break;
  case GLUT_KEY_RIGHT:
    globset.degy += globset.deginc;
    //printf("degy: %f\n", degy);
    break;
  case GLUT_KEY_PAGE_UP:
  case GLUT_KEY_F1:
    globset.dist += globset.L / 10.0;
    break;
  case GLUT_KEY_PAGE_DOWN:
  case GLUT_KEY_F2:
    globset.dist -= globset.L / 10.0;
    break;

  default:
    return;
  }
  glutPostRedisplay();
}


/* ======================== >>> degadd <<< ================================ */
void key(unsigned char k, int x, int y)
 {
   switch(k) 
     {
     case '1':
       globset.degz += globset.deginc;
       break;
     case '2':
       globset.degz -= globset.deginc;
     break;
     case '+':
       if (globset.deginc > 179.0) break;/* Maximu increment is 180.0 degrees */
       globset.deginc += 1.0; 
       //printf("degree increment: %.0f\n", deginc);
       break;
     case '-':
       if (globset.deginc < 2.0) break; /* Minimum increment is 1.0 degree */
       globset.deginc -= 1.0;
       //printf("deg increment: %.0f\n", deginc);
       break;
     case 'a':
       globset.L += 1.0;
       //printf("L=%.0f\n", L);
       setproj();
       break;
     case 'z':
       if (globset.L > 1.0)
	 globset.L -= 1.0;
       //printf("L=%.0f\n", L);
       setproj();
       break;
     case 'i':
       globset.deginc = -globset.deginc;
       //printf("deginc inverted %.0f to %.0f\n", -deginc, deginc);
       break;
     case 'r':
       globset.degx=0.0;
       globset.degy=0.0;
       globset.degz=0.0;
       //printf("Reset rotations to (0,0,0)\n");
       break;
     case 27: // escape  
       exit(0);
       break;
     case 'p':
       globset.infos = !globset.infos;
       break;
     case 'c':
       globset.bw = !globset.bw;
       //printf("bw: %d\n", bw);     
       break;
     case 'f':
       ++globset.fadeMode;
       if (globset.fadeMode > 3) globset.fadeMode = 1;
       break;
     case 'v':
       globset.axon = !globset.axon;
       myReshape(globset.Width, globset.Height);
       break;
     case 's':
       /* save the current image here */
       save_image();
       if (globset.saved_counter < 256000) 
	 globset.saved_counter++;
     break;
     default: return;

     }
   glutPostRedisplay();
}
/*  Main Loop

 *  Open window with initial window size, title bar, 
 *  RGBA display mode, and handle input events.
 */
int main(int argc, char** argv)
{
  default_pars();

#ifndef SQ_CALC_NORM
  printf("WARNING: numerical gradient will be used for SQ and it may be inaccurate\n");
#endif
  glutInit(&argc, argv);
  args(argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(globset.Width,globset.Height);
  glutCreateWindow("MOLGL by Cristiano De Michele (C) 1998-2010");
  myinit();
  loadAtomPos();
#ifdef MGL_USELIST
  if (globset.mgl_uselist)
    buildAtomsList();
#endif
  glutDisplayFunc(display);
  glutReshapeFunc(myReshape);
  glutKeyboardFunc(key);
  glutSpecialFunc(special);
  /*glutVisibilityFunc(visible);*/
  //glutPostRedisplay();
  glutMainLoop();
  return 0;             /* ANSI C requires main to return int. */
}
