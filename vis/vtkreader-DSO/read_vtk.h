#ifndef __read_vtk_H__
#define __read_vtk_H__

/*
 * This file contains functions that read in both fld.vtk and mom.vtk's that are 
 * produced by Pegasus.
 *
 * See comments below for more documentation.
 *
 */


// Field and Moment structures

// Magnetic and electric fields
struct Fields {
  float ***b1;
  float ***b2;
  float ***b3;

  float ***e1;
  float ***e2;
  float ***e3;
};

// Density, velocity (not momentum!) and pressure tensor
struct Moments {
  float ***n;

  float ***u1;
  float ***u2;
  float ***u3;

  float ***p11;
  float ***p12;
  float ***p13;
  float ***p22;
  float ***p23;
  float ***p33;
};


// Output parameter structure

struct Parameters {
  float time; // Time of simulation

  int nx; // Cells in x direction
  int ny; // Cells in y direction 
  int nz; // Cells in z direction
  
  float x0; // Leftmost point in x space
  float y0; // Leftmost point in y space
  float z0; // Leftmost point in z space

  float dx; 
  float dy;
  float dz;
};

/*
 * Read in fld.vtk file that Pegasus outputs.
 * Takes in filename and whether you want the electric field to 
 * be read (otherwise sets e-field pointers to NULL).
 *
 * It also takes pointers to the parameter/field structures where
 * you want the file to be read into.
 *
 * Allocates space in the field point, so that must be freed            
 * using freeFieldStruct.                                                           
 */
void read_fld(const char *fname, const short readEfield, struct Parameters *pOut, struct Fields *fOut);


/*
 * Read in mom.vtk file that Pegasus outputs.
 * Takes in filename and whether you want the pressure tensor to 
 * be read (otherwise sets pressure tensor pointers to NULL).
 *
 * It also takes pointers to the parameter/field structures where
 * you want the file to be read into.
 *
 * Allocates space in the field point, so that must be freed
 * using freeMomStruct.                                                           
 */
void read_mom(const char *fname, const short readPressure, struct Parameters *pOut, 
              struct Moments *mOut);

/*
 * Read in all.vtk file that Pegasus outputs.
 * Takes in filename and whether you want the pressure tensor to 
 * be read (otherwise sets pressure tensor pointers to NULL).
 *
 * It also takes pointers to the parameter/field structures where
 * you want the file to be read into.
 *
 * Allocates space in the field point, so that must be freed
 * using freeMomStruct.                                                           
 */
void read_all(const char *fname, struct Parameters *pOut,
              struct Fields *fOut, struct Moments *mOut);


// Free all pointers in Field struct. DOES NOT FREE STRUCT POINTER ITSELF.
void freeFieldStruct(struct Fields *fld);

// Free all pointers in Mom struct. DOES NOT FREE STRUCT POINTER ITSELF.
void freeMomentStruct(struct Moments *mom);
#endif
