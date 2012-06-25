/*
   VFFlow_DarcySteadyState.h
   A mixed finite elements Darcy solver based on the method presented in
   A. Masud and T.J.R. Hughes. 
      "A stabilized mixed finite element method for Darcy flow. Computer Methods 
      in Applied Mechanics and Engineering", 191(39–40):4341 – 4370, 2002.
   K.B. Nakshatrala, D.Z. Turner, K.D. Hjelmstad, and A. Masud. 
      "A stabilized mixed finite element method for Darcy flow based on a 
      multiscale decomposition of the solution". Computer Methods in Applied 
      Mechanics and Engineering, 195(33–36):4036 – 4049, 2006.    

   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/

#ifndef VFPERMFIELD_H
#define VFPERMFIELD_H


/* 
  Rename and check if all these need to be public
*/
extern PetscErrorCode CellToNodeInterpolation(DM dm, Vec node_vec, Vec cell_vec, VFCtx *ctx);
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ****displ_array, PetscReal ***vfield_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e);
extern PetscErrorCode VolumetricCrackOpening(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields);

#endif 
