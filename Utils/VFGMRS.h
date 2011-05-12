#ifndef VFGMRS_H
#define VFGMRS_H
extern PetscErrorCode VecGMRSToPetsc(PetscReal *GMRS_array,Vec x);
extern PetscErrorCode VecPetscToGMRS(Vec x,PetscReal *GMRS_array);
#endif /* VFGMRS_H */
