#include "petsc.h"
#include "VF_CInterface.h"

#undef __FUNCT__
#define __FUNCT__ "fakepmult"
/*
   fakepmult: generates a kake pressure multiplie field
*/
extern PetscErrorCode fakepmult(VFCtx *ctx) {
	PetscErrorCode ierr;
	PetscInt       i,j,k,xs,xm,ys,ym,zs,zm;
	PetscReal      ***pmultarray;
	PetscReal      r;
	
	PetscFunctionBegin;
	ierr = DAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DAVecGetArray(ctx->daScal,ctx->pmult,&pmultarray);CHKERRQ(ierr);
  for (k=zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				r = i * i / 100;
				pmultarray[k][j][i] = .3 + .1 * r;
			}
		}
	}

	ierr = DAVecRestoreArray(ctx->daScal,ctx->pmult,&pmultarray);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
