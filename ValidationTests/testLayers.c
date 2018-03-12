/*
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFPermfield.h"
#include "VFCracks.h"
#include "VFWell.h"

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  VFCtx          ctx;
  VFFields       fields;
  PetscReal      ***v_array;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ei,ej,ek;
  int            rank;
  PetscViewer    viewer;
  
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = DMView(ctx.daScalCell,PETSC_VIEWER_STDOUT_WORLD);
  
  ierr = DMDAGetInfo(ctx.daScalCell,NULL,&nx,&ny,&nz,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) ctx.fields->pmult,"layer");CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScalCell,ctx.fields->pmult,&v_array);CHKERRQ(ierr);
  for (ek = zs; ek < zs+zm; ek++) {
    printf ("[%d] ek: %d layer %d, Gc=%f\n",rank,ek,ctx.layer[ek],ctx.matprop[ctx.layer[ek]].Gc);
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
        v_array[ek][ej][ei] = ctx.layer[ek];
      }
    }
  }
  ierr = DMDAVecRestoreArray(ctx.daScalCell,ctx.fields->pmult,&v_array);CHKERRQ(ierr);

  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,"testLayers.vts",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecViewVTKDof(ctx.daScalCell,ctx.fields->pmult,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = VecView(ctx.fields->pmult,PETSC_VIEWER_STDOUT_WORLD);

  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}