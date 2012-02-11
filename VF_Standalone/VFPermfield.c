/*
   VFFlow_DarcySteadyState.c
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"

/*
   VFPermfield.c
    
   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/
#undef __FUNCT__
#define __FUNCT__ "CrackOpeningDisplacement"
extern PetscErrorCode CrackOpeningDisplacement(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
	PetscInt		ek, ej, ei;
	PetscInt		xs,xm,nx;
	PetscInt		ys,ym,ny;
	PetscInt		zs,zm,nz;
	PetscReal		hx,hy,hz;
	PetscReal		****coords_array;
	PetscReal		****displ_array;
	Vec				displ_local;
	PetscReal		***vfield_array;
	Vec				vfield_local;
	PetscReal		****perm_array;
	Vec				perm_local;
	
  	PetscFunctionBegin;
	
  ierr = DMDAGetInfo(ctx->daFlow,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daVect,&displ_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,displ_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,displ_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,displ_local,&displ_array);CHKERRQ(ierr); 
	
	ierr = DMGetLocalVector(ctx->daScal,&vfield_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,vfield_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,vfield_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,vfield_local,&vfield_array);CHKERRQ(ierr); 	
	
	/*Setting all permeability field to zero*/
	ierr = VecSet(fields->vfperm,0.);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr); 
	

	if (xs+xm == nx) xm--;
	if (ys+ym == ny) ym--;
	if (zs+zm == nz) zm--;
	
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = ComputeXYZOpening(&ctx->e3D, ei, ej, ek, hx, hy, hz, displ_array, vfield_array, perm_array);CHKERRQ(ierr);
				
	//			if(ei == nx-5 && ek == 0){
	//				printf("\n opening[%d][%d][%d] = %f\n", ek, ej, ei, perm_array[ek][ej][ei][0]);
	//			}
			
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,displ_local,&displ_array);CHKERRQ(ierr); 
	ierr = DMDAVecRestoreArray(ctx->daScal,vfield_local,&vfield_array);CHKERRQ(ierr); 	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr); 
	
	ierr = DMLocalToGlobalBegin(ctx->daVFperm,perm_local,ADD_VALUES,fields->vfperm);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daVFperm,perm_local,ADD_VALUES,fields->vfperm);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
	
#undef __FUNCT__
#define __FUNCT__ "ComputeXYZOpening"
extern PetscErrorCode ComputeXYZOpening(CartFE_Element3D *e, PetscInt ei, PetscInt ej, PetscInt ek, PetscReal hx, PetscReal hy, PetscReal hz, PetscReal ****displ_array, PetscReal ***vfield_array, PetscReal ****perm_array)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k;
	PetscInt		eg;
	PetscReal		*dx_vfield_loc;
	PetscReal		*dy_vfield_loc;
	PetscReal		*dz_vfield_loc;
	PetscReal		*udispl_loc;
	PetscReal		*vdispl_loc;
	PetscReal		*wdispl_loc;
	PetscReal		dnorm;

	PetscFunctionBegin;
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&udispl_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&vdispl_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&wdispl_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&dx_vfield_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&dy_vfield_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&dz_vfield_loc);CHKERRQ(ierr);

	perm_array[ek][ej][ei][0] = 0.;	
	perm_array[ek][ej][ei][1] = 0.;	
	perm_array[ek][ej][ei][2] = 0.;	
	perm_array[ek][ej][ei][3] = 0.;
	perm_array[ek][ej][ei][4] = 0.;
	perm_array[ek][ej][ei][5] = 0.;
	
	for (eg = 0; eg < e->ng; eg++){
		udispl_loc[eg] = 0.;
		vdispl_loc[eg] = 0.;
		wdispl_loc[eg] = 0.;
		dx_vfield_loc[eg] = 0.;
		dy_vfield_loc[eg] = 0.;
		dz_vfield_loc[eg] = 0.;
	}
	
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
					udispl_loc[eg] += displ_array[ek+k][ej+j][ei+i][0] * e->phi[k][j][i][eg];
					vdispl_loc[eg] += displ_array[ek+k][ej+j][ei+i][1] * e->phi[k][j][i][eg];
					wdispl_loc[eg] += displ_array[ek+k][ej+j][ei+i][2] * e->phi[k][j][i][eg];
					dx_vfield_loc[eg] += vfield_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][0][eg];
					dy_vfield_loc[eg] += vfield_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][1][eg];
					dz_vfield_loc[eg] += vfield_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][2][eg];
				}
			}
		}
		/* Normalizing gradient of the v-field */
		dnorm = sqrt(pow(dx_vfield_loc[eg],2)+pow(dy_vfield_loc[eg],2)+pow(dz_vfield_loc[eg],2));
		if(dnorm != 0.){
			dx_vfield_loc[eg] = dx_vfield_loc[eg]/dnorm;
			dy_vfield_loc[eg] = dy_vfield_loc[eg]/dnorm;
			dz_vfield_loc[eg] = dz_vfield_loc[eg]/dnorm;
		}

	}

	for(eg = 1; eg < e->ng; eg++){
		perm_array[ek][ej][ei][0] += (udispl_loc[eg]*dx_vfield_loc[eg] + vdispl_loc[eg]*dy_vfield_loc[eg] + wdispl_loc[eg]*dz_vfield_loc[eg])*e->weight[eg];
		perm_array[ek][ej][ei][1] += dx_vfield_loc[eg]*e->weight[eg];
		perm_array[ek][ej][ei][2] += dy_vfield_loc[eg]*e->weight[eg];
		perm_array[ek][ej][ei][3] += dz_vfield_loc[eg]*e->weight[eg];
	}
	PetscFunctionReturn(0);
}

