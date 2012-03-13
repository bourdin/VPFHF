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
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal			****coords_array;
	PetscReal			****displ_array;
	Vec					displ_local;
	PetscReal			***vfield_array;
	Vec				    vfield_local;
	PetscReal			****perm_array;
	Vec					perm_local;
    PetscReal			Vol_change_total1 = 0;
    PetscReal			Vol_change_total2 = 0;
	
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
	
	
	if (xs+xm == nx)    xm--;
	if (ys+ym == ny)    ym--;
	if (zs+zm == nz)    zm--;
	
	
    
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = ComputeXYZOpening(&ctx->e3D, ei, ej, ek, hx, hy, hz, displ_array, vfield_array, perm_array);CHKERRQ(ierr);
				
                Vol_change_total1 += perm_array[ek][ej][ei][0]; 
                Vol_change_total2 += perm_array[ek][ej][ei][3]; 
					//			if(ei == nx-5 && ek == 0){
					//				printf("\n opening[%d][%d][%d] = %f\n", ek, ej, ei, perm_array[ek][ej][ei][0]);
					//			}

			}
		}
	}
	
    printf("\n###################################################################\n\n\n");
    printf("#        Volume change calculated gradv.u and volumetric strain     \n");
    printf("#        Volume change using gradv .u = %f\t        \n", Vol_change_total1);
    printf("#        Volume change using volumetric strain = %f\t       \n\n", Vol_change_total2);
	

    
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
    PetscReal		*du_loc;
	PetscReal		*dv_loc;
	PetscReal		*dw_loc;
	

	PetscFunctionBegin;
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&udispl_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&vdispl_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&wdispl_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&dx_vfield_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&dy_vfield_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&dz_vfield_loc);CHKERRQ(ierr);	
    ierr = PetscMalloc(e->ng * sizeof(PetscReal),&du_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&dv_loc);CHKERRQ(ierr);
	ierr = PetscMalloc(e->ng * sizeof(PetscReal),&dw_loc);CHKERRQ(ierr);
    
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
        du_loc[eg] = 0.;
		dv_loc[eg] = 0.;
		dw_loc[eg] = 0.;
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
					
					du_loc[eg] += displ_array[ek+k][ej+j][ei+i][0] * e->dphi[k][j][i][0][eg];
					dv_loc[eg] += displ_array[ek+k][ej+j][ei+i][1] * e->dphi[k][j][i][1][eg];
					dw_loc[eg] += displ_array[ek+k][ej+j][ei+i][2] * e->dphi[k][j][i][2][eg];
				}
			}
		}
	}
	
	for(eg = 0; eg < e->ng; eg++){
		perm_array[ek][ej][ei][0] += (udispl_loc[eg]*dx_vfield_loc[eg] + vdispl_loc[eg]*dy_vfield_loc[eg] + wdispl_loc[eg]*dz_vfield_loc[eg])*e->weight[eg]; /* k_{11}    */
		perm_array[ek][ej][ei][1] += dx_vfield_loc[eg]*e->weight[eg]; /* k_{12}    */
		perm_array[ek][ej][ei][2] += dy_vfield_loc[eg]*e->weight[eg];   /* k_{13}    */
        perm_array[ek][ej][ei][3] += (du_loc[eg]+dv_loc[eg]+dw_loc[eg])*e->weight[eg];  /* k_{22}    */
        perm_array[ek][ej][ei][4] += dz_vfield_loc[eg]*e->weight[eg];   /* k_{23}    */
		
    }
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "NodeToCellInterpolation"
extern PetscErrorCode NodeToCellInterpolation(DM dm, Vec node_vec, Vec cell_vec)
{
	PetscErrorCode ierr;
	PetscInt		dof;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscInt        k, j, i;
	PetscInt        c;
	PetscReal   ****node_arraydof;
	PetscReal   ****cell_arraydof;
	PetscReal   ***node_array;
	PetscReal   ***cell_array;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(dm,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	printf("\nThis is the dof = %d\n", dof);
	
	
	if (dof == 1){
		printf("\nInside dof equal to one one\n");
		ierr = DMDAVecGetArray(dm, node_vec,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecGetArray(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		printf("\nInside dof greater than one\n");
		ierr = DMDAVecGetArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
		ierr = DMDAVecGetArrayDOF(dm, cell_vec,&cell_arraydof);CHKERRQ(ierr);   
	}
	
	
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
				/*	(Corner 1): ek = 0; ej = 0; ei = 0	*/
				if(ek == 0 && ej == 0 && ei == 0){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek+k][ej+j][ei+i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}	
							}
						}
					}
				}
				/*  (Corner 2):	ek = 0; ej = 0; ei = nx-1	*/
				else if(ek == 0 && ej == 0 && ei == nx-1){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek+k][ej+j][ei-i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek+k][ej+j][ei-i][c];
									}
								}	
							}
						}
					}
				}
				/*	(Corner 3):	ek = 0; ej = ny-1; ei = 0	*/
				else if(ek == 0 && ej == ny-1 && ei == 0){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek+k][ej-j][ei+i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek+k][ej-j][ei+i][c];
									}
								}	
							}
						}
					}
				}								
				/*	(Corner 4):	ek = nz-1; ej = 0; ei = 0		*/
				else if(ek == nz-1 && ej == 0 && ei == 0){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek-k][ej+j][ei+i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek-k][ej+j][ei+i][c];
									}
								}	
							}
						}
					}
				}				
				/*	(Corner 5):	ek = 0; ej = ny-1; ei = nx-1	*/
				else if(ek == 0 && ej == ny-1 && ei == nx-1){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek+k][ej-j][ei-i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek+k][ej-j][ei-i][c];
									}
								}	
							}
						}
					}
				}
				
				/*	(Corner 6):	ek = nz-1; ej = 0; ei = nx-1	*/
				else if(ek == nz-1 && ej == 0 && ei == nx-1){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek-k][ej+j][ei-i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek-k][ej+j][ei-i][c];										
									}
								}	
							}
						}
					}
				}
				/*	(Corner 7):	ek = nz-1; ej = ny-1; ei = 0	*/
				else if(ek == nz-1 && ej == ny-1 && ei == 0){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek-k][ej-j][ei+i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek-k][ej-j][ei+i][c];											
									}
								}	
							}
						}
					}
				}
				/*	(Corner 8):	ek = nz-1; ej = ny-1; ei = nx-1		*/
				else if(ek == nz-1 && ej == ny-1 && ei == nx-1){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek-k][ej-j][ei-i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek-k][ej-j][ei-i][c];
									}
								}	
							}
						}
					}
				}				
				/*	(Edge 1):	ei = 0; ek = 0		*/
				else if(ei == 0 && ej != 0 && ej != ny-1 && ek == 0){
					for(k = 0; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek+k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}
				
				/*	(Edge 2):	ei = nx-1;	ek = 0	*/
				else if(ei == nx-1 && ej != 0 && ej != ny-1 && ek == 0){
					for(k = 0; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek+k][ej+j][ei-i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek+k][ej+j][ei-i][c];
									}
								}
							}
						}
					}
				}				
				
				/*	(Edge 3):	ej = 0;		ek = 0		*/
				else if(ei != 0 && ei != nx-1 && ej == 0 && ek == 0){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek+k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}
				/*	(Edge 4):	ej = ny-1;		ek = 0		*/
				else if(ei != 0 && ei != nx-1 && ej == ny-1 && ek == 0){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek+k][ej-j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek+k][ej-j][ei+i][c];
									}
								}
							}
						}
					}
				}				
				/*	(Edge 5):	ei = 0; ek = nz-1	*/
				else if(ei == 0 && ej != 0 && ej != ny-1 && ek == nz-1){
					for(k = 0; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek-k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek-k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}
				
				/*	(Edge 6):	ei = nx-1;	ek = nz-1	*/
				else if(ei == nx-1 && ej != 0 && ej != ny-1 && ek == nz-1){
					for(k = 0; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek-k][ej+j][ei-i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek-k][ej+j][ei-i][c];
									}
								}
							}
						}
					}
				}				
				
				/*	(Edge 7):	ej = 0;		ek = nz-1	*/
				else if(ei != 0 && ei != nx-1 && ej == 0 && ek == nz-1){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek-k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek-k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}
				/*	(Edge 8):	ej = ny-1;		ek = nz-1	*/
				else if(ei != 0 && ei != nx-1 && ej == ny-1 && ek == nz-1){
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek-k][ej-j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek-k][ej-j][ei+i][c];
									}
								}
							}
						}
					}
				}				
				/*	(Edge 9):	ei = 0;		ej = 0		*/
				else if(ei == 0 && ej == 0 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek+k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}
				/*	(Edge 10):	ei = 0;		ej = ny-1	*/
				else if(ei == 0 && ej == ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek+k][ej-j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek+k][ej-j][ei+i][c];
									}
								}
							}
						}
					}
				}				
				/*	(Edge 11):	ei = nx-1;	ej = 0		*/
				else if(ei == nx-1 && ej == 0 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek+k][ej+j][ei-i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek+k][ej+j][ei-i][c];
									}
								}
							}
						}
					}
				}				
				/*	(Edge 12):	ei = nx-1;	ej = ny-1		*/
				else if(ei == nx-1 && ej == ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./12.*node_array[ek+k][ej-j][ei-i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./12.*node_arraydof[ek+k][ej-j][ei-i][c];
									}
								}
							}
						}
					}
				}	
				
				/*	(Side 1):	ei = 0		*/
				else if(ei == 0 && ej != 0 && ej != ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./18.*node_array[ek+k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./18.*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}	
				/*	(Side 2):	ei = nx-1		*/
				else if(ei == nx-1 && ej != 0 && ej != ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./18.*node_array[ek+k][ej+j][ei-i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./18.*node_arraydof[ek+k][ej+j][ei-i][c];
									}
								}
							}
						}
					}
				}
				/*	(Side 3):	ej = 0		*/
				else if(ei != 0 && ei != nx-1 && ej == 0 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./18.*node_array[ek+k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./18.*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}
				/*	(Side 4):	ej = ny-1		*/
				else if(ei != 0 && ei != nx-1 && ej == ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./18.*node_array[ek+k][ej-j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./18.*node_arraydof[ek+k][ej-j][ei+i][c];
									}
								}
							}
						}
					}
				}
				/*	(Side 5):	ek = 0		*/
				else if(ei != 0 && ei != nx-1 && ej != 0 && ej != ny-1 && ek == 0 ){
					for(k = 0; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./18.*node_array[ek+k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./18.*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}
				/*	(Side 6):	ek = nz-1	*/
				else if(ei != 0 && ei != nx-1 && ej != 0 && ej != ny-1 && ek == nz-1 ){
					for(k = 0; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./18.*node_array[ek-k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./18.*node_arraydof[ek-k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}
				/*	Interior of domain		*/
				else if (ei != 0 && ei != nx-1 && ej != 0 && ej != ny-1 && ek != 0 && ek != nz-1) {
					for(k = -1; k < 2; k++){
						for(j = -1; j < 2; j++){
							for(i = -1; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./27.*node_array[ek+k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./27.*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}								
				
			}
		}
	}

	if(dof == 1){
		ierr = DMDAVecRestoreArray(dm, node_vec,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecRestoreArray(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		ierr = DMDAVecRestoreArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
		ierr = DMDAVecRestoreArrayDOF(dm, cell_vec,&cell_arraydof);CHKERRQ(ierr);    
	}
    PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "NodeToCellInterpolation1"
extern PetscErrorCode NodeToCellInterpolation1(DM dm, Vec node_vec, Vec cell_vec)
{
	PetscErrorCode ierr;
	PetscInt		dof;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscInt        k, j, i;
	PetscInt        c;
	PetscReal   ****node_arraydof;
	PetscReal   ****cell_arraydof;
	PetscReal   ***node_array;
	PetscReal   ***cell_array;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(dm,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	printf("\nThis is the dof = %d\n", dof);
	
	
	if (dof == 1){
		printf("\nInside dof equal to one one\n");
		ierr = DMDAVecGetArrayDOF(dm, node_vec,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecGetArrayDOF(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		printf("\nInside dof greater than one\n");
		ierr = DMDAVecGetArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
		ierr = DMDAVecGetArrayDOF(dm, cell_vec,&cell_arraydof);CHKERRQ(ierr);   
	}
	
	if (xs+xm == nx) xm--;
	if (ys+ym == ny) ym--;
	if (zs+zm == nz) zm--;

	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
					for(k = 0; k < 2; k++){
						for(j = 0; j < 2; j++){
							for(i = 0; i < 2; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 0.125*node_array[ek+k][ej+j][ei+i];
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 0.125*node_arraydof[ek+k][ej+j][ei+i][c];
								}
							}	
						}
					}
				}
			}
		}
	}
	if(dof == 1){
		ierr = DMDAVecRestoreArray(dm, node_vec,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecRestoreArray(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		ierr = DMDAVecRestoreArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
		ierr = DMDAVecRestoreArrayDOF(dm, cell_vec,&cell_arraydof);CHKERRQ(ierr);    
	}
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CellToNodeInterpolation"
extern PetscErrorCode CellToNodeInterpolation(DM dm, Vec node_vec, Vec cell_vec)
{
	PetscErrorCode ierr;
	PetscInt		dof;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscInt        k, j, i;
	PetscInt        c;
	PetscReal   ****node_arraydof;
	PetscReal   ****cell_arraydof;
	PetscReal   ***node_array;
	PetscReal   ***cell_array;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(dm,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	printf("\nThis is the dof = %d\n", dof);
	
	
	if (dof == 1){
		printf("\ndof = one............................\n");
		ierr = DMDAVecGetArray(dm, node_vec,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecGetArray(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		ierr = DMDAVecGetArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
		ierr = DMDAVecGetArrayDOF(dm, cell_vec,&cell_arraydof);CHKERRQ(ierr);   
	}
	
	
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
				/*	(Corner 1): ek = 0; ej = 0; ei = 0	*/
				if(ek == 0 && ej == 0 && ei == 0){
					if(dof == 1)
					{
						cell_array[ek][ej][ei] = node_array[ek][ej][ei];
					}	
					else
					{
						for(c = 0; c < dof; c++){
						cell_arraydof[ek][ej][ei][c] = node_arraydof[ek][ej][ei][c];
						}
					}
				}
				/*  (Corner 2):	ek = 0; ej = 0; ei = nx-1	*/
				else if(ek == 0 && ej == 0 && ei == nx-1){
					if(dof == 1)
					{
						cell_array[ek][ej][ei] = node_array[ek][ej][ei-1];		/* ei-1 since nx-1 is padded with zeros*/
					}
					else
					{
						for(c = 0; c < dof; c++){
							cell_arraydof[ek][ej][ei][c] = node_arraydof[ek][ej][ei-1][c];
						}
					}
				}
				/*	(Corner 3):	ek = 0; ej = ny-1; ei = 0	*/
				else if(ek == 0 && ej == ny-1 && ei == 0){
					if(dof == 1)
					{
						cell_array[ek][ej][ei] = node_array[ek][ej-1][ei];
					}
					else
					{
						for(c = 0; c < dof; c++){
							cell_arraydof[ek][ej][ei][c] = node_arraydof[ek][ej-1][ei][c];						
						}
					}
				}								
				/*	(Corner 4):	ek = nz-1; ej = 0; ei = 0		*/
				else if(ek == nz-1 && ej == 0 && ei == 0){
					if(dof == 1)
					{
						cell_array[ek][ej][ei] = node_array[ek-1][ej][ei];
					}
					else
					{
						for(c = 0; c < dof; c++){
							cell_arraydof[ek][ej][ei][c] = node_arraydof[ek-1][ej][ei][c];
						}
					}	
				}				
				/*	(Corner 5):	ek = 0; ej = ny-1; ei = nx-1	*/
				else if(ek == 0 && ej == ny-1 && ei == nx-1){
					if(dof == 1)
					{
						cell_array[ek][ej][ei] = node_array[ek][ej-1][ei-1];
					}
					else
					{
						for(c = 0; c < dof; c++){
							cell_arraydof[ek][ej][ei][c] = node_arraydof[ek][ej-1][ei-1][c];
						}
					}	
				}
				/*	(Corner 6):	ek = nz-1; ej = 0; ei = nx-1	*/
				else if(ek == nz-1 && ej == 0 && ei == nx-1){
					if(dof == 1)
					{
						cell_array[ek][ej][ei] = node_array[ek-1][ej][ei];
					}
					else
					{
						for(c = 0; c < dof; c++){
							cell_arraydof[ek][ej][ei][c] = node_arraydof[ek-1][ej][ei][c];										
						}
					}	
				}
				/*	(Corner 7):	ek = nz-1; ej = ny-1; ei = 0	*/
				else if(ek == nz-1 && ej == ny-1 && ei == 0){
					if(dof == 1)
					{
						cell_array[ek][ej][ei] = node_array[ek-1][ej-1][ei];
					}
					else
					{
						for(c = 0; c < dof; c++){
							cell_arraydof[ek][ej][ei][c] = node_arraydof[ek-1][ej-1][ei][c];											
						}
					}	
				}
				/*	(Corner 8):	ek = nz-1; ej = ny-1; ei = nx-1		*/
				else if(ek == nz-1 && ej == ny-1 && ei == nx-1){
					if(dof == 1)
					{
						cell_array[ek][ej][ei] = node_array[ek-1][ej-1][ei-1];
					}
					else
					{
						for(c = 0; c < dof; c++){
							cell_arraydof[ek][ej][ei][c] = node_arraydof[ek-1][ej-1][ei-1][c];
						}
					}	
				}				
				/*	(Edge 1):	ei = 0; ek = 0		*/
				else if(ei == 0 && ej != 0 && ej != ny-1 && ek == 0){
					for(j = -1; j < 1; j++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek][ej+j][ei]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek][ej+j][ei][c];
								}
							}					
						}
					}
				/*	(Edge 2):	ei = nx-1;	ek = 0	*/
				else if(ei == nx-1 && ej != 0 && ej != ny-1 && ek == 0){
					for(j = -1; j < 1; j++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./1.*node_array[ek][ej+j][ei-1]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./1.*node_arraydof[ek][ej+j][ei-1][c];
								}
							}
						}
					}			
				
				/*	(Edge 3):	ej = 0;		ek = 0		*/
				else if(ei != 0 && ei != nx-1 && ej == 0 && ek == 0){
					for(i = -1; i < 1; i++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek][ej][ei+i]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek][ej][ei+i][c];
								}
							}
						}
					}
				/*	(Edge 4):	ej = ny-1;		ek = 0		*/
				else if(ei != 0 && ei != nx-1 && ej == ny-1 && ek == 0){
					for(i = -1; i < 1; i++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek][ej-1][ei+i]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek][ej-1][ei+i][c];
								}
							}
						}
					}
				/*	(Edge 5):	ei = 0; ek = nz-1	*/
				else if(ei == 0 && ej != 0 && ej != ny-1 && ek == nz-1){
					for(j = -1; j < 1; j++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek-1][ej+j][ei]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek-1][ej+j][ei][c];
								}
							}
						}
					}
				/*	(Edge 6):	ei = nx-1;	ek = nz-1	*/
				else if(ei == nx-1 && ej != 0 && ej != ny-1 && ek == nz-1){
						for(j = -1; j < 1; j++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./2.*node_array[ek-1][ej+j][ei-1]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek-1][ej+j][ei-1][c];
									}
								}
							}
						}
				/*	(Edge 7):	ej = 0;		ek = nz-1	*/
				else if(ei != 0 && ei != nx-1 && ej == 0 && ek == nz-1){
					for(i = -1; i < 1; i++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek-1][ej][ei+i]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek-1][ej][ei+i][c];
								}
							}
						}
					}
				/*	(Edge 8):	ej = ny-1;		ek = nz-1	*/
				else if(ei != 0 && ei != nx-1 && ej == ny-1 && ek == nz-1){
					for(i = -1; i < 1; i++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek-1][ej-1][ei+i]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek-1][ej-1][ei+i][c];
								}
							}
						}
					}		
				/*	(Edge 9):	ei = 0;		ej = 0		*/
				else if(ei == 0 && ej == 0 && ek != 0 && ek != nz-1){
					for(k = -1; k < 1; k++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek+k][ej][ei]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek+k][ej][ei][c];
								}
							}
						}
					}
				/*	(Edge 10):	ei = 0;		ej = ny-1	*/
				else if(ei == 0 && ej == ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 1; k++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek+k][ej-1][ei]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek+k][ej-1][ei][c];
								}
							}
						}
					}				
				/*	(Edge 11):	ei = nx-1;	ej = 0		*/
				else if(ei == nx-1 && ej == 0 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek+k][ej][ei-1]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek+k][ej][ei-1][c];
								}
							}
						}
					}		
				/*	(Edge 12):	ei = nx-1;	ej = ny-1		*/
				else if(ei == nx-1 && ej == ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 2; k++){
						if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./2.*node_array[ek+k][ej-1][ei-1]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./2.*node_arraydof[ek+k][ej-1][ei-1][c];
								}
							}
						}
					}
				/*	(Side 1):	ei = 0		*/
				else if(ei == 0 && ej != 0 && ej != ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 1; k++){
						for(j = -1; j < 1; j++){
							if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./4.*node_array[ek+k][ej+j][ei]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./4.*node_arraydof[ek+k][ej+j][ei][c];
								}
							}
						}
					}
				}
				/*	(Side 2):	ei = nx-1		*/
				else if(ei == nx-1 && ej != 0 && ej != ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 1; k++){
						for(j = -1; j < 1; j++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./4.*node_array[ek+k][ej+j][ei-1]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./4.*node_arraydof[ek+k][ej+j][ei-1][c];
									}
								}
							}
						}
					}
				/*	(Side 3):	ej = 0		*/
				else if(ei != 0 && ei != nx-1 && ej == 0 && ek != 0 && ek != nz-1){
					for(k = -1; k < 1; k++){
						for(i = -1; i < 1; i++){
							if(dof == 1)
							{
								cell_array[ek][ej][ei] += 1./4.*node_array[ek+k][ej][ei+i]; 
							}
							else
							{
								for(c = 0; c < dof; c++){
									cell_arraydof[ek][ej][ei][c] += 1./4.*node_arraydof[ek+k][ej][ei+i][c];
								}
							}
						}
					}
				}
				/*	(Side 4):	ej = ny-1		*/
				else if(ei != 0 && ei != nx-1 && ej == ny-1 && ek != 0 && ek != nz-1){
					for(k = -1; k < 1; k++){
							for(i = -1; i < 1; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./4.*node_array[ek+k][ej-1][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./4.*node_arraydof[ek+k][ej-1][ei+i][c];
									}
								}
							}
						}
					}
				/*	(Side 5):	ek = 0		*/
				else if(ei != 0 && ei != nx-1 && ej != 0 && ej != ny-1 && ek == 0 ){
						for(j = -1; j < 1; j++){
							for(i = -1; i < 1; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./4.*node_array[ek][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./4.*node_arraydof[ek][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				/*	(Side 6):	ek = nz-1	*/
				else if(ei != 0 && ei != nx-1 && ej != 0 && ej != ny-1 && ek == nz-1 ){
						for(j = -1; j < 1; j++){
							for(i = -1; i < 1; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./4.*node_array[ek-1][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./4.*node_arraydof[ek-1][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				/*	Interior of domain		*/
				else if (ei != 0 && ei != nx-1 && ej != 0 && ej != ny-1 && ek != 0 && ek != nz-1) {
					for(k = -1; k < 1; k++){
						for(j = -1; j < 1; j++){
							for(i = -1; i < 1; i++){
								if(dof == 1)
								{
									cell_array[ek][ej][ei] += 1./8.*node_array[ek+k][ej+j][ei+i]; 
								}
								else
								{
									for(c = 0; c < dof; c++){
										cell_arraydof[ek][ej][ei][c] += 1./8.*node_arraydof[ek+k][ej+j][ei+i][c];
									}
								}
							}
						}
					}
				}								
				
			}
		}
	}
	PetscReal sum = 0,  sum1 = 0;
	for(k = 0; k < zs+zm; k++){
		for(j = 0; j < ys+ym; j++){
			for(i = 0; i < xs+xm; i++){
				sum += node_array[k][j][i];
				sum1 += cell_array[k][j][i];
					//printf("\nnode[%d][%d][%d] = %f\t cell[%d][%d][%d] = %f\n",k,j,i, node_array[k][j][i],k,j,i, cell_array[k][j][i]);
			}
		}
	}
	printf("\nsum = %f\t sum1 = %f\n", sum, sum1);
	
	if(dof == 1){
		ierr = DMDAVecRestoreArray(dm, node_vec,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecRestoreArray(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		ierr = DMDAVecRestoreArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
		ierr = DMDAVecRestoreArrayDOF(dm, cell_vec,&cell_arraydof);CHKERRQ(ierr);    
	}
    PetscFunctionReturn(0);
}