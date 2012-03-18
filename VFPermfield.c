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
 Compute the total crack opening displacement (crack volume) by integrating u.\nabla v 
 over the entire domain
 
 (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
 */

#undef __FUNCT__
#define __FUNCT__ "CrackOpeningDisplacement"
/*
  Comments:
    Rename TotalCrackOpening3D or something like that. 
    Make calls to VolumetricCrackOpening3D_local below
    Make sure to call MPI_AllReduce to summ across all MPI processes
*/
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
	Vec					  displ_local;
	PetscReal			***vfield_array;
	Vec				    vfield_local;
	PetscReal			****perm_array;
	Vec					perm_local;
    PetscReal			MyVol_change_total1 = 0;
    PetscReal			MyVol_change_total2 = 0;
   
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
				
                MyVol_change_total1 += perm_array[ek][ej][ei][0]; 
                MyVol_change_total2 += perm_array[ek][ej][ei][3]; 
					//			if(ei == nx-5 && ek == 0){
					//				printf("\n opening[%d][%d][%d] = %f\n", ek, ej, ei, perm_array[ek][ej][ei][0]);
					//			}

			}
		}
	}
	ierr = MPI_Reduce(&MyVol_change_total1,&Vol_change_total1,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Reduce(&MyVol_change_total2,&Vol_change_total2,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
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
/*
  Comments:
    rename VolumetricCrackOpening_local. Return a C array containing the value of u.\nabla v at each of the gauss point of the local element
    the interface should be :
    PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *VolumetricCrackOpening_local,
                                                  PetscReal ****u_array,
                                                  PetscReal ***v_array,
                                                  PetscInt ek,
                                                  PetscInt ej,
                                                  PetscInt ei,
                                                  CartFE_Element3D *e);
    in order to match that of ElasticEnergyDensity3D_localin VFU.c
*/
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
	PetscReal		dx_vfield = 0;
	PetscReal		dy_vfield = 0;
	PetscReal		dz_vfield = 0;
	PetscReal		vol = 0;
	PetscReal		n_area;
	

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
			//		perm_array[ek][ej][ei][2] += dy_vfield_loc[eg]*e->weight[eg];   /* k_{13}    */
        perm_array[ek][ej][ei][3] += (du_loc[eg]+dv_loc[eg]+dw_loc[eg])*e->weight[eg];  /* k_{22}    */
			//   perm_array[ek][ej][ei][4] += dz_vfield_loc[eg]*e->weight[eg];   /* k_{23}    */
		dx_vfield += dx_vfield_loc[eg]*e->weight[eg]; 
		dy_vfield += dy_vfield_loc[eg]*e->weight[eg]; 
		dz_vfield += dz_vfield_loc[eg]*e->weight[eg]; 
		vol += e->weight[eg];
    }
	dx_vfield = dx_vfield/vol;
	dy_vfield = dy_vfield/vol;
	dz_vfield = dz_vfield/vol;
	perm_array[ek][ej][ei][4] = n_area = dx_vfield*hz*hy+dy_vfield*hx*hz+dz_vfield*hx*hy;
	if(n_area < 0)
		n_area = -1*n_area;
		
	perm_array[ek][ej][ei][2] = perm_array[ek][ej][ei][0]/vol;   /* k_{13}    */

	perm_array[ek][ej][ei][5] = perm_array[ek][ej][ei][0]/n_area;   /* k_{23}    */

	
    PetscFunctionReturn(0);
}


/*
  Comments:
    loop the other way around:
      for cell = 0 to num_cell
        compute int_cell cell_vec
        compute patch area
        for all vertex v of the current cell
          add vec_node[global[v]] += int_cell cell_vec at vertex v    
          add patch area [global[v]]patch area at vertex v
          
    for all vertices v
      vec_node[v] <- vec_node[v] / patch area[v]
      
    When this works, renormalize so that the cell and node vectors have the same L1 norm
*/      
#undef __FUNCT__
#define __FUNCT__ "CellToNodeInterpolation"
extern PetscErrorCode CellToNodeInterpolation(DM dm, Vec node_vec, Vec cell_vec, VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt		    dof;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscInt        k, j, i;
	PetscInt        c;
	PetscReal		hx,hy,hz;
	PetscReal		****node_arraydof;
	PetscReal		****cell_arraydof;
	PetscReal		***node_array;
	PetscReal		***cell_array;
	PetscReal		***volsum_array;
	PetscReal		****coords_array;
	PetscReal		***vol_array;
	Vec				volume;

	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(dm,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(node_vec,0.);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal, &volume);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) volume,"Volume");CHKERRQ(ierr);
	ierr = VecSet(volume,0.0);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,volume,&vol_array);CHKERRQ(ierr);
	if (dof == 1){
		ierr = DMDAVecGetArray(dm, node_vec,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecGetArray(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		ierr = DMDAVecGetArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
		ierr = DMDAVecGetArrayDOF(dm, cell_vec,&cell_arraydof);CHKERRQ(ierr);   
	}
	if (xs+xm == nx) xm--;
	if (ys+ym == ny) ym--;
	if (zs+zm == nz) zm--;	
	
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				for(k = 0; k < 2; k++){
					for(j = 0; j < 2; j++){
						for(i = 0; i < 2; i++){
							vol_array[ek+k][ej+j][ei+i] += hz*hy*hx;
							if(dof == 1)
							{
								node_array[ek+k][ej+j][ei+i] += cell_array[ek][ej][ei]*hx*hy*hz;
							}
							else
							{
								for(c = 0; c < dof; c++){
									node_arraydof[ek+k][ej+j][ei+i][c] += cell_arraydof[ek][ej][ei][c]*hx*hy*hz;
								}
							}
						}
					}
				}
			}
		}
	}
	if (xs+xm == nx-1) xm++;
	if (ys+ym == ny-1) ym++;
	if (zs+zm == nz-1) zm++;	
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
				if(dof == 1)
				{
					node_array[ek][ej][ei] = node_array[ek][ej][ei]/vol_array[ek][ej][ei];
				}
				else
				{
					for(c = 0; c < dof; c++){
					node_arraydof[ek][ej][ei][c] = cell_arraydof[ek][ej][ei][c]/vol_array[ek][ej][ei];
					}
				}
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal, volume,&vol_array);CHKERRQ(ierr);    
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