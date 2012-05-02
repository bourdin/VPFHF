/*
 VFPermfield.c
 Compute the total crack opening displacement (crack volume) by integrating u.\nabla v 
 over the entire domain
 
 (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
 */
#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"

#undef __FUNCT__
#define __FUNCT__ "VolumetricCrackOpening"
extern PetscErrorCode VolumetricCrackOpening(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal			****coords_array;
	PetscReal			****u_array;
	Vec					u_local;
	PetscReal			***v_array;
	Vec				    v_local;
    PetscReal			MyVol_change_total1 = 0;
	PetscReal			***volcrackopening_array;
	Vec					volcrackopening_local;   
	PetscReal			myCrackVolumeLocal = 0.,myCrackVolume = 0.;
	
	
  	PetscFunctionBegin;
    ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr); 
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr); 
	
	ierr = VecSet(fields->VolCrackOpening,0.);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&volcrackopening_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->VolCrackOpening,INSERT_VALUES,volcrackopening_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->VolCrackOpening,INSERT_VALUES,volcrackopening_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,volcrackopening_local,&volcrackopening_array);CHKERRQ(ierr); 
	*CrackVolume = 0.;
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
				ierr = VolumetricCrackOpening3D_local(&myCrackVolumeLocal, volcrackopening_array, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);				
                myCrackVolume += myCrackVolumeLocal; 
			}
		}
	}
	ierr = MPI_Allreduce(&myCrackVolume,CrackVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
/*	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n###################################################################\n");CHKERRQ(ierr);	
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n Crack volume calculated using gradv.u \t = %g\n", *CrackVolume);CHKERRQ(ierr);	
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
*/	
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr); 
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr); 	
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,volcrackopening_local,&volcrackopening_array);CHKERRQ(ierr); 	
	ierr = DMRestoreLocalVector(ctx->daScal,&volcrackopening_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScal,volcrackopening_local,ADD_VALUES,fields->VolCrackOpening);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScal,volcrackopening_local,ADD_VALUES,fields->VolCrackOpening);CHKERRQ(ierr);
	
	ierr = CellToNodeInterpolation(ctx->daScal, fields->FVCellndof, fields->VolCrackOpening, ctx); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumetricCrackOpening3D_local"
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ****displ_array, PetscReal ***vfield_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e)
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
	PetscReal		element_vol = 0;

	PetscFunctionBegin;
	ierr = PetscMalloc3(e->ng,PetscReal,&dx_vfield_loc,e->ng,PetscReal,&dy_vfield_loc,e->ng,PetscReal,&dz_vfield_loc);CHKERRQ(ierr);
	ierr = PetscMalloc3(e->ng,PetscReal,&udispl_loc,e->ng,PetscReal,&vdispl_loc,e->ng,PetscReal,&wdispl_loc);CHKERRQ(ierr);

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
	}
	volcrackopening_array[ek][ej][ei] = 0.;	
	*CrackVolume_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
		volcrackopening_array[ek][ej][ei] += (udispl_loc[eg]*dx_vfield_loc[eg] + vdispl_loc[eg]*dy_vfield_loc[eg] + wdispl_loc[eg]*dz_vfield_loc[eg])*e->weight[eg]; 
		element_vol += e->weight[eg];
		*CrackVolume_local += (udispl_loc[eg]*dx_vfield_loc[eg] + vdispl_loc[eg]*dy_vfield_loc[eg] + wdispl_loc[eg]*dz_vfield_loc[eg])*e->weight[eg]; 
	}
	volcrackopening_array[ek][ej][ei] = volcrackopening_array[ek][ej][ei]/element_vol;
	ierr = PetscFree3(dx_vfield_loc,dy_vfield_loc,dz_vfield_loc);CHKERRQ(ierr);
	ierr = PetscFree3(udispl_loc,udispl_loc,udispl_loc);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}     

#undef __FUNCT__
#define __FUNCT__ "CellToNodeInterpolation"
extern PetscErrorCode CellToNodeInterpolation(DM dm,Vec node_vec,Vec cell_vec,VFCtx *ctx)
{
	PetscErrorCode  ierr;
	PetscInt		dof;
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
	Vec				volume_local;
	PetscReal		nodal_sum_local = 0.;
	PetscReal		cell_sum_local = 0.;
	PetscReal		TotalNodeSum = 0.;
	PetscReal		TotalCellSum = 0.;
	Vec				node_local;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(dm,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(node_vec,0.);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(ctx->daScal, &volume);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) volume,"Volume");CHKERRQ(ierr);
	ierr = VecSet(volume,0.0);CHKERRQ(ierr);	
	ierr = DMGetLocalVector(ctx->daScal,&volume_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,volume,INSERT_VALUES,volume_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,volume,INSERT_VALUES,volume_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,volume_local,&vol_array);CHKERRQ(ierr);	
	
	ierr = DMGetLocalVector(dm,&node_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(dm,node_vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(dm,node_vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
	
	if (dof == 1){
		ierr = DMDAVecGetArray(dm, node_local,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecGetArray(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		ierr = DMDAVecGetArrayDOF(dm, node_local,&node_arraydof);CHKERRQ(ierr);    
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
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArray(ctx->daScal,volume_local,&vol_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&volume_local);CHKERRQ(ierr); 
	ierr = DMLocalToGlobalBegin(ctx->daScal,volume_local,ADD_VALUES,volume);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScal,volume_local,ADD_VALUES,volume);CHKERRQ(ierr);
	if(dof == 1){
		ierr = DMDAVecRestoreArray(dm, node_local,&node_array);CHKERRQ(ierr);  
	}		
	else
	{
		ierr = DMDAVecRestoreArrayDOF(dm, node_local,&node_arraydof);CHKERRQ(ierr); 
	}
	ierr = DMRestoreLocalVector(dm,&node_local);CHKERRQ(ierr); 
	ierr = DMLocalToGlobalBegin(dm,node_local,ADD_VALUES,node_vec);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(dm,node_local,ADD_VALUES,node_vec);CHKERRQ(ierr);
	
	if (dof == 1){
		ierr = DMDAVecGetArray(dm, node_vec,&node_array);CHKERRQ(ierr);    
	}		
	else
	{
		ierr = DMDAVecGetArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
	}
	ierr = DMDAVecGetArray(ctx->daScal,volume,&vol_array);CHKERRQ(ierr);
	if (xs+xm == nx-1) xm++;
	if (ys+ym == ny-1) ym++;
	if (zs+zm == nz-1) zm++;	
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
				if(dof == 1)
				{
					node_array[ek][ej][ei] = node_array[ek][ej][ei]/vol_array[ek][ej][ei];
					nodal_sum_local += PetscAbs(node_array[ek][ej][ei]);
					cell_sum_local += PetscAbs(cell_array[ek][ej][ei]);
				}
				else
				{
					for(c = 0; c < dof; c++){
						node_arraydof[ek][ej][ei][c] = node_arraydof[ek][ej][ei][c]/vol_array[ek][ej][ei];
						nodal_sum_local += PetscAbs(node_arraydof[ek][ej][ei][c]);
						cell_sum_local += PetscAbs(cell_arraydof[ek][ej][ei][c]);
					}
				}
			}
		}
	}	
	ierr = DMDAVecRestoreArray(ctx->daScal,volume,&vol_array);CHKERRQ(ierr);
	if(dof == 1){
		ierr = DMDAVecRestoreArray(dm, node_vec,&node_array);CHKERRQ(ierr);    
		ierr = DMDAVecRestoreArray(dm, cell_vec,&cell_array);CHKERRQ(ierr);    
	}		
	else
	{
		ierr = DMDAVecRestoreArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);    
		ierr = DMDAVecRestoreArrayDOF(dm, cell_vec,&cell_arraydof);CHKERRQ(ierr);    
	}
	ierr = MPI_Reduce(&nodal_sum_local,&TotalNodeSum,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Reduce(&cell_sum_local,&TotalCellSum,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\nNodal sum\t= %g\nCell sum\t= %g\n", TotalNodeSum, TotalCellSum);CHKERRQ(ierr);	
    PetscFunctionReturn(0);
}
