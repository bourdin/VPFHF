/*
 VFPermfield.c
 Compute the total crack opening displacement (crack volume) by integrating u.\nabla v
 over the entire domain
 
 (c) 2012	Chukwudozie, LSU
 */
#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFPermfield.h"

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
	PetscReal       ****coords_array;
	PetscReal       ****u_array;
	Vec             u_local;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscReal       myCrackVolumeLocal = 0.,myCrackVolume = 0.;
  Vec             CellVolCrackOpening;
  Vec             VolCrackOpening_local;
	PetscReal       ***volcrackopening_array;

	
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	
  
	ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daScalCell,&CellVolCrackOpening);CHKERRQ(ierr);
  ierr = VecSet(CellVolCrackOpening,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScalCell,&VolCrackOpening_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,CellVolCrackOpening,INSERT_VALUES,VolCrackOpening_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,CellVolCrackOpening,INSERT_VALUES,VolCrackOpening_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,VolCrackOpening_local,&volcrackopening_array);CHKERRQ(ierr);
  
	*CrackVolume = 0.;
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
  /*
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\n###################################################################\n");CHKERRQ(ierr);
   ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n Crack volume calculated using gradv.u \t = %g\n", *CrackVolume);CHKERRQ(ierr);
   ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
   */
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScalCell,VolCrackOpening_local,&volcrackopening_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&VolCrackOpening_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,VolCrackOpening_local,ADD_VALUES,CellVolCrackOpening);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,VolCrackOpening_local,ADD_VALUES,CellVolCrackOpening);CHKERRQ(ierr);
  
  ierr = VecSet(fields->VolCrackOpening,0.);CHKERRQ(ierr);
  ierr = VecCopy(CellVolCrackOpening,fields->pmult);CHKERRQ(ierr);
  
	ierr = CellToNodeInterpolation(fields->VolCrackOpening,CellVolCrackOpening,ctx); CHKERRQ(ierr);
  ierr = VecDestroy(&CellVolCrackOpening);CHKERRQ(ierr);

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
  PetscReal sum = 0, sum1 = 0, sum2 = 0;

	volcrackopening_array[ek][ej][ei] = 0.;
	*CrackVolume_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
    sum2 += wdispl_loc[eg]*dz_vfield_loc[eg]*e->weight[eg];
    sum1 += vdispl_loc[eg]*dy_vfield_loc[eg]*e->weight[eg];
    sum += udispl_loc[eg]*dx_vfield_loc[eg]*e->weight[eg];
		volcrackopening_array[ek][ej][ei] += (udispl_loc[eg]*dx_vfield_loc[eg] + vdispl_loc[eg]*dy_vfield_loc[eg] + wdispl_loc[eg]*dz_vfield_loc[eg])*e->weight[eg];
		element_vol += e->weight[eg];
		*CrackVolume_local += (udispl_loc[eg]*dx_vfield_loc[eg] + vdispl_loc[eg]*dy_vfield_loc[eg] + wdispl_loc[eg]*dz_vfield_loc[eg])*e->weight[eg];
	}  
  volcrackopening_array[ek][ej][ei] = volcrackopening_array[ek][ej][ei]/element_vol;
  if(volcrackopening_array[ek][ej][ei] < 0){
    volcrackopening_array[ek][ej][ei] = 0.;
  }
	ierr = PetscFree3(dx_vfield_loc,dy_vfield_loc,dz_vfield_loc);CHKERRQ(ierr);
	ierr = PetscFree3(udispl_loc,vdispl_loc,wdispl_loc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CellToNodeInterpolation"
extern PetscErrorCode CellToNodeInterpolation(Vec node_vec,Vec cell_vec,VFCtx *ctx)
{
	PetscErrorCode  ierr;
	PetscInt        dof;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscInt        k, j, i;
	PetscInt        c;
	PetscReal       hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ***vol_array;
	Vec             volume;
	Vec             volume_local;
	PetscReal       nodal_sum_local = 0.;
	PetscReal       cell_sum_local = 0.;
	PetscReal       TotalNodeSum = 0.;
	PetscReal       TotalCellSum = 0.;
	Vec             node_local;
  PetscReal       ***node_array;
	PetscReal       ****node_arraydof;
	PetscReal       ***cell_array;
	PetscReal       ****cell_arraydof;
	DM              dm;
	
	PetscFunctionBegin;
	ierr = PetscObjectQuery((PetscObject)node_vec,"DM",(PetscObject*)&dm);CHKERRQ(ierr);
	if (!dm) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
	ierr = DMDAGetInfo(dm,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal, &volume);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) volume,"Volume");CHKERRQ(ierr);
	ierr = VecSet(volume,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&volume_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,volume,INSERT_VALUES,volume_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,volume,INSERT_VALUES,volume_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,volume_local,&vol_array);CHKERRQ(ierr);
	if (dof == 1){
		ierr = DMGetLocalVector(dm,&node_local);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(dm,node_vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(dm,node_vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(dm,node_local,&node_array);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(ctx->daScalCell,cell_vec,&cell_array);CHKERRQ(ierr);
	}
	else
	{
    /*For now, no cell DA for variables with dof = 3*/
		ierr = DMGetLocalVector(dm,&node_local);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(dm,node_vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(dm,node_vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
		ierr = DMDAVecGetArrayDOF(dm, node_local,&node_arraydof);CHKERRQ(ierr);
		ierr = DMDAVecGetArrayDOF(ctx->daVFperm,cell_vec,&cell_arraydof);CHKERRQ(ierr);
	}
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
    ierr = DMRestoreLocalVector(dm,&node_local);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(dm,node_local,ADD_VALUES,node_vec);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(dm,node_local,ADD_VALUES,node_vec);CHKERRQ(ierr);
	}
	else
	{
		ierr = DMDAVecRestoreArrayDOF(dm, node_local,&node_arraydof);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm,&node_local);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(dm,node_local,ADD_VALUES,node_vec);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(dm,node_local,ADD_VALUES,node_vec);CHKERRQ(ierr);
 }
	if (dof == 1){
		ierr = DMDAVecGetArray(dm, node_vec,&node_array);CHKERRQ(ierr);
	}
	else
	{
		ierr = DMDAVecGetArrayDOF(dm, node_vec,&node_arraydof);CHKERRQ(ierr);
	}
	ierr = DMDAVecGetArray(ctx->daScal,volume,&vol_array);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
				if(dof == 1)
				{
					node_array[ek][ej][ei] = node_array[ek][ej][ei]/vol_array[ek][ej][ei];
					nodal_sum_local += PetscAbs(node_array[ek][ej][ei]);
				}
				else
				{
					for(c = 0; c < dof; c++){
						node_arraydof[ek][ej][ei][c] = node_arraydof[ek][ej][ei][c]/vol_array[ek][ej][ei];
						nodal_sum_local += PetscAbs(node_arraydof[ek][ej][ei][c]);
					}
          
				}
			}
		}
	}
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
				if(dof == 1)
				{
            cell_sum_local += PetscAbs(cell_array[ek][ej][ei]);
				}
				else
				{
					for(c = 0; c < dof; c++){
						cell_sum_local += PetscAbs(cell_arraydof[ek][ej][ei][c]);
					}          
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(ctx->daScal,volume,&vol_array);CHKERRQ(ierr);
	if(dof == 1){
		ierr = DMDAVecRestoreArray(dm,node_vec,&node_array);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(ctx->daScalCell,cell_vec,&cell_array);CHKERRQ(ierr);
	}
  else
	{
		ierr = DMDAVecRestoreArrayDOF(dm,node_vec,&node_arraydof);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,cell_vec,&cell_arraydof);CHKERRQ(ierr);
	}
  ierr = VecDestroy(&volume);CHKERRQ(ierr);
	ierr = MPI_Reduce(&nodal_sum_local,&TotalNodeSum,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Reduce(&cell_sum_local,&TotalCellSum,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumetricLeakOffRate"
extern PetscErrorCode VolumetricLeakOffRate(PetscReal *LeakOffRate, VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ****q_array;
	Vec             q_local;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscReal       ***volleakoffrate_array;
  Vec             CellVolLeakoffRate;
	Vec             VolLeakoffRate_local;
	PetscReal       myLeakOffRateLocal = 0.,myLeakOffRate = 0.;
	
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVect,&q_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->velocity,INSERT_VALUES,q_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->velocity,INSERT_VALUES,q_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,q_local,&q_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daScalCell,&CellVolLeakoffRate);CHKERRQ(ierr);
  ierr = VecSet(CellVolLeakoffRate,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScalCell,&VolLeakoffRate_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,CellVolLeakoffRate,INSERT_VALUES,VolLeakoffRate_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,CellVolLeakoffRate,INSERT_VALUES,VolLeakoffRate_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,VolLeakoffRate_local,&volleakoffrate_array);CHKERRQ(ierr);
  
	*LeakOffRate = 0.;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricCrackOpening3D_local(&myLeakOffRateLocal, volleakoffrate_array, q_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        myLeakOffRate += myLeakOffRateLocal;
			}
		}
	}
	ierr = MPI_Allreduce(&myLeakOffRate,LeakOffRate,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,q_local,&q_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&q_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  
	ierr = DMDAVecRestoreArray(ctx->daScalCell,VolLeakoffRate_local,&volleakoffrate_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&VolLeakoffRate_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,VolLeakoffRate_local,ADD_VALUES,CellVolLeakoffRate);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,VolLeakoffRate_local,ADD_VALUES,CellVolLeakoffRate);CHKERRQ(ierr);
  
  ierr = VecSet(fields->VolLeakOffRate,0.);CHKERRQ(ierr);
	ierr = CellToNodeInterpolation(fields->VolLeakOffRate,CellVolLeakoffRate,ctx); CHKERRQ(ierr);
  ierr = VecDestroy(&CellVolLeakoffRate);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumetricLeakOffRate_local"
extern PetscErrorCode VolumetricLeakOffRate_local(PetscReal *LeakoffRate_local, PetscReal ***volleakoffrate_array, PetscReal ****q_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k;
	PetscInt		eg;
	PetscReal		*dx_vfield_loc;
	PetscReal		*dy_vfield_loc;
	PetscReal		*dz_vfield_loc;
	PetscReal		*velx_loc;
	PetscReal		*vely_loc;
	PetscReal		*velz_loc;
	PetscReal		element_vol = 0;
	
	PetscFunctionBegin;
	ierr = PetscMalloc3(e->ng,PetscReal,&dx_vfield_loc,e->ng,PetscReal,&dy_vfield_loc,e->ng,PetscReal,&dz_vfield_loc);CHKERRQ(ierr);
	ierr = PetscMalloc3(e->ng,PetscReal,&velx_loc,e->ng,PetscReal,&vely_loc,e->ng,PetscReal,&velz_loc);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
		velx_loc[eg] = 0.;
		vely_loc[eg] = 0.;
		velz_loc[eg] = 0.;
		dx_vfield_loc[eg] = 0.;
		dy_vfield_loc[eg] = 0.;
		dz_vfield_loc[eg] = 0.;
	}
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
					velx_loc[eg] += q_array[ek+k][ej+j][ei+i][0] * e->phi[k][j][i][eg];
					vely_loc[eg] += q_array[ek+k][ej+j][ei+i][1] * e->phi[k][j][i][eg];
					velz_loc[eg] += q_array[ek+k][ej+j][ei+i][2] * e->phi[k][j][i][eg];
					dx_vfield_loc[eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][0][eg];
					dy_vfield_loc[eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][1][eg];
					dz_vfield_loc[eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][2][eg];
				}
			}
		}
	}
	volleakoffrate_array[ek][ej][ei] = 0.;
	*LeakoffRate_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
		volleakoffrate_array[ek][ej][ei] += (velx_loc[eg]*dx_vfield_loc[eg] + vely_loc[eg]*dy_vfield_loc[eg] + velz_loc[eg]*dz_vfield_loc[eg])*e->weight[eg];
		element_vol += e->weight[eg];
		*LeakoffRate_local += (velx_loc[eg]*dx_vfield_loc[eg] + vely_loc[eg]*dy_vfield_loc[eg] + velz_loc[eg]*dz_vfield_loc[eg])*e->weight[eg];
	}
	volleakoffrate_array[ek][ej][ei] = volleakoffrate_array[ek][ej][ei]/element_vol;
	ierr = PetscFree3(dx_vfield_loc,dy_vfield_loc,dz_vfield_loc);CHKERRQ(ierr);
	ierr = PetscFree3(velx_loc,vely_loc,velz_loc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Permeabilityfield"
extern PetscErrorCode Permeabilityfield(PetscReal *COD_local, PetscReal ***volcrackopening_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element2D *e, FACE face)
{
	PetscErrorCode ierr;
	PetscInt       i,j,k,g;
	PetscReal		*volopening_loc;
	
	PetscFunctionBegin;
	ierr = PetscMalloc(e->ng*sizeof(PetscReal),&volopening_loc);CHKERRQ(ierr);
	for (g = 0; g < e->ng; g++){
		volopening_loc[g] = 0.;
	}
	switch (face) {
		case X0:
			for (k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							volopening_loc[g] += e->phi[k][j][i][g] * volcrackopening_array[ek+k][ej+j][ei+i];
						}
					}
				}
			}
			break;
		case Y0:
			for (k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							volopening_loc[g] += e->phi[k][j][i][g] * volcrackopening_array[ek+k][ej+j][ei+i];
						}
					}
				}
			}
			break;
		case Z0:
			for (k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							volopening_loc[g] += e->phi[k][j][i][g] * volcrackopening_array[ek+k][ej+j][ei+i];
						}
					}
				}
			}
			break;
		case Z1:
			break;
		case Y1:
			break;
		case X1:
			break;
		default:
			/*	SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be one of X0,X1,Y0,Y1,Z0,Z1\n");	*/
			break;
	}
	*COD_local = 0.;
	for(g = 0; g < e->ng; g++){
		*COD_local += e->weight[g] * volopening_loc[g];
	}
	ierr = PetscFree(volopening_loc);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PermeabilityUpDate"
extern PetscErrorCode PermeabilityUpDate(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscReal       hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ***cod_array;
	Vec             cod_local;
	PetscReal       ****perm_array;
	Vec             perm_local;
	PetscReal       ****u_array;
	Vec             u_local;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscReal       ***pmult_array;
	Vec             pmult_local;
	Vec             COD;
	PetscReal       myCOD = 0.;
	PetscInt        k, j, i,c;
	PetscInt        num_int_cell;
	PetscReal       Ele_v_ave = 0.;
	
	PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
	ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	
	ierr = VecSet(fields->vfperm,0.);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);
  
	ierr = DMCreateGlobalVector(ctx->daScalCell,&COD);CHKERRQ(ierr);
	ierr = VecSet(COD,0.);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScalCell,&cod_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,COD,INSERT_VALUES,cod_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,COD,INSERT_VALUES,cod_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,cod_local,&cod_array);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricCrackOpening3D_local(&myCOD, cod_array, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
			}
		}
	}
	PetscReal	vmult = 0;
	PetscReal  maxperm = 0;
	num_int_cell  = 2;
	ierr = PetscOptionsInt("-int_cells","\n\tNumber of cells for integration","",num_int_cell,&num_int_cell,PETSC_NULL);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				Ele_v_ave = 0.;
				for (k = 0; k < 2; k++) {
					for (j = 0; j < 2; j++) {
						for (i = 0; i < 2; i++) {
							Ele_v_ave += v_array[ek+k][ej+j][ei+i];
						}
					}
				}
				Ele_v_ave = Ele_v_ave/8.;
        
				if(Ele_v_ave > ctx->vfprop.irrevtol){
					vmult = 1.;
					vmult = Ele_v_ave;
				}
				else{
            //					vmult = 0.;
					vmult = Ele_v_ave;
				}
				for(c = (-1*num_int_cell); c < num_int_cell+1; c++){
					if(ei+c >= 0 && ei+c <= nx-1){
						hx = coords_array[ek][ej][ei+c+1][0]-coords_array[ek][ej][ei+c][0];
						perm_array[ek][ej][ei][1] += cod_array[ek][ej][ei+c]*hx*(1.-vmult*vmult);
					}
					if(ej+c >= 0 && ej+c <= ny-1){
						hy = coords_array[ek][ej+c+1][ei][1]-coords_array[ek][ej+c][ei][1];
						perm_array[ek][ej][ei][2] += cod_array[ek][ej+c][ei]*hy*(1.-vmult*vmult);
					}
					if(ek+c >= 0 && ek+c <= nz-1){
						hz = coords_array[ek+c+1][ej][ei][2]-coords_array[ek+c][ej][ei][2];
						perm_array[ek][ej][ei][0] += cod_array[ek+c][ej][ei]*hz*(1.-vmult*vmult);
						perm_array[ek][ej][ei][0] = cod_array[ek][ej][ei];
					}
				}				
				maxperm = perm_array[ek][ej][ei][0];
				for(c = 1; c < 2; c++){
					if (perm_array[ek][ej][ei][c] < maxperm) {
						maxperm = perm_array[ek][ej][ei][c];
					}
					perm_array[ek][ej][ei][5] = maxperm;
					maxperm = 0.;
				}
        
			}
		}
	}
	PetscReal sum0 = 0, sum1 = 0, sum2 = 0;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				sum0 = sum0 + perm_array[ek][ej][ei][0]*(hy*hz);
				sum1 = sum1 + perm_array[ek][ej][ei][1]*(hx*hz);
				sum2 = sum2 + perm_array[ek][ej][ei][2]*(hx*hy);
			}
		}
	}
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daVFperm,perm_local,INSERT_VALUES,fields->vfperm);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daVFperm,perm_local,INSERT_VALUES,fields->vfperm);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecDestroy(&COD);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SimplePermeabilityUpDate"
extern PetscErrorCode SimplePermeabilityUpDate(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscReal       ****coords_array;	
	PetscReal       ****perm_array;
	Vec             perm_local;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscInt        k, j, i,c;
	PetscReal       Ele_v_ave = 0.;
	
	PetscFunctionBegin;
    ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	
	ierr = VecSet(fields->vfperm,0.);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);
  
	PetscReal  maxperm = ctx->vfprop.permmax;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				Ele_v_ave = 0.;
				for (k = 0; k < 2; k++) {
					for (j = 0; j < 2; j++) {
						for (i = 0; i < 2; i++) {
							Ele_v_ave += v_array[ek+k][ej+j][ei+i];
						}
					}
				}
				Ele_v_ave = Ele_v_ave/8.;
				if(Ele_v_ave < 0.) {
				  for (c = 0; c < 3; c++) {
				    perm_array[ek][ej][ei][c] = maxperm;
				  }
				}
				else {
				  for (c = 0; c < 3; c++) {
				    perm_array[ek][ej][ei][c] = maxperm*(1-Ele_v_ave);
                  }
                }				  
			}
		}
	}

	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daVFperm,perm_local,INSERT_VALUES,fields->vfperm);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daVFperm,perm_local,INSERT_VALUES,fields->vfperm);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TrialFunctionCompute1"
extern PetscErrorCode TrialFunctionCompute1(PetscReal *FunctionValue, VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei, c, c1;
	PetscReal       hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ****u_array;
	Vec             u_local;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscReal       ***pmult_array;
	Vec             pmult_local;
  PetscReal       sum = 0;
  Vec             COD;
	PetscReal       ***cod_array;
	Vec             cod_local;
	PetscReal       myCOD = 0.;
	PetscInt        num_int_cell;

	PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
	ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daScalCell,&COD);CHKERRQ(ierr);
	ierr = VecSet(COD,0.);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScalCell,&cod_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,COD,INSERT_VALUES,cod_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,COD,INSERT_VALUES,cod_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,cod_local,&cod_array);CHKERRQ(ierr);

  for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricCrackOpening3D_local(&myCOD, cod_array, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
			}
		}
	}
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScalCell,cod_local,&cod_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&cod_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,cod_local,ADD_VALUES,COD);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,cod_local,ADD_VALUES,COD);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daScalCell,&cod_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,COD,INSERT_VALUES,cod_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,COD,INSERT_VALUES,cod_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,cod_local,&cod_array);CHKERRQ(ierr);
  
  num_int_cell  = 2;
	ierr = PetscOptionsInt("-int_cells","\n\tNumber of cells for integration","",num_int_cell,&num_int_cell,PETSC_NULL);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {

				for(c = (-1*num_int_cell); c < num_int_cell+1; c++){
          for(c1 = (-1*num_int_cell); c1 < num_int_cell+1; c1++){

            if( (ei+c >= 0 && ei+c <= nx-1) && (ej+c1 >= 0 && ej+c1 <= ny-1)){
              hx = coords_array[ek][ej+c1][ei+c+1][0]-coords_array[ek][ej+c1][ei+c][0];
              hy = coords_array[ek][ej+c1+1][ei+c][1]-coords_array[ek][ej+c1][ei+c][1];
              pmult_array[ek][ej][ei] = hx*hy*pow(cod_array[ek][ej+c1][ei+c],3);
            }
					}
				}
//        ierr = PetscPrintf(PETSC_COMM_WORLD," ek = %d, ej = %d, ei %= %d, pmult = %e \t cod = %e\n", ek,ej,ei,pmult_array[ek][ej][ei],cod_array[ek][ej][ei]);CHKERRQ(ierr);

			}
		}
	}

  ierr = DMDAVecRestoreArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);

  
  for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
        hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        sum += pmult_array[ek][ej][ei]*hx*hy*hz;
			}
		}
	}  
//  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n Volume integral computed from sum  = %g\n", sum);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScalCell,cod_local,&cod_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&cod_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,cod_local,ADD_VALUES,COD);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,cod_local,ADD_VALUES,COD);CHKERRQ(ierr);
  ierr = VecDestroy(&COD);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TrialFunctionCompute"
extern PetscErrorCode TrialFunctionCompute(PetscReal *FunctionValue, VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscReal       hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ****u_array;
	Vec             u_local;
	PetscReal       ***v_array;
	Vec             v_local;
  PetscReal       myFunctionLocal = 0.,myFunction = 0.;

	PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
	ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricFunction_local(&myFunctionLocal, PETSC_NULL, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        myFunction += myFunctionLocal;    
			}
		}
	}
	ierr = MPI_Allreduce(&myFunction,FunctionValue,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumetricFunction_local"
extern PetscErrorCode VolumetricFunction_local(PetscReal *Function_local, PetscReal ***value_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e)
{
	PetscErrorCode ierr;
	PetscInt        i, j, k,c;
	PetscInt        eg;
  PetscReal       *dv_elem[3],*u_elem[3],*n_elem[3],*n_mag_elem;

	PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&n_mag_elem);CHKERRQ(ierr);

  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
      n_elem[c][eg] = 0;
    }
    n_mag_elem[eg] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
      n_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[0][eg] = dv_elem[0][eg]/n_mag_elem[eg];
      n_elem[1][eg] = dv_elem[1][eg]/n_mag_elem[eg];
      n_elem[2][eg] = dv_elem[2][eg]/n_mag_elem[eg];
    }
    if((ierr = PetscIsInfOrNanScalar(n_elem[0][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[1][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[2][eg])) )
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = n_mag_elem[eg] = 0;
    }
  }
  *Function_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
      *Function_local += 4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*n_mag_elem[eg]*e->weight[eg];  
  }
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],n_mag_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}










#undef __FUNCT__
#define __FUNCT__ "VF_PermeabilityUpDate"
extern PetscErrorCode VF_PermeabilityUpDate(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei;
	PetscReal       hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ***cod_array;
	Vec             cod_local;
	PetscReal       ****perm_array;
	Vec             perm_local;
	PetscReal       ****u_array;
	Vec             u_local;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscReal       ***pmult_array;
	Vec             pmult_local;
	Vec             COD;
	PetscReal       myCOD = 0.;
	PetscInt        k, j, i,c;
	PetscInt        num_int_cell;
	PetscReal       Ele_v_ave = 0.;
	
	PetscFunctionBegin;
    ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
	ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	
	ierr = VecSet(fields->vfperm,0.);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);
  
	ierr = DMCreateGlobalVector(ctx->daScalCell,&COD);CHKERRQ(ierr);
	ierr = VecSet(COD,0.);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScalCell,&cod_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,COD,INSERT_VALUES,cod_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,COD,INSERT_VALUES,cod_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,cod_local,&cod_array);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricCrackOpening3D_local(&myCOD, cod_array, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
			}
		}
	}
	PetscReal	vmult = 0;
	PetscReal  maxperm = 10;
	PetscReal  minperm = 1.e-1;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				Ele_v_ave = 0.;
				for (k = 0; k < 2; k++) {
					for (j = 0; j < 2; j++) {
						for (i = 0; i < 2; i++) {
							Ele_v_ave += v_array[ek+k][ej+j][ei+i];
						}
					}
				}
				Ele_v_ave = Ele_v_ave/8.;
        
				if(Ele_v_ave > 0.2){
                    for(c = 0; c < 2; c++)
                      perm_array[ek][ej][ei][c] = minperm;
				}
				else{
                  for(c = 0; c < 2; c++)
                    perm_array[ek][ej][ei][c] = maxperm;
				}
			}
		}
	}

	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daVFperm,perm_local,INSERT_VALUES,fields->vfperm);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daVFperm,perm_local,INSERT_VALUES,fields->vfperm);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecDestroy(&COD);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

