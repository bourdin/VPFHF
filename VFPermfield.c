/*
 VFPermfield.c
 Compute the total crack opening displacement (crack volume) by integrating u.\nabla v
 over the entire domain
 
 (c) 2012 - 2014 	Chukwudozie, LSU, B. Bourdin, LSU
 */
#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFPermfield.h"

#undef __FUNCT__
#define __FUNCT__ "VolumetricCrackOpening3D_local"
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ****displ_array, PetscReal ***vfield_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
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
	ierr = PetscMalloc3(e->ng,&dx_vfield_loc,e->ng,&dy_vfield_loc,e->ng,&dz_vfield_loc);CHKERRQ(ierr);
	ierr = PetscMalloc3(e->ng,&udispl_loc,e->ng,&vdispl_loc,e->ng,&wdispl_loc);CHKERRQ(ierr);
  
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
  
	ierr = PetscFree3(dx_vfield_loc,dy_vfield_loc,dz_vfield_loc);CHKERRQ(ierr);
	ierr = PetscFree3(udispl_loc,vdispl_loc,wdispl_loc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

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
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricCrackOpening3D_local(&myCrackVolumeLocal, volcrackopening_array, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        myCrackVolume += myCrackVolumeLocal;
			}
		}
	}
	ierr = MPI_Allreduce(&myCrackVolume,CrackVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
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
  
	ierr = CellToNodeInterpolation(fields->VolCrackOpening,CellVolCrackOpening,ctx); CHKERRQ(ierr);
  ierr = VecDestroy(&CellVolCrackOpening);CHKERRQ(ierr);
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_IntegrateDisplacement_local"
extern PetscErrorCode VF_IntegrateDisplacement_local(Vec Uv,Vec Uc,VFCtx *ctx)
{
	PetscErrorCode  ierr;
	PetscInt        xs,xm,nx;
	PetscInt        ys,ym,ny;
	PetscInt        zs,zm,nz;
	PetscInt        ek, ej, ei, dof;
	PetscInt        k, j, i, n;
	PetscInt        c;
  Vec             uv_local;
  Vec             uc_local;
  PetscReal       ****uv_array;
	PetscReal       ****uc_array;
	DM              dm;
	
	PetscFunctionBegin;
  ierr = VecGetDM(Uv,&dm);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dm,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(dm,&uv_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,Uv,INSERT_VALUES,uv_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,Uv,INSERT_VALUES,uv_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(dm, uv_local,&uv_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(dm,&uc_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,Uc,INSERT_VALUES,uc_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,Uc,INSERT_VALUES,uc_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(dm, uc_local,&uc_array);CHKERRQ(ierr);
  n = 15;
  
  for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys + ym; ej++) {
			for (ei = xs; ei < xs + xm; ei++) {
        for(c = 0; c < 3; c++){
          uv_array[ek][ej][ei][c] = 0;
          uv_array[ek][ej][ei][c] = uc_array[ek][ej][ei][c];
        }
        for(k = 1; k < n; k++){
          if(ek-k >= 0 && ek+k <= nz-1){
            uv_array[ek][ej][ei][2] += (uc_array[ek+k][ej][ei][2]+uc_array[ek-k][ej][ei][2]);
          }
        }
        for(j = 1; j < n; j++){
          if(ej-j >= 0 && ej+j <= ny-1){
            uv_array[ek][ej][ei][1] += (uc_array[ek][ej+j][ei][1]+uc_array[ek][ej-j][ei][1]);
          }
        }
        for(i = 1; i < n; i++){
          if(ei-i >= 0 && ei+i <= nx-1){
            uv_array[ek][ej][ei][0] += (uc_array[ek][ej][ei+i][0]+uc_array[ek][ej][ei-i][0]);
          }
        }
      }
    }
  }
  
  ierr = DMDAVecRestoreArrayDOF(dm, uv_local,&uv_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&uv_local);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,uv_local,INSERT_VALUES,Uv);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,uv_local,INSERT_VALUES,Uv);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(dm, uc_local,&uc_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&uc_local);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,uc_local,INSERT_VALUES,Uc);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,uc_local,INSERT_VALUES,Uc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeUcdotGradVlocal"
extern PetscErrorCode ComputeUcdotGradVlocal(PetscReal *cod, PetscReal *grad_elem, PetscReal *n_elem, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFEElement3D *s)
{
	PetscInt        k, j, i;
	PetscInt        c;
	PetscReal       dv_elem[3],u_elem[3];
  PetscReal       dv_mag_elem = 0;
  
	PetscFunctionBegin;
  for (c = 0; c < 3; c++){
    dv_elem[c] = 0;
    u_elem[c] = 0;
  }
  for (k = 0; k < s->nphiz; k++) {
    for (j = 0; j < s->nphiy; j++) {
      for (i = 0; i < s->nphix; i++) {
        for (c = 0; c < 3; c++){
          dv_elem[c] += v_array[ek+k][ej+j][ei+i] * s->dphi[k][j][i][c];
          u_elem[c] += u_array[ek+k][ej+j][ei+i][c]*s->phi[k][j][i];
        }
      }
    }
  }
  dv_mag_elem = sqrt((pow(dv_elem[0],2))+(pow(dv_elem[1],2))+(pow(dv_elem[2],2)));
  for(c = 0; c < 3; c++){
    n_elem[c] = dv_elem[c]/dv_mag_elem;
  }
  if((PetscIsInfOrNanScalar(n_elem[0])) || (PetscIsInfOrNanScalar(n_elem[1])) || (PetscIsInfOrNanScalar(n_elem[2])))
  {
    n_elem[0] = n_elem[1] = n_elem[2] = dv_mag_elem = 0;
  }
  *cod = 0;
  for (c = 0; c < 3; c++){
    *cod += dv_elem[c]*u_elem[c];
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateFractureWidth"
extern PetscErrorCode UpdateFractureWidth(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscInt		    xs1,xm1;
	PetscInt		    ys1,ym1;
	PetscInt		    zs1,zm1;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ****u_array;
	Vec             u_local;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscReal       grad_cc[3] = {0,0,0};
	PetscReal       coorda_array[3]= {0,0,0};
	PetscReal       coordb_array[3]= {0,0,0};
	PetscReal       coordc_array[3]= {0,0,0};
  Vec             w_local;
	PetscReal       ***w_array;
  PetscReal       len;
  PetscReal       cod[2] = {0,0}, lc, cod_in;
  PetscReal		    hwx,hwy,hwz;
  PetscInt		    ekk, ejj, eii;
	PetscInt		    ek1, ej1, ei1,c;
  PetscReal       ave_V = 0;
  Vec             pmult_local;
	PetscReal       ***pmult_array;
  Vec             coordinates;
  Vec             coords_local;
  PetscReal       coords1[3]= {0,0,0};
	PetscReal       coords2[3]= {0,0,0};
	PetscReal       tlent1;
	PetscReal       tlent2;
  PetscReal       n_cc[3] = {0,0,0};
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daWScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daWScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetGhostCorners(ctx->daWScalCell,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daWScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daWScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daWScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daWScal,v_local,&v_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daWVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daWVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daWVect,fields->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daWVect,u_local,&u_array);CHKERRQ(ierr);
  
  ierr = DMGetCoordinates(ctx->daScal,&coordinates);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daWVect,&coords_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daWVect,coordinates,INSERT_VALUES,coords_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daWVect,coordinates,INSERT_VALUES,coords_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daWVect,coords_local,&coords_array);CHKERRQ(ierr);
  
  ierr = VecSet(fields->widthc,0.);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScalCell,&w_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,fields->widthc,INSERT_VALUES,w_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,fields->widthc,INSERT_VALUES,w_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,w_local,&w_array);CHKERRQ(ierr);
  
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
        hwx = hx/2.;
        hwy = hy/2.;
        hwz = hz/2.;
        ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        ierr = CartFEElement3DInit(&ctx->s3D,hwx,hwy,hwz,hx,hy,hz);CHKERRQ(ierr);
        ierr = ComputeAverageVlocal(&ave_V, v_array, ek, ej, ei, &ctx->s3D);
        ierr = ComputeCellCenterGradV_local(grad_cc, v_array, ek, ej, ei, &ctx->s3D);
        w_array[ek][ej][ei] = 0;
        if(ave_V < 1.0 && ave_V > 0.0){
          len = ctx->IntLengthRes;
          ierr = ComputeUcdotGradVlocal(&cod[0], grad_cc, n_cc, u_array, v_array, ek, ej, ei, &ctx->s3D);
          cod_in = cod[0];
          coorda_array[0] = coordb_array[0] = coordc_array[0] = (coords_array[ek][ej][ei+1][0]+coords_array[ek][ej][ei][0])/2.;
          coorda_array[1] = coordb_array[1] = coordc_array[1] = (coords_array[ek][ej+1][ei][1]+coords_array[ek][ej][ei][1])/2.;
          coorda_array[2] = coordb_array[2] = coordc_array[2] = (coords_array[ek+1][ej][ei][2]+coords_array[ek][ej][ei][2])/2.;
          for(c = 0; c < 3; c++){
            coordc_array[c] = coordb_array[c]+grad_cc[c]*len;
          }
          n_cc[2] = n_cc[1] = n_cc[0] = 0;
          lc = sqrt((pow(coordc_array[0]-coorda_array[0],2))+(pow(coordc_array[1]-coorda_array[1],2))+(pow(coordc_array[2]-coorda_array[2],2)));
          while(lc < ctx->WidthIntLength && ave_V < 1.0 && (grad_cc[0]*n_cc[0]+grad_cc[1]*n_cc[1]+grad_cc[2]*n_cc[2] >= 0.)){
            for (ek1 = zs1; ek1 < zs1+zm1; ek1++) {
              for (ej1 = ys1; ej1 < ys1+ym1; ej1++) {
                for (ei1 = xs1; ei1 < xs1+xm1; ei1++) {
                  if(
                     ((coords_array[ek1][ej1][ei1+1][0] >= coordc_array[0]) && (coords_array[ek1][ej1][ei1][0] <= coordc_array[0] ))    &&
                     ((coords_array[ek1][ej1+1][ei1][1] >= coordc_array[1]) && (coords_array[ek1][ej1][ei1][1] <= coordc_array[1] ))    &&
                     ((coords_array[ek1+1][ej1][ei1][2] >= coordc_array[2]) && (coords_array[ek1][ej1][ei1][2] <= coordc_array[2] ))
                     )
                  {
                    ekk = ek1;ejj = ej1;eii = ei1;
                  }
                }
              }
            }
            hx = coords_array[ekk][ejj][eii+1][0]-coords_array[ekk][ejj][eii][0];
            hy = coords_array[ekk][ejj+1][eii][1]-coords_array[ekk][ejj][eii][1];
            hz = coords_array[ekk+1][ejj][eii][2]-coords_array[ekk][ejj][eii][2];
            hwx = coordc_array[0]-coords_array[ekk][ejj][eii][0];
            hwy = coordc_array[1]-coords_array[ekk][ejj][eii][1];
            hwz = coordc_array[2]-coords_array[ekk][ejj][eii][2];
            ierr = CartFEElement3DInit(&ctx->s3D,hwx,hwy,hwz,hx,hy,hz);CHKERRQ(ierr);
            ierr = ComputeAverageVlocal(&ave_V, v_array, ekk, ejj, eii, &ctx->s3D);
            ierr = ComputeUcdotGradVlocal(&cod[1], grad_cc, n_cc, u_array, v_array, ekk, ejj, eii, &ctx->s3D);
            ierr = VFCartFEElement1DInit(&ctx->e1D,len);CHKERRQ(ierr);
            ierr = IntegrateUcdotGradVlocal(&w_array[ek][ej][ei],cod, &ctx->e1D);
            if((n_cc[0] == 0) && (n_cc[1] == 0) && (n_cc[2] == 0)){
              n_cc[0] = grad_cc[0];
              n_cc[1] = grad_cc[1];
              n_cc[2] = grad_cc[2];
            }
            for(c = 0; c < 3; c++){
              coords1[c] = coordc_array[c]+n_cc[c]*len;
              coords2[c] = coordc_array[c]-n_cc[c]*len;
            }
            tlent1 = sqrt((pow(coords1[0]-coorda_array[0],2))+(pow(coords1[1]-coorda_array[1],2))+(pow(coords1[2]-coorda_array[2],2)));
            tlent2 = sqrt((pow(coords2[0]-coorda_array[0],2))+(pow(coords2[1]-coorda_array[1],2))+(pow(coords2[2]-coorda_array[2],2)));
            
            if(tlent1 > tlent2){
              lc += sqrt((pow(coordc_array[0]-coords1[0],2))+(pow(coordc_array[1]-coords1[1],2))+(pow(coordc_array[2]-coords1[2],2)));
              for(c = 0; c < 3; c++){
                coordc_array[c] = coords1[c];
                n_cc[c] = n_cc[c];
              }
            }
            else{
              lc += sqrt((pow(coordc_array[0]-coords2[0],2))+(pow(coordc_array[1]-coords2[1],2))+(pow(coordc_array[2]-coords2[2],2)));
              for(c = 0; c < 3; c++)
                coordc_array[c] = coords2[c];
              n_cc[c] = -1*n_cc[c];
            }
            cod[0] = cod[1];
          }
          n_cc[2] = n_cc[1] = n_cc[0] = 0;
          lc = 0;
          ave_V = 0;
          cod[0] = cod_in;
          for(c = 0; c < 3; c++){
            coordc_array[c] = coorda_array[c]-grad_cc[c]*len;
          }
          lc = sqrt((pow(coordc_array[0]-coorda_array[0],2))+(pow(coordc_array[1]-coorda_array[1],2))+(pow(coordc_array[2]-coorda_array[2],2)));
          while(lc < ctx->WidthIntLength && ave_V < 1.0 && (-grad_cc[0]*n_cc[0]-grad_cc[1]*n_cc[1]-grad_cc[2]*n_cc[2] >= 0.)){
            for (ek1 = zs1; ek1 < zs1+zm1; ek1++) {
              for (ej1 = ys1; ej1 < ys1+ym1; ej1++) {
                for (ei1 = xs1; ei1 < xs1+xm1; ei1++) {
                  if(
                     ((coords_array[ek1][ej1][ei1+1][0] >= coordc_array[0]) && (coords_array[ek1][ej1][ei1][0] <= coordc_array[0] ))    &&
                     ((coords_array[ek1][ej1+1][ei1][1] >= coordc_array[1]) && (coords_array[ek1][ej1][ei1][1] <= coordc_array[1] ))    &&
                     ((coords_array[ek1+1][ej1][ei1][2] >= coordc_array[2]) && (coords_array[ek1][ej1][ei1][2] <= coordc_array[2] ))
                     )
                  {
                    ekk = ek1;ejj = ej1;eii = ei1;
                  }
                }
              }
            }
            hx = coords_array[ekk][ejj][eii+1][0]-coords_array[ekk][ejj][eii][0];
            hy = coords_array[ekk][ejj+1][eii][1]-coords_array[ekk][ejj][eii][1];
            hz = coords_array[ekk+1][ejj][eii][2]-coords_array[ekk][ejj][eii][2];
            hwx = coordc_array[0]-coords_array[ekk][ejj][eii][0];
            hwy = coordc_array[1]-coords_array[ekk][ejj][eii][1];
            hwz = coordc_array[2]-coords_array[ekk][ejj][eii][2];
            ierr = CartFEElement3DInit(&ctx->s3D,hwx,hwy,hwz,hx,hy,hz);CHKERRQ(ierr);
            ierr = ComputeAverageVlocal(&ave_V, v_array, ekk, ejj, eii, &ctx->s3D);
            ierr = ComputeUcdotGradVlocal(&cod[1], grad_cc, n_cc, u_array, v_array, ekk, ejj, eii, &ctx->s3D);
            ierr = VFCartFEElement1DInit(&ctx->e1D,len);CHKERRQ(ierr);
            ierr = IntegrateUcdotGradVlocal(&w_array[ek][ej][ei],cod, &ctx->e1D);
            if((n_cc[0] == 0) && (n_cc[1] == 0) && (n_cc[2] == 0)){
              n_cc[0] = grad_cc[0];
              n_cc[1] = grad_cc[1];
              n_cc[2] = grad_cc[2];
            }
            for(c = 0; c < 3; c++){
              coords1[c] = coordc_array[c]+n_cc[c]*len;
              coords2[c] = coordc_array[c]-n_cc[c]*len;
            }
            tlent1 = sqrt((pow(coords1[0]-coorda_array[0],2))+(pow(coords1[1]-coorda_array[1],2))+(pow(coords1[2]-coorda_array[2],2)));
            tlent2 = sqrt((pow(coords2[0]-coorda_array[0],2))+(pow(coords2[1]-coorda_array[1],2))+(pow(coords2[2]-coorda_array[2],2)));
            if(tlent1 > tlent2){
              lc += sqrt((pow(coordc_array[0]-coords1[0],2))+(pow(coordc_array[1]-coords1[1],2))+(pow(coordc_array[2]-coords1[2],2)));
              for(c = 0; c < 3; c++){
                coordc_array[c] = coords1[c];
                n_cc[c] = n_cc[c];
              }
            }
            else{
              lc += sqrt((pow(coordc_array[0]-coords2[0],2))+(pow(coordc_array[1]-coords2[1],2))+(pow(coordc_array[2]-coords2[2],2)));
              for(c = 0; c < 3; c++){
                coordc_array[c] = coords2[c];
                n_cc[c] = -1*n_cc[c];
              }
            }
            cod[0] = cod[1];
          }
          if(w_array[ek][ej][ei] < 0){
            w_array[ek][ej][ei] = 0;
          }
        }
        pmult_array[ek][ej][ei] = pow(w_array[ek][ej][ei],3)/12.;
			}
		}
  }
  ierr = DMDAVecRestoreArrayDOF(ctx->daWVect,coords_local,&coords_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daWVect,&coords_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScalCell,w_local,&w_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&w_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,w_local,ADD_VALUES,fields->widthc);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,w_local,ADD_VALUES,fields->widthc);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daWVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daWVect,&u_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daWScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daWScal,&v_local);CHKERRQ(ierr);
  
  ierr = VecSet(fields->width,0.);CHKERRQ(ierr);
  ierr = CellToNodeInterpolation(fields->width,fields->widthc,ctx); CHKERRQ(ierr);
  
	ierr = DMDAVecRestoreArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeAverageVlocal"
extern PetscErrorCode ComputeAverageVlocal(PetscReal *v_elem, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFEElement3D *s)
{
	PetscInt        k, j, i;
  
	PetscFunctionBegin;
  *v_elem = 0;
  for (k = 0; k < s->nphiz; k++) {
    for (j = 0; j < s->nphiy; j++) {
      for (i = 0; i < s->nphix; i++) {
        *v_elem += v_array[ek+k][ej+j][ei+i] * s->phi[k][j][i];
      }
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IntegrateUcdotGradVlocal"
extern PetscErrorCode IntegrateUcdotGradVlocal(PetscReal *w_ave, PetscReal *cod, VFCartFEElement1D *e)
{
  PetscErrorCode ierr;
	PetscInt        i, eg;
	PetscReal       *cod_elem;
  
	PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&cod_elem);CHKERRQ(ierr);
  for (eg = 0; eg <  e->ng; eg++){
    cod_elem[eg] = 0;
  }
  for (i = 0; i < e->nphix; i++) {
    for (eg = 0; eg <  e->ng; eg++){
      cod_elem[eg] += cod[i]*e->phi[0][0][i][eg];
    }
  }
  for (eg = 0; eg < e->ng; eg++) {
    *w_ave += cod_elem[eg]*e->weight[eg]/2.;
  }
  ierr = PetscFree(cod_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeCellCenterGradV_local"
extern PetscErrorCode ComputeCellCenterGradV_local(PetscReal *grad_cc, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFEElement3D *s)
{
  PetscInt        k, j, i;
	PetscInt        c;
  PetscReal       v_mag;
  
	PetscFunctionBegin;
  for (c = 0; c < 3; c++){
    grad_cc[c] = 0;
  }
  for (k = 0; k < s->nphiz; k++) {
    for (j = 0; j < s->nphiy; j++) {
      for (i = 0; i < s->nphix; i++) {
        for (c = 0; c < 3; c++){
          grad_cc[c] += v_array[ek+k][ej+j][ei+i] * s->dphi[k][j][i][c];
        }
      }
    }
  }
  v_mag = sqrt(pow(grad_cc[0],2)+pow(grad_cc[1],2)+pow(grad_cc[2],2));
  for (c = 0; c < 3; c++){
    grad_cc[c] = grad_cc[c]/v_mag;
  }
  
  
  if((PetscIsInfOrNanScalar(grad_cc[0])) || (PetscIsInfOrNanScalar(grad_cc[1])) || (PetscIsInfOrNanScalar(grad_cc[2])) )
  {
    grad_cc[0] = grad_cc[1] = grad_cc[2] = v_mag = 0;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumetricCrackOpening3D_localCC"
extern PetscErrorCode VolumetricCrackOpening3D_localCC(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ***udotn_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
  PetscReal   *n_elem[3],*dv_mag_elem;
  PetscReal   *u_elem[3],*dv_elem[3];
	PetscReal		element_vol = 0;
  
	PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,&dv_elem[0],
                      e->ng,&dv_elem[1],
                      e->ng,&dv_elem[2],
                      e->ng,&u_elem[0],
                      e->ng,&u_elem[1],
                      e->ng,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,&n_elem[0],e->ng,&n_elem[1],e->ng,&n_elem[2],e->ng,&dv_mag_elem);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
    }
    dv_mag_elem[eg] = 0;
	}
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
				}
			}
		}
	}
	volcrackopening_array[ek][ej][ei] = 0.;
	for(eg = 0; eg < e->ng; eg++){
    element_vol += e->weight[eg];
    for (c = 0; c < 3; c++){
      volcrackopening_array[ek][ej][ei] += u_elem[c][eg]*dv_elem[c][eg]*e->weight[eg];
    }
	}
  if(CrackVolume_local != PETSC_NULL){
    *CrackVolume_local = 0.;
    *CrackVolume_local = volcrackopening_array[ek][ej][ei];
    
  }
  volcrackopening_array[ek][ej][ei] = volcrackopening_array[ek][ej][ei]/element_vol;
  
  if(udotn_array != PETSC_NULL){
    udotn_array[ek][ej][ei] = 0;
    for (eg = 0; eg < e->ng; eg++){
      dv_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
      for(c = 0; c < 3; c++){
        n_elem[c][eg] = PetscAbs(dv_elem[c][eg]);
      }
      if((PetscIsInfOrNanScalar(n_elem[0][eg])) || (PetscIsInfOrNanScalar(n_elem[1][eg])) || (PetscIsInfOrNanScalar(n_elem[2][eg])) || (PetscIsInfOrNanScalar(dv_mag_elem[eg])) )
      {
        n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = dv_mag_elem[eg] = 0;
      }
      for(c = 0; c < 3; c++){
        
        udotn_array[ek][ej][ei] += u_elem[c][eg]*n_elem[c][eg]*e->weight[eg];
      }
    }
    udotn_array[ek][ej][ei] = udotn_array[ek][ej][ei]/element_vol;
  }
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],dv_mag_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumetricStrainVolume_local"
extern PetscErrorCode VolumetricStrainVolume_local(PetscReal *VolStrainVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, VFMatProp *matprop, PetscReal ****u_diff_array, PetscReal ***v_array)
{
  
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
	PetscReal    *v_elem,*du_elem[3],beta;
  
	PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,&v_elem,
                      e->ng,&du_elem[0],
                      e->ng,&du_elem[1],
                      e->ng,&du_elem[2]);CHKERRQ(ierr);
  beta  = matprop->beta;
  for (eg = 0; eg < e->ng; eg++){
    v_elem[eg] = 0.;
    for (c = 0; c < 3; c++){
      du_elem[c][eg] = 0.;
    }
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
          for (c = 0; c < 3; c++){
            du_elem[c][eg] += u_diff_array[ek+k][ej+j][ei+i][c]*e->dphi[k][j][i][c][eg];
          }
        }
      }
    }
  }
	*VolStrainVolume_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      *VolStrainVolume_local += 1*beta*(pow(v_elem[eg],2))*du_elem[c][eg]*e->weight[eg];
    }
	}
	ierr = PetscFree4(v_elem,du_elem[0],du_elem[1],du_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFCheckVolumeBalance"
extern PetscErrorCode VFCheckVolumeBalance(PetscReal *ModulusVolume, PetscReal *DivVolume, PetscReal *SurfVolume, PetscReal *SumWellRate,PetscReal *SumSourceRate,PetscReal *VolStrainVolume,VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei,i;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ****vel_array;
  PetscReal       ***v_array;
  Vec             v_local;
	Vec             vel_local;
  Vec             Pressure_diff;
  Vec             press_diff_local;
  PetscReal       ***press_diff_array;
  PetscReal       ***src_array;
	Vec             src_local;
  PetscReal       ****u_diff_array;
	Vec             U_diff,u_diff_local;
  PetscReal      ***one_array;
  Vec            one_local;
  Vec            Ones;
  PetscReal       mymodVolumeLocal = 0.,mymodVolume = 0.;
  PetscReal       mydivVolumeLocal = 0.,mydivVolume = 0.;
  PetscReal       mysurfVolumeLocal = 0.,mysurfVolume = 0.;
  PetscReal       mysourceVolumeLocal = 0.,mysourceVolume = 0.;
  PetscReal       mystrainVolumeLocal = 0.,mystrainVolume = 0.;
  FACE           face;
  PetscReal       timestepsize = 0;
  PetscReal      ***m_inv_array;
  Vec            m_inv_local;
  
  PetscFunctionBegin;
  timestepsize     = ctx->timevalue;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&src_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->Source,INSERT_VALUES,src_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->Source,INSERT_VALUES,src_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,src_local,&src_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVect,&vel_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->velocity,INSERT_VALUES,vel_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->velocity,INSERT_VALUES,vel_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,vel_local,&vel_array);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->daScal,&Pressure_diff);CHKERRQ(ierr);
  ierr = VecSet(Pressure_diff,0.);CHKERRQ(ierr);
  ierr = VecAXPY(Pressure_diff,-1.0,ctx->pressure_old);CHKERRQ(ierr);
  ierr = VecAXPY(Pressure_diff,1.0,fields->pressure);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&press_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,Pressure_diff,INSERT_VALUES,press_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,Pressure_diff,INSERT_VALUES,press_diff_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,press_diff_local,&press_diff_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daVect,&U_diff);CHKERRQ(ierr);
  ierr = VecSet(U_diff,0.);CHKERRQ(ierr);
  ierr = VecAXPY(U_diff,-1.0,ctx->U_old);CHKERRQ(ierr);
  ierr = VecAXPY(U_diff,1.0,fields->U);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daVect,&u_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,U_diff,INSERT_VALUES,u_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,U_diff,INSERT_VALUES,u_diff_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_diff_local,&u_diff_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daScal,&Ones);CHKERRQ(ierr);
  ierr = VecSet(Ones,1.0);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScal,&one_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,Ones,INSERT_VALUES,one_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,Ones,INSERT_VALUES,one_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,one_local,&one_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScalCell,&m_inv_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScalCell,ctx->M_inv,INSERT_VALUES,m_inv_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScalCell,ctx->M_inv,INSERT_VALUES,m_inv_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScalCell,m_inv_local,&m_inv_array);CHKERRQ(ierr);
  
	*ModulusVolume = 0.;
	*DivVolume = 0.;
  *SurfVolume = 0;
  *SumWellRate = 0;
  *SumSourceRate = 0;
  *VolStrainVolume = 0;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = ModulusVolume_local(&mymodVolumeLocal, ek, ej, ei,&ctx->e3D,m_inv_array[ek][ej][ei],press_diff_array,v_array);CHKERRQ(ierr);
				ierr = DivergenceVolume_local(&mydivVolumeLocal, ek, ej, ei, &ctx->e3D,vel_array, v_array);CHKERRQ(ierr);
        if(ctx->hasFluidSources){
          ierr = SourceVolume_local(&mysourceVolumeLocal, ek, ej, ei, &ctx->e3D, src_array, v_array);CHKERRQ(ierr);
        }
        if(ctx->FlowDisplCoupling){
          ierr = VolumetricStrainVolume_local(&mystrainVolumeLocal, ek, ej, ei, &ctx->e3D,&ctx->matprop[ctx->layer[ek]],u_diff_array,v_array);CHKERRQ(ierr);
        }
        mymodVolume += mymodVolumeLocal;
        mysourceVolume += timestepsize*mysourceVolumeLocal;
        mydivVolume += timestepsize*mydivVolumeLocal;
        mystrainVolume += mystrainVolumeLocal;
        if(ei == 0){
          face = X0;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += timestepsize*mysurfVolumeLocal;
        }
        if(ei == nx-1){
          face = X1;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += timestepsize*mysurfVolumeLocal;
        }
        if(ej == 0){
          face = Y0;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += timestepsize*mysurfVolumeLocal;
        }
        if(ej == ny-1){
          face = Y1;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += timestepsize*mysurfVolumeLocal;
        }
        if(ek == 0){
          face = Z0;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += timestepsize*mysurfVolumeLocal;
        }
        if(ek == nz-1){
          face = Z1;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += timestepsize*mysurfVolumeLocal;
        }
			}
		}
	}
  *SumWellRate = 0;
  for (i = 0; i < ctx->numWells; i++) {
    if(ctx->well[i].type == INJECTOR){
      *SumWellRate = *SumWellRate + timestepsize*ctx->well[i].Qw;
    }
    else{
      *SumWellRate = *SumWellRate - timestepsize*ctx->well[i].Qw;
    }
  }
	ierr = MPI_Allreduce(&mymodVolume,ModulusVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mydivVolume,DivVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mysurfVolume,SurfVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mysourceVolume,SumSourceRate,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mystrainVolume,VolStrainVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,press_diff_local,&press_diff_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&press_diff_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,src_local,&src_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&src_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,vel_local,&vel_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&vel_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_diff_local,&u_diff_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&u_diff_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,one_local,&one_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&one_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScalCell,m_inv_local,&m_inv_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScalCell,&m_inv_local);CHKERRQ(ierr);
  ierr = VecDestroy(&Ones);CHKERRQ(ierr);
  ierr = VecDestroy(&Pressure_diff);CHKERRQ(ierr);
  ierr = VecDestroy(&U_diff);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModulusVolume_local"
extern PetscErrorCode ModulusVolume_local(PetscReal *ModVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, PetscReal m_inv, PetscReal ***press_diff_array,PetscReal ***v_array)
{
  
	PetscErrorCode ierr;
	PetscInt		i, j, k;
	PetscInt		eg;
	PetscReal		*v_elem,*press_diff_elem;
  
	PetscFunctionBegin;
  
  ierr = PetscMalloc2(e->ng,&v_elem,e->ng,&press_diff_elem);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
    v_elem[eg] = 0;
    press_diff_elem[eg] = 0;
	}
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
          v_elem[eg] += v_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
          press_diff_elem[eg] += press_diff_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
				}
			}
		}
	}
	*ModVolume_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
		*ModVolume_local += (pow(v_elem[eg],2))*m_inv*press_diff_elem[eg]*e->weight[eg];
	}
	ierr = PetscFree2(v_elem,press_diff_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SurfaceFluxVolume_local"
extern PetscErrorCode SurfaceFluxVolume_local(PetscReal *mysurfVolumeLocal,PetscInt ek, PetscInt ej, PetscInt ei, FACE face, VFCartFEElement2D *e, PetscReal ****vel_array, PetscReal ***v_array)
{
  
	PetscErrorCode ierr;
	PetscInt		i, j, k;
	PetscInt		g;
	PetscReal   *flux_elem,*v_elem;
  
	PetscFunctionBegin;
  ierr = PetscMalloc2(e->ng,&v_elem,e->ng,&flux_elem);CHKERRQ(ierr);
	for (g = 0; g < e->ng; g++){
    v_elem[g] = 0;
    flux_elem[g] = 0;
	}
  *mysurfVolumeLocal = 0.;
  switch (face) {
		case X0:
			for (k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							v_elem[g] += e->phi[i][j][k][g]*v_array[ek+k][ej+j][ei+i];
							flux_elem[g] += -1*e->phi[i][j][k][g]*vel_array[ek+k][ej+j][ei+i][0];
						}
					}
				}
			}
			break;
    case X1:
			for (k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							v_elem[g] += e->phi[i][j][k][g]*v_array[ek+k][ej+j][ei+1];
							flux_elem[g] += e->phi[i][j][k][g]*vel_array[ek+k][ej+j][ei+1][0];
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
              v_elem[g] += e->phi[j][k][i][g]*v_array[ek+k][ej+j][ei+i];
              flux_elem[g] += -1*e->phi[j][k][i][g]*vel_array[ek+k][ej+j][ei+i][1];
						}
					}
				}
			}
			break;
    case Y1:
			for (k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
              v_elem[g] += e->phi[j][k][i][g]*v_array[ek+k][ej+1][ei+i];
              flux_elem[g] += e->phi[j][k][i][g]*vel_array[ek+k][ej+1][ei+i][1];
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
              v_elem[g] += e->phi[k][j][i][g]*v_array[ek+k][ej+i][ei+i];
              flux_elem[g] += -1*e->phi[k][j][i][g]*vel_array[ek+k][ej+i][ei+i][2];
						}
					}
				}
			}
      break;
    case Z1:
			for (k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
              v_elem[g] += e->phi[k][j][i][g]*v_array[ek+1][ej+i][ei+i];
              flux_elem[g] += e->phi[k][j][i][g]*vel_array[ek+1][ej+i][ei+i][2];
						}
					}
				}
			}
  }
  for(g = 0; g < e->ng; g++){
    *mysurfVolumeLocal += (pow(v_elem[g],2))*flux_elem[g]*e->weight[g];
  }
  ierr = PetscFree2(v_elem,flux_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DivergenceVolume_local"
extern PetscErrorCode DivergenceVolume_local(PetscReal *DivVolume_local, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, PetscReal ****vel_array,PetscReal ***v_array)
{
  
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
	PetscReal		*dvel_elem[3],*v_elem;
	PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,&dvel_elem[0],e->ng,&dvel_elem[1],e->ng,&dvel_elem[2],e->ng,&v_elem);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
    v_elem[eg] = 0;
    for (c = 0; c < 3; c++){
      dvel_elem[c][eg] = 0;
    }
	}
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
          v_elem[eg] += v_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
          for(c = 0; c < 3; c++){
            dvel_elem[c][eg] += vel_array[ek+k][ej+j][ei+i][c]*e->dphi[k][j][i][c][eg];
          }
				}
			}
		}
	}
	*DivVolume_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      *DivVolume_local += (pow(v_elem[eg],2))*dvel_elem[c][eg]*e->weight[eg];
    }
	}
	ierr = PetscFree4(dvel_elem[0],dvel_elem[1],dvel_elem[2],v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SourceVolume_local"
extern PetscErrorCode SourceVolume_local(PetscReal *SrcVolume_local, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e,PetscReal ***src_array, PetscReal ***v_array)
{
  
	PetscErrorCode ierr;
	PetscInt		i, j, k;
	PetscInt		eg;
	PetscReal		*src_elem,*v_elem;
  
	PetscFunctionBegin;
  ierr = PetscMalloc2(e->ng,&v_elem,e->ng,&src_elem);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
    v_elem[eg] = 0;
    src_elem[eg] = 0;
	}
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
          v_elem[eg] += v_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
          src_elem[eg] += src_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
				}
			}
		}
	}
	*SrcVolume_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
		*SrcVolume_local += (pow(v_elem[eg],2))*src_elem[eg]*e->weight[eg];
	}
	ierr = PetscFree2(v_elem,src_elem);CHKERRQ(ierr);
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
	Vec             node_local;
  PetscReal       ***node_array;
	PetscReal       ****node_arraydof;
	PetscReal       ***cell_array;
	PetscReal       ****cell_arraydof;
	DM              dm;
	
	PetscFunctionBegin;
  ierr = VecGetDM(node_vec,&dm);CHKERRQ(ierr);
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
	else if (dof == 3)
	{
		ierr = DMGetLocalVector(dm,&node_local);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(dm,node_vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(dm,node_vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
		ierr = DMDAVecGetArrayDOF(dm, node_local,&node_arraydof);CHKERRQ(ierr);
		ierr = DMDAVecGetArrayDOF(ctx->daVectCell,cell_vec,&cell_arraydof);CHKERRQ(ierr);
	}
  else
  {
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
  else if (dof == 3)
	{
		ierr = DMDAVecRestoreArrayDOF(dm,node_vec,&node_arraydof);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArrayDOF(ctx->daVectCell,cell_vec,&cell_arraydof);CHKERRQ(ierr);
	}
  else
	{
		ierr = DMDAVecRestoreArrayDOF(dm,node_vec,&node_arraydof);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,cell_vec,&cell_arraydof);CHKERRQ(ierr);
	}
  ierr = VecDestroy(&volume);CHKERRQ(ierr);
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
  PetscReal       timestepsize;
  
  PetscFunctionBegin;
  timestepsize = ctx->timevalue;
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
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricCrackOpening3D_localCC(&myLeakOffRateLocal, volleakoffrate_array,PETSC_NULL, q_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        myLeakOffRate += timestepsize*myLeakOffRateLocal;
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
extern PetscErrorCode VolumetricLeakOffRate_local(PetscReal *LeakoffRate_local, PetscReal ***volleakoffrate_array, PetscReal ****q_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
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
	ierr = PetscMalloc3(e->ng,&dx_vfield_loc,e->ng,&dy_vfield_loc,e->ng,&dz_vfield_loc);CHKERRQ(ierr);
	ierr = PetscMalloc3(e->ng,&velx_loc,e->ng,&vely_loc,e->ng,&velz_loc);CHKERRQ(ierr);
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
#define __FUNCT__ "VolumetricFunction_local"
extern PetscErrorCode VolumetricFunction_local(PetscReal *Function_local, PetscReal ***value_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
{
	PetscErrorCode ierr;
	PetscInt        i, j, k,c;
	PetscInt        eg;
  PetscReal       *dv_elem[3],*u_elem[3],*n_elem[3],*n_mag_elem;
  
	PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,&dv_elem[0],e->ng,&dv_elem[1],e->ng,&dv_elem[2],e->ng,&u_elem[0],e->ng,&u_elem[1],e->ng,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,&n_elem[0],e->ng,&n_elem[1],e->ng,&n_elem[2],e->ng,&n_mag_elem);CHKERRQ(ierr);
  
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
#define __FUNCT__ "VolumetricFractureWellRate_local"
extern PetscErrorCode VolumetricFractureWellRate_local(PetscReal *InjVolume_local, PetscReal ***regrate_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k,c;
	PetscInt		eg;
  PetscReal    *dv_elem[3],*regrate_elem;
  
	PetscFunctionBegin;
	ierr = PetscMalloc4(e->ng,&dv_elem[0],e->ng,&dv_elem[1],e->ng,&dv_elem[2],e->ng,&regrate_elem);CHKERRQ(ierr);
  
	for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
    }
		regrate_elem[eg] = 0.;
	}
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
          }
          regrate_elem[eg] += regrate_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
        }
      }
    }
  }
  *InjVolume_local = 0.;
  
  for(eg = 0; eg < e->ng; eg++){
    *InjVolume_local += regrate_elem[eg]*(sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2))))*e->weight[eg];
  }
	ierr = PetscFree4(dv_elem[0],dv_elem[1],dv_elem[2],regrate_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumetricFractureWellRate"
extern PetscErrorCode VolumetricFractureWellRate(PetscReal *InjectedVolume, VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ***regrate_array;
	Vec             regrate_local;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscReal       myInjVolumeRateLocal = 0.,myInjVolumeRate = 0.;
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&regrate_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->RegFracWellFlowRate,INSERT_VALUES,regrate_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->RegFracWellFlowRate,INSERT_VALUES,regrate_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,regrate_local,&regrate_array);CHKERRQ(ierr);
  
	*InjectedVolume = 0.;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricFractureWellRate_local(&myInjVolumeRateLocal, regrate_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        myInjVolumeRate += myInjVolumeRateLocal;
			}
		}
	}
	ierr = MPI_Allreduce(&myInjVolumeRate,InjectedVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,regrate_local,&regrate_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&regrate_local);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_IntegrateOnBoundary_local"
extern PetscErrorCode VF_IntegrateOnBoundary_local(PetscReal *SumnIntegral_local,PetscReal ***node_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,VFCartFEElement2D *e)
{
  PetscErrorCode ierr;
	PetscInt		i, j, k, g;
  PetscReal		*node_elem;
  
  PetscFunctionBegin;
	ierr = PetscMalloc(e->ng*sizeof(PetscReal),&node_elem);CHKERRQ(ierr);
  for (g = 0; g < e->ng; g++){
		node_elem[g] = 0.;
	}
  *SumnIntegral_local = 0.;
  switch (face) {
		case X0:
			for (k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							node_elem[g] += e->phi[i][j][k][g]*node_array[ek+k][ej+j][ei+i];
						}
					}
				}
			}
			break;
    case X1:
			for (k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							node_elem[g] += e->phi[i][j][k][g]*node_array[ek+k][ej+j][ei+1];
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
              node_elem[g] += e->phi[j][k][i][g]*node_array[ek+k][ej+j][ei+i];
						}
					}
				}
			}
			break;
    case Y1:
			for (k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
              node_elem[g] += e->phi[j][k][i][g]*node_array[ek+k][ej+1][ei+i];
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
              node_elem[g] += e->phi[k][j][i][g]*node_array[ek+k][ej+i][ei+i];
						}
					}
				}
			}
      break;
    case Z1:
			for (k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
              node_elem[g] += e->phi[k][j][i][g]*node_array[ek+1][ej+i][ei+i];
						}
					}
				}
			}
  }
  for(g = 0; g < e->ng; g++){
    *SumnIntegral_local += node_elem[g]*e->weight[g];
  }
  ierr = PetscFree(node_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_IntegrateOnBoundary"
extern PetscErrorCode VF_IntegrateOnBoundary(PetscReal *SumnIntegral,Vec vec, FACE face, VFCtx *ctx)
{
  PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
  Vec             node_local;
  PetscReal       ***node_array;
  PetscInt        dof;
	DM              dm;
  PetscReal       SumnIntegralLocal = 0.,mySumnIntegral = 0.;
  
  PetscFunctionBegin;
	ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dm,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&node_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,vec,INSERT_VALUES,node_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(dm,node_local,&node_array);CHKERRQ(ierr);
  
  
	*SumnIntegral = 0.;
  switch (face) {
    case X0:
      ei = 0;
      for (ek = zs; ek < zs+zm; ek++) {
        for (ej = ys; ej < ys+ym; ej++) {
          hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
          hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
          ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          ierr = VF_IntegrateOnBoundary_local(&SumnIntegralLocal,node_array,ek,ej,ei,face,&ctx->e2D);CHKERRQ(ierr);
          mySumnIntegral += SumnIntegralLocal;
        }
      }
      break;
    case X1:
      ei = nx-1;
      for (ek = zs; ek < zs+zm; ek++) {
        for (ej = ys; ej < ys+ym; ej++) {
          hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
          hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
          ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          ierr = VF_IntegrateOnBoundary_local(&SumnIntegralLocal,node_array,ek,ej,ei,face,&ctx->e2D);CHKERRQ(ierr);
          mySumnIntegral += SumnIntegralLocal;
        }
      }
      break;
    case Y0:
      ej = 0;
      for (ek = zs; ek < zs+zm; ek++) {
        for (ei = xs; ei < xs+xm; ei++) {
          hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
          hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
          ierr = VF_IntegrateOnBoundary_local(&SumnIntegralLocal,node_array,ek,ej,ei,face,&ctx->e2D);CHKERRQ(ierr);
          mySumnIntegral += SumnIntegralLocal;
        }
      }
      break;
    case Y1:
      ej = ny-1;
      for (ek = zs; ek < zs+zm; ek++) {
        for (ei = xs; ei < xs+xm; ei++) {
          hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
          hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
          ierr = VF_IntegrateOnBoundary_local(&SumnIntegralLocal,node_array,ek,ej,ei,face,&ctx->e2D);CHKERRQ(ierr);
          mySumnIntegral += SumnIntegralLocal;
        }
      }
      break;
    case Z0:
      ek = 0;
      for (ej = ys; ej < ys+ym; ej++) {
        for (ei = xs; ei < xs+xm; ei++) {
          hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
          hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
          ierr = VF_IntegrateOnBoundary_local(&SumnIntegralLocal,node_array,ek,ej,ei,face,&ctx->e2D);CHKERRQ(ierr);
          mySumnIntegral += SumnIntegralLocal;
        }
      }
      break;
    case Z1:
      ek = nz-1;
      for (ej = ys; ej < ys+ym; ej++) {
        for (ei = xs; ei < xs+xm; ei++) {
          hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
          hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
          ierr = VF_IntegrateOnBoundary_local(&SumnIntegralLocal,node_array,ek,ej,ei,face,&ctx->e2D);CHKERRQ(ierr);
          mySumnIntegral += SumnIntegralLocal;
        }
      }
      break;
  }
  ierr = MPI_Allreduce(&mySumnIntegral,SumnIntegral,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(dm, node_local,&node_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&node_local);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_ComputeRegularizedFracturePressure_local"
extern PetscErrorCode VF_ComputeRegularizedFracturePressure_local(PetscReal ***press_c_array, PetscReal ***press_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
  PetscReal   *pre_elem;
  PetscReal   *u_elem[3],*dv_elem[3];
	PetscReal		element_vol = 0;
  
	PetscFunctionBegin;
  ierr = PetscMalloc7(e->ng,&dv_elem[0],e->ng,&dv_elem[1],e->ng,&dv_elem[2],
                      e->ng,&u_elem[0],e->ng,&u_elem[1],e->ng,&u_elem[2],
                      e->ng,&pre_elem);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0.;
      u_elem[c][eg] = 0;
    }
    pre_elem[eg] = 0;
	}
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
          for (c = 0; c < 3; c++){
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
          }
          pre_elem[eg] += press_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
				}
			}
		}
	}
	press_c_array[ek][ej][ei] = 0.;
	for(eg = 0; eg < e->ng; eg++){
    for (c = 0; c < 3; c++){
      press_c_array[ek][ej][ei] += pre_elem[eg]*u_elem[c][eg]*dv_elem[c][eg]*e->weight[eg];
      element_vol += e->weight[eg];
    }
  }
  press_c_array[ek][ej][ei] = press_c_array[ek][ej][ei]/element_vol;
	ierr = PetscFree7(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2],pre_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_ComputeRegularizedFracturePressure"
extern PetscErrorCode VF_ComputeRegularizedFracturePressure(VFCtx *ctx, VFFields *fields)
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
  PetscReal       ***press_array;
  Vec             press_local;
  Vec             Pressure_cell;
  Vec             press_c_local;
  PetscReal       ***press_c_array;
  
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
  ierr = DMCreateGlobalVector(ctx->daScalCell,&Pressure_cell);CHKERRQ(ierr);
  ierr = VecSet(Pressure_cell,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScalCell,&press_c_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,Pressure_cell,INSERT_VALUES,press_c_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,Pressure_cell,INSERT_VALUES,press_c_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,press_c_local,&press_c_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&press_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->pressure,INSERT_VALUES,press_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->pressure,INSERT_VALUES,press_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,press_local,&press_array);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VF_ComputeRegularizedFracturePressure_local(press_c_array, press_array, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
			}
		}
	}
  ierr = DMDAVecRestoreArray(ctx->daScal,press_local,&press_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&press_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScalCell,press_c_local,&press_c_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&press_c_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,press_c_local,ADD_VALUES,Pressure_cell);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,press_c_local,ADD_VALUES,Pressure_cell);CHKERRQ(ierr);
  ierr = VecSet(fields->fracpressure,0.);CHKERRQ(ierr);
  ierr = CellToNodeInterpolation(fields->fracpressure,Pressure_cell,ctx); CHKERRQ(ierr);
  ierr = VecDestroy(&Pressure_cell);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumeFromWidth"
extern PetscErrorCode VolumeFromWidth(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ***v_array;
	Vec             v_local;
	PetscReal       myCrackVolumeLocal = 0.,myCrackVolume = 0.;
  PetscReal      ***w_array;
  Vec            w_local;
  
	
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScalCell,&w_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScalCell,fields->widthc,INSERT_VALUES,w_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScalCell,fields->widthc,INSERT_VALUES,w_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScalCell,w_local,&w_array);CHKERRQ(ierr);
  
	*CrackVolume = 0.;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumeFromWidth_local(&myCrackVolumeLocal, w_array[ek][ej][ei], v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        myCrackVolume += myCrackVolumeLocal;
			}
		}
	}
	ierr = MPI_Allreduce(&myCrackVolume,CrackVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScalCell,w_local,&w_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&w_local);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VolumeFromWidth_local"
extern PetscErrorCode VolumeFromWidth_local(PetscReal *CrackVolume_local, PetscReal w, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
  PetscReal   *dv_elem[3],*dv_mag_elem;
  
  
	PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,&dv_elem[0],
                      e->ng,&dv_elem[1],
                      e->ng,&dv_elem[2],
                      e->ng,&dv_mag_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for (c = 0; c < 3; c++){
      dv_elem[c][eg] = 0.;
    }
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c= 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
    dv_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
  }
	*CrackVolume_local = 0.;
	for(eg = 0; eg < e->ng; eg++){
		*CrackVolume_local += w*dv_mag_elem[eg]*e->weight[eg];
	}
  ierr = PetscFree4(dv_elem[0],dv_elem[1],dv_elem[2],dv_mag_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "UpdatePermeablitysingMultipliers"
extern PetscErrorCode UpdatePermeablitysingMultipliers(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
	PetscReal       ***v_array;
	Vec             v_local;
  PetscReal       ave_V = 0;
  Vec             pmult_local;
	PetscReal       ***pmult_array;
  PetscReal       ****perm_array;
  Vec             perm_local;
  PetscReal       ****perm1_array;
  Vec             perm1_local;
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daVFperm,&perm1_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daVFperm,ctx->Perm,INSERT_VALUES,perm1_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daVFperm,ctx->Perm,INSERT_VALUES,perm1_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm1_local,&perm1_array);CHKERRQ(ierr);
  
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
        ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        ierr = CartFEElement3DInit(&ctx->s3D,hx/2.,hy/2.,hz/2.,hx,hy,hz);CHKERRQ(ierr);
        ierr = ComputeAverageVlocal(&ave_V, v_array, ek, ej, ei, &ctx->s3D);
        if(ave_V < ctx->pmult_vtol){
          perm_array[ek][ej][ei][0] = perm1_array[ek][ej][ei][0]*pmult_array[ek][ej][ei];
          perm_array[ek][ej][ei][1] = perm1_array[ek][ej][ei][1]*pmult_array[ek][ej][ei];
          perm_array[ek][ej][ei][2] = perm1_array[ek][ej][ei][2]*pmult_array[ek][ej][ei];
        }
      }
    }
  }
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ctx->daVFperm,perm_local,INSERT_VALUES,fields->vfperm);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ctx->daVFperm,perm_local,INSERT_VALUES,fields->vfperm);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm1_local,&perm1_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVFperm,&perm1_local);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ctx->daVFperm,perm1_local,INSERT_VALUES,ctx->Perm);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ctx->daVFperm,perm1_local,INSERT_VALUES,ctx->Perm);CHKERRQ(ierr);
  
	ierr = DMDAVecRestoreArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,pmult_local,INSERT_VALUES,fields->pmult);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}
