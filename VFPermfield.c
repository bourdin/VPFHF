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
#define __FUNCT__ "VolumetricCrackOpening"
extern PetscErrorCode VolumetricCrackOpening(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    k, j, i, c;
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
  Vec             Cellvfield;
  Vec             v_c_local;
	PetscReal       ***v_c_array;
  Vec             pmult_local;
	PetscReal       ***pmult_array;
  PetscReal       gradx,grady,gradz;
  PetscReal       sumx,sumy,sumz;
  Vec             Udotn;
  Vec             udotn_local;
	PetscReal       ***udotn_array;
  Vec             Uc;
  Vec             uc_local;
	PetscReal       ****uc_array;
  Vec             uv_local;
	PetscReal       ****uv_array;
  
  
  
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
  ierr = DMCreateGlobalVector(ctx->daScalCell,&Cellvfield);CHKERRQ(ierr);
  
  ierr = VecSet(Cellvfield,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScalCell,&v_c_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,Cellvfield,INSERT_VALUES,v_c_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,Cellvfield,INSERT_VALUES,v_c_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,v_c_local,&v_c_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daScalCell,&Udotn);CHKERRQ(ierr);
  ierr = VecSet(Udotn,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScalCell,&udotn_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,Udotn,INSERT_VALUES,udotn_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,Udotn,INSERT_VALUES,udotn_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,udotn_local,&udotn_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daVectCell,&Uc);CHKERRQ(ierr);
  ierr = VecSet(fields->Uc,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVectCell,&uc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVectCell,fields->Uc,INSERT_VALUES,uc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVectCell,fields->Uc,INSERT_VALUES,uc_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVectCell,uc_local,&uc_array);CHKERRQ(ierr);
  
 	*CrackVolume = 0.;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        ierr = VF_ComputeAverageVField_local(v_c_array, v_array, ek, ej, ei, &ctx->e3D);
        if(v_c_array[ek][ej][ei] < 0.04){
          v_c_array[ek][ej][ei] = 0;
        }
				ierr = VolumetricCrackOpening3D_local(&myCrackVolumeLocal, volcrackopening_array, udotn_array, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        myCrackVolume += myCrackVolumeLocal;
        ierr = VF_FractureDisplacement_local(uc_array, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        
			}
		}
	}
	ierr = MPI_Allreduce(&myCrackVolume,CrackVolume,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScalCell,VolCrackOpening_local,&volcrackopening_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&VolCrackOpening_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,VolCrackOpening_local,ADD_VALUES,CellVolCrackOpening);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,VolCrackOpening_local,ADD_VALUES,CellVolCrackOpening);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScalCell,v_c_local,&v_c_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&v_c_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,v_c_local,ADD_VALUES,Cellvfield);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,v_c_local,ADD_VALUES,Cellvfield);CHKERRQ(ierr);
  
  
  
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
  
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVectCell,uc_local,&uc_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVectCell,&uc_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daVectCell,uc_local,INSERT_VALUES,fields->Uc);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daVectCell,uc_local,INSERT_VALUES,fields->Uc);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScalCell,udotn_local,&udotn_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&udotn_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,udotn_local,INSERT_VALUES,Udotn);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,udotn_local,INSERT_VALUES,Udotn);CHKERRQ(ierr);
  
  
  
 
  ierr = VecSet(fields->VolCrackOpening,0.);CHKERRQ(ierr);
	ierr = CellToNodeInterpolation(fields->VolCrackOpening,CellVolCrackOpening,ctx); CHKERRQ(ierr);
  
  ierr = VF_IntegrateDisplacement_local(fields->Uv,fields->Uc,ctx);CHKERRQ(ierr);

  ierr = VecSet(fields->Ud,0.);CHKERRQ(ierr);
	ierr = CellToNodeInterpolation(fields->Ud,fields->Uv,ctx); CHKERRQ(ierr);
  
  
  ierr = VecSet(fields->Ue,0.);CHKERRQ(ierr);
	ierr = CellToNodeInterpolation(fields->Ue,fields->Uc,ctx); CHKERRQ(ierr);

  
  
  
  
  ierr = DMGetLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);

  
  
  
  
  ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->Ud,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->Ud,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);


  
  for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
                
				ierr = VolumetricCrackOpening3D_local(PETSC_NULL, pmult_array, PETSC_NULL, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
/*
        if((ek == 48 || ek == 47 || ek == 46 || ek == 45 || ek == 44 || ek == 43|| ek == 42|| ek == 41|| ek == 40|| ek == 39|| ek == 38|| ek == 37) && ej == 0 && ei == 50){
          ierr = PetscPrintf(PETSC_COMM_WORLD," ek = %d, ej = %d, ei %= %d, pmult = %e \t ux = %e \t uy = %e \t uz = %e \t v = %e\n", ek,ej,ei,pmult_array[ek][ej][ei],u_array[ek][ej][ei][0],u_array[ek][ej][ei][1],u_array[ek][ej][ei][2],v_array[ek][ej][ei]);CHKERRQ(ierr);

        
        }
*/
			}
		}
	}
  
  
  
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);  

  ierr = DMDAVecRestoreArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScalCell,pmult_local,ADD_VALUES,fields->pmult);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScalCell,pmult_local,ADD_VALUES,fields->pmult);CHKERRQ(ierr);
  
  
  
  ierr = VecSet(fields->velocity,0.);CHKERRQ(ierr);
	ierr = CellToNodeInterpolation(fields->theta,fields->pmult,ctx); CHKERRQ(ierr);
  
  
  ierr = VecSet(fields->pressure,0.);CHKERRQ(ierr);
  ierr = VecCopy(fields->theta,fields->pressure);CHKERRQ(ierr);
  ierr = VecAYPX(fields->pressure,-1,fields->VolCrackOpening);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = VecDestroy(&CellVolCrackOpening);CHKERRQ(ierr);
  ierr = VecDestroy(&Udotn);CHKERRQ(ierr);
  ierr = VecDestroy(&Uc);CHKERRQ(ierr);
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
	PetscReal       hx,hy,hz;
  Vec             uv_local;
  Vec             uc_local;
  PetscReal       ****uv_array;
	PetscReal       ****uc_array;
	DM              dm;
  PetscReal       sum;
	
	PetscFunctionBegin;
	ierr = PetscObjectQuery((PetscObject)Uv,"DM",(PetscObject*)&dm);CHKERRQ(ierr);
	if (!dm) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
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
  n = 10;
  
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
#define __FUNCT__ "VF_FractureDisplacement_local"
extern PetscErrorCode VF_FractureDisplacement_local(PetscReal ****uc_array,PetscReal ****u_array,PetscReal ***v_array,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,g,c;
  PetscReal      *uc_elem[3],*n_elem[3],*v_elem;
  PetscReal      element_vol,v_mag_elem = 0;
  
  PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&uc_elem[0],e->ng,PetscReal,&uc_elem[1],e->ng,PetscReal,&uc_elem[2]);CHKERRQ(ierr);
  for (g = 0; g < e->ng; g++) {
    v_elem[g] = 0.;
    for(c = 0; c < 3; c++){
      n_elem[c][g] = 0.;
      uc_elem[c][g] = 0;
    }
  }
  for (k = 0; k < e->nphix; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphiz; i++) {
        for (g = 0; g < e->ng; g++) {
          v_elem[g] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][g];
          for(c = 0; c < 3; c++){
            n_elem[c][g] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][g];
            uc_elem[c][g] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][g];
            
          }
        }
      }
    }
  }
  element_vol = 0;
  for (g = 0; g < e->ng; g++){
    v_mag_elem = sqrt((pow(n_elem[0][g],2))+(pow(n_elem[1][g],2))+(pow(n_elem[2][g],2)));
    for(c = 0; c < 3; c++){
      n_elem[c][g] = n_elem[c][g]/v_mag_elem;
    }
    if((PetscIsInfOrNanScalar(n_elem[0][g])) || (PetscIsInfOrNanScalar(n_elem[1][g])) || (PetscIsInfOrNanScalar(n_elem[2][g])) )
    {
      n_elem[0][g] = n_elem[1][g] = n_elem[2][g] = 0;
    }
    element_vol += e->weight[g];
  }
  for(c = 0; c < 3; c++){
    uc_array[ek][ej][ei][c] = 0;
    for(g = 0; g < e->ng; g++){
      n_elem[c][g] = PetscAbs(n_elem[c][g]);
      uc_array[ek][ej][ei][c] += uc_elem[c][g]*(1-v_elem[g])*(1-v_elem[g])*n_elem[c][g]*e->weight[g];
//      uc_array[ek][ej][ei][c] += uc_elem[c][g]*(1-v_elem[g])*(1-v_elem[g])*e->weight[g];
    }
    uc_array[ek][ej][ei][c] = uc_array[ek][ej][ei][c]/element_vol/4;
	}
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],v_elem);CHKERRQ(ierr);
  ierr = PetscFree3(uc_elem[0],uc_elem[1],uc_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumetricCrackOpening3D_local"
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ***udotn_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
  PetscReal   *n_elem[3],*dv_mag_elem;
  PetscReal   *u_elem[3],*dv_elem[3];
	PetscReal		element_vol = 0;
  
	PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],
                      e->ng,PetscReal,&dv_elem[1],
                      e->ng,PetscReal,&dv_elem[2],
                      e->ng,PetscReal,&u_elem[0],
                      e->ng,PetscReal,&u_elem[1],
                      e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&dv_mag_elem);CHKERRQ(ierr);
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
          //        n_elem[c][eg] = PetscAbs(dv_elem[c][eg])/dv_mag_elem[eg];
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
#define __FUNCT__ "VolumetricStrainVolume_local"
extern PetscErrorCode VolumetricStrainVolume_local(PetscReal *VolStrainVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, VFFlowProp flowpropty, PetscReal ****u_diff_array, PetscReal ***v_array)
{
  
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
	PetscReal    *v_elem,*du_elem[3],alphabiot,timestepsize;
  
	PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,PetscReal,&v_elem,
                      e->ng,PetscReal,&du_elem[0],
                      e->ng,PetscReal,&du_elem[1],
                      e->ng,PetscReal,&du_elem[2]);CHKERRQ(ierr);
  timestepsize     = flowpropty.timestepsize;
  alphabiot  = flowpropty.alphabiot;
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
      *VolStrainVolume_local += -1*alphabiot*(pow(v_elem[eg],2))*du_elem[c][eg]*e->weight[eg]/timestepsize;
    }
	}
	ierr = PetscFree4(du_elem[0],du_elem[1],du_elem[2],v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_ComputeAverageVField_local"
extern PetscErrorCode VF_ComputeAverageVField_local(PetscReal ***v_c_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e)
{
	PetscInt		i, j, k;
  
	PetscFunctionBegin;
  v_c_array[ek][ej][ei] = 0;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        v_c_array[ek][ej][ei] += v_array[ek+k][ej+j][ei+i]/(e->nphiz*e->nphiy*e->nphix);
        
      }
    }
  }
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "VF_ComputeAverageGradient"
extern PetscErrorCode VF_ComputeAverageGradient(PetscReal *grad_array, PetscInt dof, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
	PetscReal		*dv_elem[3],*v_mag_elem;
  PetscReal		element_vol = 0;
  
	PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&v_mag_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for (c = 0; c < 3; c++){
      dv_elem[c][eg] = 0.;
    }
    v_mag_elem[eg] = 0.;
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
    v_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    dv_elem[dof][eg] = dv_elem[dof][eg]/v_mag_elem[eg];
    if((ierr = PetscIsInfOrNanScalar(dv_elem[dof][eg])))
    {
      dv_elem[dof][eg] = 0;
    }
    dv_elem[dof][eg] = dv_elem[dof][eg];
  }
	*grad_array = 0.;
	for(eg = 0; eg < e->ng; eg++){
		element_vol += e->weight[eg];
		*grad_array += dv_elem[dof][eg]*e->weight[eg];
	}
  *grad_array = *grad_array/element_vol;
  ierr = PetscFree4(dv_elem[0],dv_elem[1],dv_elem[2],v_mag_elem);CHKERRQ(ierr);
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
  PetscReal       mymodVolumeLocal = 0.,mymodVolume = 0.;
  PetscReal       mydivVolumeLocal = 0.,mydivVolume = 0.;
  PetscReal       mysurfVolumeLocal = 0.,mysurfVolume = 0.;
  PetscReal       mysourceVolumeLocal = 0.,mysourceVolume = 0.;
  PetscReal       mystrainVolumeLocal = 0.,mystrainVolume = 0.;
  FACE           face;
  
  PetscFunctionBegin;
  
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
  ierr = VecAXPY(U_diff,1.0,ctx->U);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daVect,&u_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,U_diff,INSERT_VALUES,u_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,U_diff,INSERT_VALUES,u_diff_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_diff_local,&u_diff_array);CHKERRQ(ierr);
  
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
				ierr = ModulusVolume_local(&mymodVolumeLocal, ek, ej, ei, &ctx->e3D,ctx->flowprop,press_diff_array,v_array);CHKERRQ(ierr);
				ierr = DivergenceVolume_local(&mydivVolumeLocal, ek, ej, ei, &ctx->e3D,vel_array, v_array);CHKERRQ(ierr);
        if(ctx->hasFluidSources){
          ierr = SourceVolume_local(&mysourceVolumeLocal, ek, ej, ei, &ctx->e3D, src_array, v_array);CHKERRQ(ierr);
        }
        if(ctx->FlowDisplCoupling){
          ierr = VolumetricStrainVolume_local(&mystrainVolumeLocal, ek, ej, ei, &ctx->e3D,ctx->flowprop,u_diff_array,v_array);CHKERRQ(ierr);
        }
        mymodVolume += mymodVolumeLocal;
        mysourceVolume += mysourceVolumeLocal;
        mydivVolume += mydivVolumeLocal;
        mystrainVolume += mystrainVolumeLocal;
        if(ei == 0){
          face = X0;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += mysurfVolumeLocal;
        }
        if(ei == nx-1){
          face = X1;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += mysurfVolumeLocal;
        }
        if(ej == 0){
          face = Y0;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += mysurfVolumeLocal;
        }
        if(ej == ny-1){
          face = Y1;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += mysurfVolumeLocal;
        }
        if(ek == 0){
          face = Z0;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += mysurfVolumeLocal;
        }
        if(ek == nz-1){
          face = Z1;
          ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
          ierr = SurfaceFluxVolume_local(&mysurfVolumeLocal,ek,ej,ei,face,&ctx->e2D,vel_array,v_array);
          mysurfVolume += mysurfVolumeLocal;
        }
			}
		}
	}
  *SumWellRate = 0;
  for (i = 0; i < ctx->numWells; i++) {
    if(ctx->well[i].type == INJECTOR){
      *SumWellRate = *SumWellRate + ctx->well[i].Qw;
    }
    else{
      *SumWellRate = *SumWellRate - ctx->well[i].Qw;
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
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModulusVolume_local"
extern PetscErrorCode ModulusVolume_local(PetscReal *ModVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, VFFlowProp flowpropty, PetscReal ***press_diff_array,PetscReal ***v_array)
{
  
	PetscErrorCode ierr;
	PetscInt		i, j, k;
	PetscInt		eg;
	PetscReal		M_inv,timestepsize,*v_elem,*press_diff_elem;
  
	PetscFunctionBegin;
  
  ierr = PetscMalloc2(e->ng,PetscReal,&v_elem,e->ng,PetscReal,&press_diff_elem);CHKERRQ(ierr);
  M_inv     = flowpropty.M_inv;
  timestepsize     = flowpropty.timestepsize;
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
		*ModVolume_local += (pow(v_elem[eg],2))*M_inv*press_diff_elem[eg]*e->weight[eg]/timestepsize;
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
  ierr = PetscMalloc2(e->ng,PetscReal,&v_elem,e->ng,PetscReal,&flux_elem);CHKERRQ(ierr);
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
  ierr = PetscMalloc4(e->ng,PetscReal,&dvel_elem[0],e->ng,PetscReal,&dvel_elem[1],e->ng,PetscReal,&dvel_elem[2],e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
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
  ierr = PetscMalloc2(e->ng,PetscReal,&v_elem,e->ng,PetscReal,&src_elem);CHKERRQ(ierr);
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
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricCrackOpening3D_local(&myLeakOffRateLocal, volleakoffrate_array,PETSC_NULL, q_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
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
extern PetscErrorCode Permeabilityfield(PetscReal *COD_local, PetscReal ***volcrackopening_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement2D *e, FACE face)
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
/*
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
 */
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
	PetscInt        k, j, i;
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
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				ierr = VolumetricCrackOpening3D_local(&myCOD, cod_array, PETSC_NULL, u_array, v_array, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
			}
		}
	}
	PetscReal	vmult = 0;
	PetscReal  maxperm = 10;
	PetscReal  minperm = 0.01;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-perm",&minperm,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-permmax",&maxperm,PETSC_NULL);CHKERRQ(ierr);
  
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
        
				if (Ele_v_ave > 0.2){
            //                    for(c = 0; c < 2; c++)
            //                      perm_array[ek][ej][ei][c] = minperm;
          perm_array[ek][ej][ei][0] = minperm;
          perm_array[ek][ej][ei][2] = minperm;
				}
				else {
            //                  for(c = 0; c < 2; c++)
            //                    perm_array[ek][ej][ei][c] = maxperm;
          perm_array[ek][ej][ei][0] = maxperm;
          perm_array[ek][ej][ei][2] = maxperm;
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

#undef __FUNCT__
#define __FUNCT__ "VolumetricFractureWellRate_local"
extern PetscErrorCode VolumetricFractureWellRate_local(PetscReal *InjVolume_local, PetscReal ***regrate_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k,c;
	PetscInt		eg;
  PetscReal    *dv_elem[3],*regrate_elem;
  
	PetscFunctionBegin;
	ierr = PetscMalloc4(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&regrate_elem);CHKERRQ(ierr);
  
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
  ierr = PetscObjectQuery((PetscObject)vec,"DM",(PetscObject*)&dm);CHKERRQ(ierr);
	if (!dm) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
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

