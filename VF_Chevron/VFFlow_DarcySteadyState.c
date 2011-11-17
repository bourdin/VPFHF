/*
   VFFlow_DarcySteadyState.c
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFFlow_DarcySteadyState.h"

#undef __FUNCT__
#define __FUNCT__ "VFFlow_DarcySteadyState"

extern PetscErrorCode VFFlow_DarcySteadyState(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode				ierr;
	PetscViewer			viewer;
	PetscInt					xs,xm,ys,ym,zs,zm;
	PetscInt					i,j,k;
	PetscInt					its;
	KSPConvergedReason			reason;
	PetscReal					****VelnPress_array;
	PetscReal					***Press_array;
	PetscFunctionBegin;
	ierr = DAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = GetFlowProp(&ctx->flowprop, ctx->units, ctx->resprop);CHKERRQ(ierr);
	ierr = SETFlowBC(&ctx->bcFlow[0], ctx->flowcase);CHKERRQ(ierr);
	ierr = FlowMatnVecAssemble(ctx->KVelP, ctx->RHSVelP, fields, ctx);CHKERRQ(ierr);
	ierr = KSPSolve(ctx->kspVelP,ctx->RHSVelP,fields->VelnPress);CHKERRQ(ierr);
	ierr = KSPGetConvergedReason(ctx->kspVelP,&reason);CHKERRQ(ierr);
	if (reason < 0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] kspVelP diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
	} else {
		ierr = KSPGetIterationNumber(ctx->kspVelP,&its);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"      kspVelP converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
	}
	/*The next few lines equate the values of pressure calculated from the flow solver, to the pressure defined in da=daScal*/
	ierr = DAVecGetArray(ctx->daScal,fields->pressure,&Press_array);CHKERRQ(ierr);
	ierr = DAVecGetArrayDOF(ctx->daFlow,fields->VelnPress,&VelnPress_array);CHKERRQ(ierr);
	for (k = zs; k < zs + zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				Press_array[k][j][i] = VelnPress_array[k][j][i][3];
//				printf("p[%d][%d][%d] = %f\n",k, j, i, Press_array[k][j][i]);
			}
		}
	}
	ierr = DAVecRestoreArray(ctx->daScal,fields->pressure,&Press_array);CHKERRQ(ierr);
	ierr = DAVecRestoreArrayDOF(ctx->daFlow,fields->VelnPress,&VelnPress_array);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SETFlowBC"
/* 
 Fake flow solver for VF_Chevron.c test
 */
extern PetscErrorCode SETFlowBC(FLOWBC *BC, FlowCases flowcase)
{
	PetscInt		i, c;
	PetscFunctionBegin;
	for (i = 0; i < 6; i++) {
		for (c = 0; c < 4; c++) {
			BC[c].face[i] = NOBC;
		}
	}
	for (i = 0; i < 12; i++) {
		for (c = 0; c < 4; c++) {
			BC[c].edge[i] = NOBC;
		}
	}
	for (i = 0; i < 8; i++) {
		for (c = 0; c < 4; c++) {
			BC[c].vertex[i] = NOBC;
		}
	}
	switch (flowcase) {
		case ALLPRESSUREBC:
			for (i = 0; i < 6; i++) {
				BC[3].face[i] = PRESSURE;
			}
			for (i = 0; i < 12; i++) {
				BC[3].edge[i] = PRESSURE;
			}
			for (i = 0; i < 8; i++) {
				BC[3].vertex[i] = PRESSURE;
			}
			break;
		case ALLNORMALFLOWBC:
			for (i = 0; i < 6; i++) {
				for (c = 0; c < 3; c++) {
					BC[c].face[i] = VELOCITY;
				}
			}
			for (i = 0; i < 12; i++) {
				for (c = 0; c < 3; c++) {
					BC[c].edge[i] = VELOCITY;
				}
			}
			for (i = 0; i < 8; i++) {
				for (c = 0; c < 3; c++) {
					BC[c].vertex[i] = VELOCITY;
				}
			}
			break;
		default:
			SETERRQ2(PETSC_ERR_USER,"ERROR: [%s] unknown FLOWCASE %i .\n",__FUNCT__,flowcase);
			break;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplyWellFlowBC"
/* 
 Fake flow solver for VF_Chevron.c test
 */
extern PetscErrorCode VecApplyWellFlowBC(Vec RHS, VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt		xs,xm,nx;
	PetscInt		ys,ym,ny;
	PetscInt		zs,zm,nz;
	PetscInt		i,j,k;
	PetscInt		dim, dof;
	PetscReal		****coords_array;
	PetscReal		****RHS_array;
	PetscReal		alpha_c;
	DA				da;
	PetscReal		loc_src[3];					/*source location*/
	PetscReal		loc_sink[3];				/*sink location*/
	PetscReal		flowrate;
	
						/*we'd assume for now that the source and sinks are located at the following nodes co
						*/	
	PetscFunctionBegin;
	ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

	flowrate = ctx->flowrate;
	alpha_c = ctx->flowprop.alpha;

	
	ierr = PetscObjectQuery((PetscObject) RHS,"DA",(PetscObject *) &da); CHKERRQ(ierr);
    if (!da) SETERRQ(PETSC_ERR_ARG_WRONG,"Vector not generated from a DA");
	
	ierr = DAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,&dof,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DAVecGetArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
	
	loc_src[0] = coords_array[1][1][1][0];
	loc_src[1] = coords_array[1][1][1][1];
	loc_src[2] = coords_array[1][1][1][2];
	
	loc_sink[0] = coords_array[nz-2][ny-2][nx-2][0];
	loc_sink[1] = coords_array[nz-2][ny-2][nx-2][1];
	loc_sink[2] = coords_array[nz-2][ny-2][nx-2][2];
	
	for (k = zs; k < zs + zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				if(coords_array[k][j][i][0] == loc_src[0] && coords_array[k][j][i][1] == loc_src[1] && coords_array[k][j][i][2] == loc_src[2])
				{
					RHS_array[k][j][i][3] += -flowrate/alpha_c;
				}
				if(coords_array[k][j][i][0] == loc_sink[0] && coords_array[k][j][i][1] == loc_sink[1] && coords_array[k][j][i][2] == loc_sink[2])
				{
					RHS_array[k][j][i][3] += flowrate/alpha_c;
				}
			}
		}
	}
	ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DAVecRestoreArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplyFlowBC"
extern PetscErrorCode VecApplyFlowBC(Vec RHS, FLOWBC *BC, ResProp resprop)
{
	PetscErrorCode ierr;
	PetscInt		xs,xm,nx;
	PetscInt		ys,ym,ny;
	PetscInt		zs,zm,nz;
	PetscInt		dim, dof;
	PetscInt		i,j,k,c;
	PetscReal		PressFaceBC_Values[6];
	PetscReal		PressEdgeBC_Values[8];
	PetscReal		PressVertexBC_Values[12];
	DA				da;
	PetscReal		VelFaceBC_Values[6][3];
	PetscReal		VelEdgeBC_Values[12][3];
	PetscReal		VelVertexBC_Values[8][3];
	PetscReal		****RHS_array;
	
	PetscFunctionBegin;
	for(i = 0; i < 6; i++) PressFaceBC_Values[i] = 0.;//resprop.Pinit;
	for(i = 0; i < 8; i++) PressEdgeBC_Values[i] = 0.;//resprop.Pinit;
	for(i = 0; i < 12; i++) PressVertexBC_Values[i] = 0.;//resprop.Pinit;

	ierr = PetscObjectQuery((PetscObject) RHS,"DA",(PetscObject *) &da); CHKERRQ(ierr);
    if (!da) SETERRQ(PETSC_ERR_ARG_WRONG,"Vector not generated from a DA");
	
	ierr = DAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,&dof,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DAVecGetArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);

	for(c = 0; c < dof; c++){
		if(xs == 0){
			i = 0;
			//Realignment of partitioning to avoid y & z direction edges 
			if (ys+ym == ny) ym--;
			if (zs+zm == nz) zm--;
			if (ys == 0) ys++;
			if (zs == 0) zs++;
			if(BC[c].face[X0] == PRESSURE){
				for(k = zs; k < zs+zm; k++){
					for(j = ys; j < ys+ym; j++){
						RHS_array[k][j][i][c] = PressFaceBC_Values[0];
					}
				}
			}
			if(BC[c].face[X0] == VELOCITY){
				for(k = zs; k < zs+zm; k++){
					for(j = ys; j < ys+ym; j++){
						RHS_array[k][j][i][c] = VelFaceBC_Values[0][c];
					}
				}
			}
			if (ys == 1) ys--;
			if (zs == 1) zs--;
			if (ys+ym == ny-1) ym++;
			if (zs+zm == nz-1) zm++;
		}
		if(xs+xm == nx){
			i = nx-1;
			//Realignment of partitioning to avoid y & z direction edges 
			if (ys+ym == ny) ym--;
			if (zs+zm == nz) zm--;
			if (ys == 0) ys++;
			if (zs == 0) zs++;	
			if(BC[c].face[X1] == PRESSURE){
				for(k = zs; k < zs+zm; k++){
					for(j = ys; j < ys+ym; j++){
						RHS_array[k][j][i][c] = PressFaceBC_Values[1];
					}
				}
			}
			if(BC[c].face[X1] == VELOCITY){
				for(k = zs; k < zs+zm; k++){
					for(j = ys; j < ys+ym; j++){
						RHS_array[k][j][i][c] = VelFaceBC_Values[1][c];
					}
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;
			if (zs == 1) zs--;
			if (ys+ym == ny-1) ym++;
			if (zs+zm == nz-1) zm++;	
		}
		if(ys == 0){
			j = 0;
			//Realignment of partitioning to avoid y & z direction edges 
			if (xs+xm == nx) xm--;
			if (zs+zm == nz) zm--;
			if (xs == 0) xs++;
			if (zs == 0) zs++;	
			if(BC[c].face[Y0] == PRESSURE){
				for(k = zs; k < zs+zm; k++){
					for(i = xs; i < xs+xm; i++){
						RHS_array[k][j][i][c] = PressFaceBC_Values[2];
					}
				}
			}
			if(BC[c].face[Y0] == VELOCITY){
				for(k = zs; k < zs+zm; k++){
					for(i = xs; i < xs+xm; i++){
						RHS_array[k][j][i][c] = VelFaceBC_Values[2][c];
					}
				}
			}
			if (xs == 1) xs--;
			if (zs == 1) zs--;
			if (xs+xm == nx-1) xm++;
			if (zs+zm == nz-1) zm++;
		}
		if(ys+ym == ny){
			j = ny-1;
			//Realignment of partitioning to avoid y & z direction edges 
			if (xs+xm == nx) xm--;
			if (zs+zm == nz) zm--;
			if (xs == 0) xs++;
			if (zs == 0) zs++;	
			if(BC[c].face[Y1] == PRESSURE){
				for(k = zs; k < zs+zm; k++){
					for(i = xs; i < xs+xm; i++){
						RHS_array[k][j][i][c] = PressFaceBC_Values[3];
					}
				}
			}
			if(BC[c].face[Y1] == VELOCITY){
				for(k = zs; k < zs+zm; k++){
					for(i = xs; i < xs+xm; i++){
						RHS_array[k][j][i][c] = VelFaceBC_Values[3][c];
					}
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;
			if (zs == 1) zs--;
			if (xs+xm == nx-1) xm++;
			if (zs+zm == nz-1) zm++;	
		}	
		if(zs == 0){
			k = 0;
			//Realignment of partitioning to avoid y & z direction edges 
			if (xs+xm == nx) xm--;
			if (ys+ym == ny) ym--;
			if (xs == 0) xs++;
			if (ys == 0) ys++;	
			if(BC[c].face[Z0] == PRESSURE){
				for(j = ys; j < ys+ym; j++){
					for(i = xs; i < xs+xm; i++){
						RHS_array[k][j][i][c] = PressFaceBC_Values[4];
					}
				}
			}
			if(BC[c].face[Z0] == VELOCITY){
				for(j = ys; j < ys+ym; j++){
					for(i = xs; i < xs+xm; i++){
						RHS_array[k][j][i][c] = VelFaceBC_Values[4][c];
					}
				}
			}
			if (xs == 1) xs--;
			if (ys == 1) ys--;
			if (xs+xm == nx-1) xm++;
			if (ys+ym == ny-1) ym++;
		}
		if(zs+zm == nz){
			k = nz-1;
			//Realignment of partitioning to avoid y & z direction edges 
			if (xs+xm == nx) xm--;
			if (ys+ym == ny) ym--;
			if (xs == 0) xs++;
			if (ys == 0) ys++;	
			if(BC[c].face[Z1] == PRESSURE){
				for(j = ys; j < ys+ym; j++){
					for(i = xs; i < xs+xm; i++){
						RHS_array[k][j][i][c] = PressFaceBC_Values[5];
					}
				}
			}
			if(BC[c].face[Z1] == VELOCITY){
				for(j = ys; j < ys+ym; j++){
					for(i = xs; i < xs+xm; i++){
						RHS_array[k][j][i][c] = VelFaceBC_Values[5][c];
					}
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;
			if (ys == 1) ys--;
			if (xs+xm == nx-1) xm++;
			if (ys+ym == ny-1) ym++;	
		}
		if (xs == 0 && zs == 0) {
			k = 0; i = 0;
			//Realignment of partitioning to avoid corner 
			if (ys+ym == ny) ym--;	
			if (ys == 0) ys++;	
			if(BC[c].edge[X0Z0] == PRESSURE){
				for(j = ys; j < ys+ym; j++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[0];
				}
			}
			if(BC[c].edge[X0Z0] == VELOCITY){
				for(j = ys; j < ys+ym; j++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[0][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
				if (ys == 1) ys--;	
				if (ys+ym == ny-1) ym++;
		}
		if (xs+xm == nx && zs == 0) {
			k = 0; i = nx-1;
			//Realignment of partitioning to avoid corner 
			if (ys+ym == ny) ym--;	
			if (ys == 0) ys++;	
			if(BC[c].edge[X1Z0] == PRESSURE){
				for(j = ys; j < ys+ym; j++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[1];
				}
			}
			if(BC[c].edge[X1Z0] == VELOCITY){
				for(j = ys; j < ys+ym; j++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[1][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;	
			if (ys+ym == ny-1) ym++;
		}
		if (ys == 0 && zs == 0) {
			k = 0; j = 0;
			//Realignment of partitioning to avoid corner 
			if (xs+xm == nx) xm--;	
			if (xs == 0) xs++;	
			if(BC[c].edge[Y0Z0] == PRESSURE){
				for(i = xs; i < xs+xm; i++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[2];
				}
			}
			if(BC[c].edge[Y0Z0] == VELOCITY){
				for(i = xs; i < xs+xm; i++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[2][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;	
			if (xs+xm == nx-1) xm++;
		}
		if (ys+ym == ny && zs == 0) {
			k = 0; j = 0;
			//Realignment of partitioning to avoid corner 
			if (xs+xm == nx) xm--;	
			if (xs == 0) xs++;	
			if(BC[c].edge[Y1Z0] == PRESSURE){
				for(i = xs; i < xs+xm; i++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[3];
				}
			}
			if(BC[c].edge[Y1Z0] == VELOCITY){
				for(i = xs; i < xs+xm; i++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[3][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;	
			if (xs+xm == nx-1) xm++;
		}
		if (xs == 0 && zs+zm == nz) {
			k = nz-1; i = 0;
			//Realignment of partitioning to avoid corner 
			if (ys+ym == ny) ym--;	
			if (ys == 0) ys++;	
			if(BC[c].edge[X0Z1] == PRESSURE){
				for(j = ys; j < ys+ym; j++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[4];
				}
			}
			if(BC[c].edge[X0Z1] == VELOCITY){
				for(j = ys; j < ys+ym; j++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[4][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;	
			if (ys+ym == ny-1) ym++;
		}
		if (xs+xm == nx && zs+zm == nz) {
			k = nz-1; i = nx-1;
			//Realignment of partitioning to avoid corner 
			if (ys+ym == ny) ym--;	
			if (ys == 0) ys++;	
			if(BC[c].edge[X1Z1] == PRESSURE){
				for(j = ys; j < ys+ym; j++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[5];
				}
			}
			if(BC[c].edge[X1Z1] == VELOCITY){
				for(j = ys; j < ys+ym; j++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[5][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;	
			if (ys+ym == ny-1) ym++;
		}		
		if (ys == 0 && zs+zm == nz) {
			k = nz-1; j = 0;
			//Realignment of partitioning to avoid corner 
			if (xs+xm == nx) xm--;	
			if (xs == 0) xs++;	
			if(BC[c].edge[Y0Z1] == PRESSURE){
				for(i = xs; i < xs+xm; i++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[6];
				}
			}
			if(BC[c].edge[Y0Z1] == VELOCITY){
				for(i = xs; i < xs+xm; i++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[6][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;	
			if (xs+xm == nx-1) xm++;
		}
		if (ys+ym == ny && zs+zm == nz) {
			k = nz-1; j = ny-1;
			//Realignment of partitioning to avoid corner 
			if (xs+xm == nx) xm--;	
			if (xs == 0) xs++;	
			if(BC[c].edge[Y1Z1] == PRESSURE){
				for(i = xs; i < xs+xm; i++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[7];
				}
			}
			if(BC[c].edge[Y1Z1] == VELOCITY){
				for(i = xs; i < xs+xm; i++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[7][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;	
			if (xs+xm == nx-1) xm++;
		}
		if (xs == 0 && ys == 0) {
			j = 0; i = 0;
			//Realignment of partitioning to avoid corner 
			if (zs+zm == nz) zm--;	
			if (zs == 0) zs++;	
			if(BC[c].edge[X0Y0] == PRESSURE){
				for(k = zs; k < zs+zm; k++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[8];
				}
			}
			if(BC[c].edge[X0Y0] == VELOCITY){
				for(k = zs; k < zs+zm; k++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[8][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (zs == 1) zs--;	
			if (zs+zm == nz-1) zm++;
		}		
		if (xs == 0 && ys+ym == ny) {
			j = ny-1; i = 0;
			//Realignment of partitioning to avoid corner 
			if (zs+zm == nz) zm--;	
			if (zs == 0) zs++;	
			if(BC[c].edge[X0Y1] == PRESSURE){
				for(k = zs; k < zs+zm; k++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[9];
				}
			}
			if(BC[c].edge[X0Y1] == VELOCITY){
				for(k = zs; k < zs+zm; k++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[9][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (zs == 1) zs--;	
			if (zs+zm == nz) zm++;
		}		
		if (xs+xm == nx && ys == 0) {
			j = 0; i = nx-1;
			//Realignment of partitioning to avoid corner 
			if (zs+zm == nz) zm--;	
			if (zs == 0) zs++;	
			if(BC[c].edge[X1Y0] == PRESSURE){
				for(k = zs; k < zs+zm; k++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[10];
				}
			}
			if(BC[c].edge[X1Y0] == VELOCITY){
				for(k = zs; k < zs+zm; k++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[10][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (zs == 1) zs--;	
			if (zs+zm == nz-1) zm++;
		}		
		if (xs+xm == nx && ys+ym == ny) {
			j = ny-1; i = nx-1;
			//Realignment of partitioning to avoid corner 
			if (zs+zm == nz) zm--;	
			if (zs == 0) zs++;	
			if(BC[c].edge[X1Y1] == PRESSURE){
				for(k = zs; k < zs+zm; k++){
					RHS_array[k][j][i][c] = PressEdgeBC_Values[11];
				}
			}
			if(BC[c].edge[X1Y1] == VELOCITY){
				for(k = zs; k < zs+zm; k++){
					RHS_array[k][j][i][c] = VelEdgeBC_Values[11][c];
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (zs == 1) zs--;	
			if (zs+zm == nz-1) zm++;
		}		
		if (xs == 0 && ys == 0 && zs == 0) {
			k = 0; j = 0; i = 0;
			if(BC[c].edge[X0Y0Z0] == PRESSURE){
				RHS_array[k][j][i][c] = PressVertexBC_Values[0];
			}
			if(BC[c].edge[X0Y0Z0] == VELOCITY){
				RHS_array[k][j][i][c] = VelVertexBC_Values[0][c];
			}
		}		
		if (xs+xm == nx && ys == 0 && zs == 0) {
			k = 0; j = 0; i = nx-1;
			if(BC[c].edge[X1Y0Z0] == PRESSURE){
				RHS_array[k][j][i][c] = PressVertexBC_Values[1];
			}
			if(BC[c].edge[X1Y0Z0] == VELOCITY){
				RHS_array[k][j][i][c] = VelVertexBC_Values[1][c];
			}
		}		
		if (xs == 0 && ys+ym == ny && zs == 0) {
			k = 0; j = ny-1; i = 0;
			if(BC[c].edge[X0Y1Z0] == PRESSURE){
				RHS_array[k][j][i][c] = PressVertexBC_Values[2];
			}
			if(BC[c].edge[X0Y1Z0] == VELOCITY){
				RHS_array[k][j][i][c] = VelVertexBC_Values[2][c];
			}
		}			
		if (xs+xm == nx && ys+ym == ny && zs == 0) {
			k = 0; j = ny-1; i = nx-1;
			if(BC[c].edge[X1Y1Z0] == PRESSURE){
				RHS_array[k][j][i][c] = PressVertexBC_Values[3];
			}
			if(BC[c].edge[X1Y1Z0] == VELOCITY){
				RHS_array[k][j][i][c] = VelVertexBC_Values[3][c];
			}
		}	
		if (xs == 0 && ys == 0 && zs+zm == nz) {
			k = nz-1; j = 0; i = 0;
			if(BC[c].edge[X0Y0Z1] == PRESSURE){
				RHS_array[k][j][i][c] = PressVertexBC_Values[4];
			}
			if(BC[c].edge[X0Y0Z1] == VELOCITY){
				RHS_array[k][j][i][c] = VelVertexBC_Values[4][c];
			}
		}	
		if (xs+xm == nx && ys == 0 && zs+zm == nz) {
			k = nz-1; j = 0; i = nx-1;
			if(BC[c].edge[X1Y0Z1] == PRESSURE){
				RHS_array[k][j][i][c] = PressVertexBC_Values[5];
			}
			if(BC[c].edge[X1Y0Z1] == VELOCITY){
				RHS_array[k][j][i][c] = VelVertexBC_Values[5][c];
			}
		}	
		if (xs == 0 && ys+ym == ny && zs+zm == nz) {
			k = nz-1; j = ny-1; i = 0;
			if(BC[c].edge[X0Y1Z1] == PRESSURE){
				RHS_array[k][j][i][c] = PressVertexBC_Values[6];
			}
			if(BC[c].edge[X0Y1Z1] == VELOCITY){
				RHS_array[k][j][i][c] = VelVertexBC_Values[6][c];
			}
		}			
		if (xs+xm == nx && ys+ym == ny && zs+zm == nz) {
			k = nz-1; j = ny-1; i = nx-1;
			if(BC[c].edge[X1Y1Z1] == PRESSURE){
				RHS_array[k][j][i][c] = PressVertexBC_Values[7];
			}
			if(BC[c].edge[X1Y1Z1] == VELOCITY){
				RHS_array[k][j][i][c] = VelVertexBC_Values[7][c];
			}
		}			
	}
	ierr = DAVecRestoreArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatApplyFlowBC"
extern PetscErrorCode MatApplyFlowBC(Mat K, DA da, FLOWBC *BC)
{
	PetscErrorCode ierr;
	PetscInt		xs,xm,nx;
	PetscInt		ys,ym,ny;
	PetscInt		zs,zm,nz;
	PetscInt		dim, dof;
	PetscReal		one = 1.;
	PetscInt		i,j,k,c;
	MatStencil     *row;
	PetscFunctionBegin;
	
	ierr = DAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,&dof,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = PetscMalloc(1 * sizeof(MatStencil),&row);CHKERRQ(ierr);	
	
	for(c = 0; c < dof; c++){
		//Realignment of partitioning to avoid y & z direction edges 
//1.	
		if(xs == 0 && BC[c].face[X0] != NOBC){
			if (ys+ym == ny) ym--;
			if (zs+zm == nz) zm--;
			if (ys == 0) ys++;
			if (zs == 0) zs++;
			for(k = zs; k < zs+zm; k++){
				for(j = ys; j < ys+ym; j++){
					row[0].i = 0; row[0].j = j; row[0].k = k; row[0].c = c; 
					ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;
			if (zs == 1) zs--;
			if (ys+ym == ny-1) ym++;
			if (zs+zm == nz-1) zm++;
		}
//2.
		if(xs + xm == nx && BC[c].face[X1] != NOBC){
			if (ys+ym == ny) ym--;
			if (zs+zm == nz) zm--;
			if (ys == 0) ys++;
			if (zs == 0) zs++;
			for(k = zs; k < zs+zm; k++){
				for(j = ys; j < ys+ym; j++){
					row[0].i = nx-1; row[0].j = j; row[0].k = k; row[0].c = c; 
					ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
				}
			}
		//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;
			if (zs == 1) zs--;
			if (ys+ym == ny-1) ym++;
			if (zs+zm == nz-1) zm++;
		}
//3
		if(ys == 0 && BC[c].face[Y0] != NOBC){
			//Realignment of partitioning to avoid x & z direction edges 
			if (xs+xm == nx) xm--;
			if (zs+zm == nz) zm--;
			if (xs == 0) xs++;
			if (zs == 0) zs++;	
			for(k = zs; k < zs+zm; k++){
				for(i = xs; i < xs+xm; i++){
					row[0].i = i; row[0].j = 0; row[0].k = k; row[0].c = c; 
					ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;
			if (zs == 1) zs--;	
			if (xs+xm == nx-1) xm++;
			if (zs+zm == nz-1) zm++;
		}
//4
		if(ys + ym == ny && BC[c].face[Y1] != NOBC){
			//Realignment of partitioning to avoid x & z direction edges 
			if (xs+xm == nx) xm--;
			if (zs+zm == nz) zm--;
			if (xs == 0) xs++;
			if (zs == 0) zs++;	
			for(k = zs; k < zs+zm; k++){
				for(i = xs; i < xs+xm; i++){
					row[0].i = i; row[0].j = ny-1; row[0].k = k; row[0].c = c; 
					ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;
			if (zs == 1) zs--;
			if (xs+xm == nx-1) xm++;
			if (zs+zm == nz-1) zm++;	
		}
//5
		if(zs == 0 && BC[c].face[Z0] != NOBC){
			//Realignment of partitioning to avoid x & z direction edges 
			if (ys+ym == ny) ym--;
			if (xs+xm == nx) xm--;
			if (ys == 0) ys++;
			if (xs == 0) xs++;	
			for(j = ys; j < ys+ym; j++){
				for(i = xs; i < xs+xm; i++){
					row[0].i = i; row[0].j = j; row[0].k = 0; row[0].c = c; 
					ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;
			if (xs == 1) xs--;	
			if (ys+ym == ny-1) ym++;
			if (xs+xm == nx-1) xm++;
		}
//6
		if(zs + zm == nz && BC[c].face[Z1] != NOBC){
			//Realignment of partitioning to avoid x & z direction edges 
			if (ys+ym == ny) ym--;
			if (xs+xm == nx) xm--;
			if (ys == 0) ys++;
			if (xs == 0) xs++;	
			for(j = ys; j < ys+ym; j++){
				for(i = xs; i < xs+xm; i++){
					row[0].i = i; row[0].j = j; row[0].k = nz-1; row[0].c = c; 
					ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
				}
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;
			if (xs == 1) xs--;
			if (ys+ym == ny-1) ym++;
			if (xs+xm == nx-1) xm++;	
		}
//1.		
		if (xs == 0 && zs == 0 && BC[c].edge[X0Z0] != NOBC) { 
			//Realignment of partitioning to avoid corners 
			if (ys+ym == ny) ym--;	
			if (ys == 0) ys++;	
			for(j = ys; j < ys+ym; j++){
				row[0].i = 0; row[0].j = j; row[0].k = 0; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;	
			if (ys+ym == ny-1) ym++;
		}
//2.		
		if (xs+xm == nx && zs == 0  && BC[c].edge[X1Z0] != NOBC) { 
			//Realignment of partitioning to avoid corners 
			if (ys+ym == ny) ym--;	
			if (ys == 0) ys++;	
			for(j = ys; j < ys+ym; j++){
				row[0].i = nx-1; row[0].j = j; row[0].k = 0; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;	
			if (ys+ym == ny-1) ym++;
		}		
//3.		Realignment of partitioning to avoid corner 
		if (ys == 0 && zs == 0 && BC[c].edge[Y0Z0] != NOBC) { 
			if (xs+xm == nx) xm--;	
			if (xs == 0) xs++;	
			for(i = xs; i < xs+xm; i++){
				row[0].i = i; row[0].j = 0; row[0].k = 0; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;	
			if (xs+xm == nx-1) xm++;
		}		
//4.		Realignment of partitioning to avoid corner 
		if (ys+ym == ny && zs == 0 && BC[c].edge[Y1Z0] != NOBC) { 
			if (xs+xm == nx) xm--;	
			if (xs == 0) xs++;	
			for(i = xs; i < xs+xm; i++){
				row[0].i = i; row[0].j = ny-1; row[0].k = 0; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;	
			if (xs+xm == nx-1) xm++;
		}		
//5.		
		if (xs == 0 && zs+zm== nz  && BC[c].edge[X0Z1] != NOBC) { 
			//Realignment of partitioning to avoid corners 
			if (ys+ym == ny) ym--;	
			if (ys == 0) ys++;	
			for(j = ys; j < ys+ym; j++){
				row[0].i = 0; row[0].j = j; row[0].k = nz-1; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;	
			if (ys+ym == ny-1) ym++;
		}		
//6.		
		if (xs+xm == nx && zs+zm== nz  && BC[c].edge[X1Z1] != NOBC) { 
			//Realignment of partitioning to avoid corners 
			if (ys+ym == ny) ym--;	
			if (ys == 0) ys++;	
			for(j = ys; j < ys+ym; j++){
				row[0].i = nx-1; row[0].j = j; row[0].k = nz-1; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (ys == 1) ys--;	
			if (ys+ym == ny-1) ym++;
		}
//7.		
		if (ys == 0 && zs+zm == nz && BC[c].edge[Y0Z1] != NOBC) { 
			if (xs+xm == nx) xm--;	
			if (xs == 0) xs++;	
			for(i = xs; i < xs+xm; i++){
				row[0].i = i; row[0].j = 0; row[0].k = nz-1; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;	
			if (xs+xm == nx-1) xm++;
		}
//8.		
		if (ys+ym == ny && zs+zm == nz  && BC[c].edge[Y1Z1] != NOBC) { 
			if (xs+xm == nx) xm--;	
			if (xs == 0) xs++;	
			for(i = xs; i < xs+xm; i++){
				row[0].i = i; row[0].j = ny-1; row[0].k = nz-1; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (xs == 1) xs--;	
			if (xs+xm == nx-1) xm++;
		}	
//9.	
		if (xs == 0 && ys == 0  && BC[c].edge[X0Y0] != NOBC) {
			//Realignment of partitioning to avoid corner 
			if (zs+zm == nz) zm--;	
			if (zs == 0) zs++;	
			for(k = zs; k < zs+zm; k++){
				row[0].i = 0; row[0].j = 0; row[0].k = k; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (zs == 1) zs--;	
			if (zs+zm == nz-1) zm++;
		}
//10.
		if (xs == 0 && ys+ym== ny  && BC[c].edge[X0Y1] != NOBC) { 
			//Realignment of partitioning to avoid corner 
			if (zs+zm == nz) zm--;	
			if (zs == 0) zs++;	
			for(k = zs; k < zs+zm; k++){
				row[0].i = 0; row[0].j = ny-1; row[0].k = k; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (zs == 1) zs--;	
			if (zs+zm == nz-1) zm++;
		}
//11.		
		if (xs+xm == nx && ys == 0  && BC[c].edge[X1Y0] != NOBC) { 
			//Realignment of partitioning to avoid corner 
			if (zs+zm == nz) zm--;	
			if (zs == 0) zs++;	
			for(k = zs; k < zs+zm; k++){
				row[0].i = nx-1; row[0].j = 0; row[0].k = k; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (zs == 1) zs--;	
			if (zs+zm == nz-1) zm++;
		}
//12.		
		if (xs+xm == nx && ys+ym== ny && BC[c].edge[X1Y1] != NOBC) { 
			//Realignment of partitioning to avoid corner 
			if (zs+zm == nz) zm--;	
			if (zs == 0) zs++;	
			for(k = zs; k < zs+zm; k++){
				row[0].i = nx-1; row[0].j = ny-1; row[0].k = k; row[0].c = c; 
				ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
			}
			//Alignment to proper partitioning preparatory for realignment
			if (zs == 1) zs--;	
			if (zs+zm == nz-1) zm++;
		}
		if(xs == 0 && ys == 0 && zs == 0 && BC[c].vertex[X0Y0Z0] != NOBC){
			row[0].i = 0; row[0].j = 0; row[0].k = 0; row[0].c = c; 
			ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
		}
		if(xs == 0 && ys == 0 && zs+zm == nz && BC[c].vertex[X0Y0Z1] != NOBC){
			row[0].i = 0; row[0].j = 0; row[0].k = nz-1; row[0].c = c; 
			ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
		}
		if(xs == 0 && ys+ym == ny && zs == 0 && BC[c].vertex[X0Y1Z0] != NOBC){
			row[0].i = 0; row[0].j = ny-1; row[0].k = 0; row[0].c = c; 
			ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
		}
		if(xs == 0 && ys+ym == ny && zs+zm == nz && BC[c].vertex[X0Y1Z1] != NOBC){
			row[0].i = 0; row[0].j = ny-1; row[0].k = nz-1; row[0].c = c; 
			ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
		}				
		if(xs+xm == nx && ys == 0 && zs == 0 && BC[c].vertex[X1Y0Z0] != NOBC){
			row[0].i = nx-1; row[0].j = 0; row[0].k = 0; row[0].c = c; 
			ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
		}
		if(xs+xm == nx && ys == 0 && zs+zm == nz && BC[c].vertex[X1Y0Z1] != NOBC){
			row[0].i = nx-1; row[0].j = 0; row[0].k = nz-1; row[0].c = c; 
			ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
		}
		if(xs+xm == nx && ys+ym == ny && zs == 0 && BC[c].vertex[X1Y1Z0] != NOBC){
			row[0].i = nx-1; row[0].j = ny-1; row[0].k = 0; row[0].c = c; 
			ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
		}
		if(xs+xm == nx && ys+ym == ny && zs+zm == nz && BC[c].vertex[X1Y1Z1] != NOBC){
			row[0].i = nx-1; row[0].j = ny-1; row[0].k = nz-1; row[0].c = c; 
			ierr = MatZeroRowsStencil(K,1,row,one);CHKERRQ(ierr);
		}		
	}
	
	ierr = PetscFree(row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowMatnVecAssemble"
extern PetscErrorCode FlowMatnVecAssemble(Mat K, Vec RHS, VFFields * fields, VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt		xs,xm,nx;
	PetscInt		ys,ym,ny;
	PetscInt		zs,zm,nz;
	PetscInt		ek, ej, ei;
	PetscInt		i, j, k, l;
	PetscInt		veldof = 3;
	PetscInt		c;
//	PetscReal		****perm_array;
	PetscReal		****coords_array;
	PetscReal		****RHS_array;
	PetscReal		*RHS_local;
	Vec				RHS_localVec;
//	Vec				perm_local;
	PetscReal		hx,hy,hz;
	PetscReal		*KA_local, *KB_local;
	PetscReal		beta_c, mu, gx, gy, gz;
	PetscInt		nrow = ctx->e3D.nphix * ctx->e3D.nphiy * ctx->e3D.nphiz;
	MatStencil		*row, *row1;
	PetscFunctionBegin;
	beta_c = ctx->flowprop.beta;
	mu = ctx->flowprop.mu;
	gx = ctx->flowprop.g[0];
	gy = ctx->flowprop.g[1];
	gz = ctx->flowprop.g[2];
	
//	printf("/n#################### %f = ################/n", beta_c/mu);
	
	ierr = DAGetInfo(ctx->daFlow,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	// This line ensures that the number of cells is one less than the number of nodes. Force processing of cells to stop once the second to the last node is processed
	if (xs+xm == nx) xm--;
	if (ys+ym == ny) ym--;
	if (zs+zm == nz) zm--;
	
	ierr = MatZeroEntries(K);CHKERRQ(ierr);
	ierr = VecSet(RHS,0.);CHKERRQ(ierr);
	
	//Get coordinates from daVect since ctx->coordinates was created as an object in daVect
	ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DAGetLocalVector(ctx->daFlow,&RHS_localVec);CHKERRQ(ierr);
	ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
	ierr = DAVecGetArrayDOF(ctx->daFlow, RHS_localVec,&RHS_array);CHKERRQ(ierr);    

	/*
	 get permeability values/mult
	 */
	/*
	ierr = DAGetLocalVector(ctx->daFlow,&perm_local);CHKERRQ(ierr);
	ierr = DAGlobalToLocalBegin(ctx->daFlow,fields->pmult,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DAGlobalToLocalEnd(ctx->daFlow,fields->pmult,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DAVecGetArray(ctx->daFlow,perm_local,&perm_array);CHKERRQ(ierr);   
	*/
	
	/*
	 get local mat and RHS
	 */
	ierr = PetscMalloc(nrow * nrow * sizeof(PetscReal),&KA_local);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow * nrow * sizeof(PetscReal),&KB_local);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow * sizeof(PetscReal),&RHS_local);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow * sizeof(MatStencil),&row);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow * sizeof(MatStencil),&row1);CHKERRQ(ierr);
	
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				/*This computes the local contribution of the global A matrix*/
				ierr = FLow_MatA(KA_local, &ctx->e3D, ek, ej, ei);CHKERRQ(ierr);
				for(c = 0; c < veldof; c++){
					/*This computes the local contribution of the global B matrix for each direction i.e x, y and z*/
					ierr = FLow_MatB(KB_local, &ctx->e3D, ek, ej, ei, c);CHKERRQ(ierr);
					for(l = 0, k = 0; k < ctx->e3D.nphiz; k++){
						for(j = 0; j < ctx->e3D.nphiy; j++){
							for(i = 0; i < ctx->e3D.nphix; i++, l++){
								row[l].i = ei+i; row[l].j = ej+j; row[l].k = ek+k; row[l].c = c;
								row1[l].i = ei+i; row1[l].j = ej+j; row1[l].k = ek+k; row1[l].c = 3;								
							}
						}
					}
					ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KA_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(K,nrow,row1,nrow,row,KB_local,ADD_VALUES);CHKERRQ(ierr);
					for(l = 0; l < nrow*nrow;l++) KB_local[l] = KB_local[l] * beta_c/mu;												/*proper value for the transpose. Still need to multiply with the permeability when available*/
			//		for(l = 0; l < nrow*nrow;l++) KB_local[l] = KB_local[l]*1.127;//*beta_c/mu;								/*proper value for the transpose. Still need to multiply with the permeability when available*/
					ierr = MatSetValuesStencil(K,nrow,row,nrow,row1,KB_local,ADD_VALUES);CHKERRQ(ierr);						/* Adding transpose of KBx_local*/
				}
				/*Used same KA_local to access the local components of the D global matrix*/
				ierr = FLow_MatD(KA_local, &ctx->e3D, ek, ej, ei, ctx->flowprop);CHKERRQ(ierr);
				for(l = 0, k = 0; k < ctx->e3D.nphiz; k++){
					for(j = 0; j < ctx->e3D.nphiy; j++){
						for(i = 0; i < ctx->e3D.nphix; i++, l++){
							row[l].i = ei+i; row[l].j = ej+j; row[l].k = ek+k; row[l].c = 3;
						}
					}
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KA_local,ADD_VALUES);CHKERRQ(ierr);
				/*Assembling the righthand side vector f*/
				for(c = 0; c < veldof; c++){
					ierr = FLow_Vecf(RHS_local, &ctx->e3D, ek, ej, ei, c, ctx->flowprop);CHKERRQ(ierr);
					for(l = 0, k = 0; k < ctx->e3D.nphiz; k++){
						for(j = 0; j < ctx->e3D.nphiy; j++){
							for(i = 0; i < ctx->e3D.nphix; i++, l++){
								RHS_array[ek+k][ej+j][ei+i][c] += RHS_local[l];
							}
						}
					}
				}
				/*Assembling the righthand side vector g*/
				ierr = FLow_Vecg(RHS_local, &ctx->e3D,  ek, ej, ei, ctx->flowprop);CHKERRQ(ierr);
				//	for(i = 0; i < 8; i++)
				//		printf("a[%d] = %f\n", i, RHS_local[i]);
				for(l = 0, k = 0; k < ctx->e3D.nphiz; k++){
					for(j = 0; j < ctx->e3D.nphiy; j++){
						for(i = 0; i < ctx->e3D.nphix; i++, l++){
							RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
						}
					}
				}
			}
		}
	}
	
	if (xs+xm == nx) xm++;
	if (ys+ym == ny) ym++;
	if (zs+zm == nz) zm++;
	
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	PetscViewer viewer;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrixb.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = MatView(K,viewer);CHKERRQ(ierr);
	
	ierr = MatApplyFlowBC(K,ctx->daFlow,&ctx->bcFlow[0]);CHKERRQ(ierr);
		
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		
	ierr = DAVecRestoreArrayDOF(ctx->daFlow,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DALocalToGlobalBegin(ctx->daFlow,RHS_localVec,RHS);CHKERRQ(ierr);
	ierr = DALocalToGlobalEnd(ctx->daFlow,RHS_localVec,RHS);CHKERRQ(ierr);
	
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Vectorb.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = VecView(RHS,viewer);CHKERRQ(ierr);
	
	ierr = VecApplyFlowBC(RHS,&ctx->bcFlow[0], ctx->resprop);CHKERRQ(ierr);printf("Flowsolver worked\n");
	ierr = VecApplyWellFlowBC(RHS, ctx);
	ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = PetscFree2(row,row1);CHKERRQ(ierr);
	ierr = PetscFree3(KA_local,KB_local,RHS_local);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FLow_Vecg"
extern PetscErrorCode FLow_Vecg(PetscReal *Kg_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, FlowProp flowpropty)
{
	PetscInt		i, j, k, l;
	PetscInt		eg;
	PetscReal		beta_c, mu, rho, gamma_c, gx, gy, gz;
	
	PetscFunctionBegin;
	beta_c = flowpropty.beta;
	gamma_c = flowpropty.gamma;
	rho = flowpropty.rho;
	mu = flowpropty.mu;
	gx = flowpropty.g[0];
	gy = flowpropty.g[1];
	gz = flowpropty.g[2];
	for(l = 0, k = 0; k < e->nphiz; k++){
		for(j = 0; j < e->nphiy; j++){
			for(i = 0; i < e->nphix; i++, l++){
				Kg_local[l] = 0.;
				for(eg = 0; eg < e->ng; eg++){
					/* Need to multiply by the permability when it is available*/
					Kg_local[l] += -0.5 * rho * gamma_c *beta_c /mu * (gx * e->dphi[k][j][i][0][eg] + gy * e->dphi[k][j][i][1][eg] + gz * e->dphi[k][j][i][2][eg]) * e->weight[eg];
				}
			}
		}
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FLow_Vecf"
extern PetscErrorCode FLow_Vecf(PetscReal *Kf_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, FlowProp flowpropty)
{
	PetscInt		i, j, k, l;
	PetscInt		eg;
	PetscReal		beta_c, mu, rho, gamma_c, g;
	
	PetscFunctionBegin;
	beta_c = flowpropty.beta;
	gamma_c = flowpropty.gamma;
	rho = flowpropty.rho;
	mu = flowpropty.mu;
	g = flowpropty.g[c];

	for(l = 0, k = 0; k < e->nphiz; k++){
		for(j = 0; j < e->nphiy; j++){
			for(i = 0; i < e->nphix; i++, l++){
				Kf_ele[l] = 0.;
				for(eg = 0; eg < e->ng; eg++){
						/* Need to multiply by the permability when it is available*/
					Kf_ele[l] += 0.5 * g * gamma_c * beta_c * rho / mu * e->phi[k][j][i][eg] * e->weight[eg];
				}
			}
		}
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FLow_MatD"
extern PetscErrorCode FLow_MatD(PetscReal *Kd_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, FlowProp flowpropty)
{
	PetscInt		i, j, k, l;
	PetscInt		ii, jj, kk;
	PetscInt		eg;
	PetscReal		beta_c, mu;
	
	PetscFunctionBegin;
	beta_c = flowpropty.beta;
	mu = flowpropty.mu;
	for(l = 0, k = 0; k < e->nphiz; k++){
		for(j = 0; j < e->nphiy; j++){
			for(i = 0; i < e->nphix; i++){
				for(kk = 0; kk < e->nphiz; kk++){
					for(jj = 0; jj < e->nphiy; jj++){
						for(ii = 0; ii < e->nphix; ii++, l++){
							Kd_ele[l] = 0.;
							for(eg = 0; eg < e->ng; eg++){
								/* Need to multiply by the permability when it is available*/
								Kd_ele[l] += -0.5 * beta_c/mu * (e->dphi[k][j][i][0][eg] * e->dphi[kk][jj][ii][0][eg]
								+ e->dphi[k][j][i][1][eg] * e->dphi[kk][jj][ii][1][eg]
								+ e->dphi[k][j][i][2][eg] * e->dphi[kk][jj][ii][2][eg]) * e->weight[eg];
							}
						}
					}
				}
			}
		}
	}	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FLow_MatB"
extern PetscErrorCode FLow_MatB(PetscReal *KB_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c)
{
	PetscInt		i, j, k, l;
	PetscInt		ii, jj, kk;
	PetscInt		eg;
	PetscFunctionBegin;
	for(l = 0, k = 0; k < e->nphiz; k++){
		for(j = 0; j < e->nphiy; j++){
			for(i = 0; i < e->nphix; i++){
				for(kk = 0; kk < e->nphiz; kk++){
					for(jj = 0; jj < e->nphiy; jj++){
						for(ii = 0; ii < e->nphix; ii++, l++){
							KB_ele[l] = 0.;
							for(eg = 0; eg < e->ng; eg++){
								KB_ele[l] += -(e->phi[k][j][i][eg] * e->dphi[kk][jj][ii][c][eg] 
											   + 0.5 * e->dphi[k][j][i][c][eg] * e->phi[kk][jj][ii][eg]) * e->weight[eg];
							}
						}
					}
				}
			}
		}
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FLow_MatA"
extern PetscErrorCode FLow_MatA(PetscReal *A_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei)
{
	PetscInt		i, j, k, l;
	PetscInt		ii, jj, kk;
	PetscInt		eg;
	
	PetscFunctionBegin;
	for(l = 0, k = 0; k < e->nphiz; k++){
		for(j = 0; j < e->nphiy; j++){
			for(i = 0; i < e->nphix; i++){
				for(kk = 0; kk < e->nphiz; kk++){
					for(jj = 0; jj < e->nphiy; jj++){
						for(ii = 0; ii < e->nphix; ii++, l++){
							A_local[l] = 0;
							for(eg = 0; eg < e->ng; eg++){
								A_local[l] += 0.5 * e->phi[k][j][i][eg] * e->phi[kk][jj][ii][eg] * e->weight[eg];
							}
						}
					}
				}
			}
		}
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowSolverInitialize"
extern PetscErrorCode FlowSolverInitialize(VFCtx *ctx)
{
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"","");CHKERRQ(ierr);
	{
		ctx->flowrate = 50.;														/*Flow rate, assumed to be consistent with whatever flow units used for now*/
		ctx->units = FieldUnits;
		ierr = PetscOptionsEnum("-flowunits","\n\tFlow solver","",FlowUnitName,(PetscEnum)ctx->units,(PetscEnum*)&ctx->units,PETSC_NULL);CHKERRQ(ierr);
		ctx->flowcase = ALLPRESSUREBC;
		ierr = PetscOptionsEnum("-flow boundary conditions","\n\tFlow solver","",FlowBC_Case,(PetscEnum)ctx->flowcase,(PetscEnum*)&ctx->flowcase,PETSC_NULL);CHKERRQ(ierr);
	}
	ierr = PetscOptionsEnd();CHKERRQ(ierr);

	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
	if (comm_size == 1) {
		ierr = DAGetMatrix(ctx->daFlow,MATSEQAIJ,&ctx->KVelP);CHKERRQ(ierr);
	} else {
		ierr = DAGetMatrix(ctx->daFlow,MATMPIAIJ,&ctx->KVelP);CHKERRQ(ierr);
	}
	ierr = MatSetOption(ctx->KVelP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = DACreateGlobalVector(ctx->daFlow,&ctx->RHSVelP);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) ctx->RHSVelP,"RHS of flow solver");CHKERRQ(ierr);
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&ctx->kspVelP);CHKERRQ(ierr);
	
	ierr = KSPSetTolerances(ctx->kspVelP,1.e-8,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSetOperators(ctx->kspVelP,ctx->KVelP,ctx->KVelP,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ctx->kspVelP,PETSC_TRUE);CHKERRQ(ierr);
	ierr = KSPAppendOptionsPrefix(ctx->kspVelP,"VelP_");CHKERRQ(ierr);
	ierr = KSPSetType(ctx->kspVelP,KSPBCGS);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ctx->kspVelP);CHKERRQ(ierr);
	ierr = KSPGetPC(ctx->kspVelP,&ctx->pcVelP);CHKERRQ(ierr);
	ierr = PCSetType(ctx->pcVelP,PCSOR);CHKERRQ(ierr);
	ierr = PCSetFromOptions(ctx->pcVelP);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowSolverFinalize"
extern PetscErrorCode FlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	
	PetscErrorCode ierr;
	PetscFunctionBegin;
    ierr = DADestroy(ctx->daFlow);CHKERRQ(ierr);
	ierr = VecDestroy(fields->VelnPress);CHKERRQ(ierr);
	ierr = KSPDestroy(ctx->kspVelP);CHKERRQ(ierr);
	ierr = MatDestroy(ctx->KVelP);CHKERRQ(ierr);
	ierr = VecDestroy(ctx->RHSVelP);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetFlowProp"
extern PetscErrorCode GetFlowProp(FlowProp *flowprop, 	FlowUnit flowunit, ResProp resprop)
{
	PetscFunctionBegin;
	switch (flowunit) {
		case FieldUnits:
			flowprop->mu = resprop.visc;							/*viscosity in cp*/
			flowprop->rho = 62.428*resprop.fdens;					/*density in lb/ft^3*/
			flowprop->cf = resprop.cf;								/*compressibility in psi^{-1}*/
			flowprop->beta = 1.127;									/*flow rate conversion constant*/
			flowprop->gamma = 2.158e-4;								/*pressue conversion constant*/
			flowprop->alpha = 5.615;								/*volume conversion constatnt*/
			flowprop->g[0] = 0.;									/*x-componenet of gravity. unit is ft/s^2*/
			flowprop->g[1] = 0.;									/*y-component of gravity. unit is ft/s^2*/
			flowprop->g[2] = 32.17;									/*z-component of gravity. unit is ft/s^2*/
			break;
		case MetricUnits:
			flowprop->mu = 0.001*resprop.visc;						/*viscosity in Pa.s*/
			flowprop->rho = 1000*resprop.fdens;						/*density in kg/m^3*/
			flowprop->cf = 1.450e-4*resprop.cf;						/*compressibility in Pa^{-1}*/
			flowprop->beta = 86.4e-6;								/*flow rate conversion constant*/
			flowprop->gamma = 1e-3;									/*pressue conversion constant*/
			flowprop->alpha = 1;									/*volume conversion constatnt*/
			flowprop->g[0] = 0.;									/*x-componenet of gravity. unit is m/s^2*/
			flowprop->g[1] = 0.;									/*y-component of gravity. unit is m/s^2*/
			flowprop->g[2] = 9.81;									/*z-component of gravity. unit is m/s^2*/
			break;			
		default:
			SETERRQ2(PETSC_ERR_USER,"ERROR: [%s] unknown FLOWCASE %i .\n",__FUNCT__,flowunit);
			break;
	}
	PetscFunctionReturn(0);
}



