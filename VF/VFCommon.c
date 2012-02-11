#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFCracks.h"
#include "VFWell.h"
#include "../Utils/xdmf.h"

#undef __FUNCT__
#define __FUNCT__ "VFLogInitialize"
/*
  VFLogInitialize

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFLogInitialize(VFLog *vflog)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogStageRegister("I/O operations",&vflog->VF_IOStage);CHKERRQ(ierr);

  ierr = PetscLogStageRegister("U assembly",&vflog->VF_UAssemblyStage);CHKERRQ(ierr);
  ierr = PetscCookieRegister("U Mat local",&vflog->VF_MatULocalCookie);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("U Mat local",vflog->VF_MatULocalCookie,&vflog->VF_MatULocalEvent);CHKERRQ(ierr);
  ierr = PetscCookieRegister("U VecView local",&vflog->VF_VecULocalCookie);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("U Vec local",vflog->VF_VecULocalCookie,&vflog->VF_VecULocalEvent);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("U solver",&vflog->VF_USolverStage);CHKERRQ(ierr);

  ierr = PetscLogStageRegister("V assembly",&vflog->VF_VAssemblyStage);CHKERRQ(ierr);
  ierr = PetscCookieRegister("V Mat local",&vflog->VF_MatVLocalCookie);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("V Mat local",vflog->VF_MatVLocalCookie,&vflog->VF_MatVLocalEvent);CHKERRQ(ierr);
  ierr = PetscCookieRegister("V Vec local",&vflog->VF_VecVLocalCookie);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("V Vec local",vflog->VF_VecVLocalCookie,&vflog->VF_VecVLocalEvent);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("V solver",&vflog->VF_VSolverStage);CHKERRQ(ierr);
  
  ierr = PetscLogStageRegister("Energy",&vflog->VF_EnergyStage);CHKERRQ(ierr);
  ierr = PetscCookieRegister("Energy",&vflog->VF_EnergyLocalCookie);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Energy",vflog->VF_EnergyLocalCookie,&vflog->VF_EnergyLocalEvent);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFCtxGet"
/*
  VFCtxGet

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFCtxGet(VFCtx *ctx)
{
  PetscErrorCode      ierr;
  PetscInt            nopt;
  PetscTruth          flg;
  int                 i; 
  PetscTruth          hashelp;
  PetscReal          *buffer;
  char                prefix[PETSC_MAX_PATH_LEN+1];

  PetscFunctionBegin;
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-help",&hashelp,&flg);CHKERRQ(ierr);
  
  /* 
    Get options and register help messge
  */
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\n\nVF: general options:","");CHKERRQ(ierr);
  {
    ctx->verbose = 0;
    ierr = PetscOptionsInt("-verbose","\n\tDisplay debug informations about the computation\t","",ctx->verbose,&ctx->verbose,PETSC_NULL);CHKERRQ(ierr);
    ctx->mode = FRACTURE;
    ierr = PetscOptionsEnum("-mode","\n\tType of simulation","",VFModeName,(PetscEnum)ctx->mode,(PetscEnum*)&ctx->mode,PETSC_NULL);CHKERRQ(ierr);
    nopt = 6;
    ierr = PetscMalloc(nopt * sizeof(PetscReal),&buffer); CHKERRQ(ierr);
    for (i = 0;i < nopt;i++) {
      buffer[i]=0.;
    }
    ctx->hasInsitu = PETSC_FALSE;
    ierr = PetscOptionsRealArray("-insitumin","\n\tIn-situ stresses in the lower z-section.\n\tUse Voigt notations (s11, s22, s33, s23, s13, s12)","",buffer,&nopt,&flg);CHKERRQ(ierr);    
    if (nopt > 6 && !hashelp) {
      SETERRQ2(PETSC_ERR_USER,"ERROR: Expecting at most 6 component of the insitu stresses, got %i in %s\n",nopt,__FUNCT__);
    }
    for (i = 0;i < 6;i++) {
      ctx->insitumin[i]=buffer[i];
    }
    ctx->hasInsitu = flg;
    
    nopt = 6;
    ierr = PetscOptionsRealArray("-insitumax","\n\tIn-situ stresses in the upper z-section, will re-use insitumin if omitted","",buffer,&nopt,&flg);CHKERRQ(ierr);    
    if (nopt == 0) {
      for (i = 0;i < 6;i++) {
        ctx->insitumax[i]=ctx->insitumin[i];
      }
    } else {
      if (nopt > 6 && !hashelp) {
        SETERRQ2(PETSC_ERR_USER,"ERROR: Expecting at most 6 component of the insitu stresses, got %i in %s\n",nopt,__FUNCT__);
      }
      for (i = 0;i < 6;i++) {
        ctx->insitumax[i]=buffer[i];
      }
    }
    if (ctx->hasInsitu || flg) ctx->hasInsitu = PETSC_TRUE;
    
    ctx->nlayer = 1;
    ierr = PetscOptionsInt("-nlayer","\n\tNumber of layers","",ctx->nlayer,&ctx->nlayer,PETSC_NULL);CHKERRQ(ierr);
    nopt = ctx->nlayer-1;
    ierr = PetscMalloc((ctx->nlayer) * sizeof(PetscReal),&ctx->layersep);CHKERRQ(ierr);
    ctx->layersep[0] = -1e+30;
    ctx->layersep[1] = 0.;    
    ierr = PetscOptionsRealArray("-layersep","\n\tComma separated list of (nlayer-1) layer interfaces","",&ctx->layersep[1],&nopt,PETSC_NULL);CHKERRQ(ierr);    
    if (ctx->nlayer !=1 && nopt != ctx->nlayer-1 && !hashelp) {
      SETERRQ3(PETSC_ERR_USER,"ERROR: Expecting %i layers separators, got only %i in %s\n",ctx->nlayer+1,nopt,__FUNCT__);
    }
    ierr = PetscSNPrintf(ctx->prefix,FILENAME_MAX,"TEST");CHKERRQ(ierr);
    ierr = PetscOptionsString("-p","\n\tfile prefix","",ctx->prefix,ctx->prefix,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);

    ctx->altmintol  = 1.e-4;
    ierr = PetscOptionsReal("-altmintol","\n\tTolerance for alternate minimizations algorithm","",ctx->altmintol,&ctx->altmintol,PETSC_NULL);CHKERRQ(ierr);
    ctx->altminmaxit= 10000;
    ierr = PetscOptionsInt("-altminmaxit","\n\tMaximum number of altername minimizations iterations","",ctx->altminmaxit,&ctx->altminmaxit,PETSC_NULL);CHKERRQ(ierr);
    ctx->preset = SYMXY;
    ierr = PetscOptionsEnum("-preset","\n\tPreset simulation type","",VFPresetName,(PetscEnum)ctx->preset,(PetscEnum*)&ctx->preset,PETSC_NULL);CHKERRQ(ierr);
    ctx->unilateral = UNILATERAL_NONE;
    ierr = PetscOptionsEnum("-unilateral","\n\tType of unilateral conditions","",VFUnilateralName,(PetscEnum)ctx->unilateral,(PetscEnum*)&ctx->unilateral,PETSC_NULL);CHKERRQ(ierr);
    ctx->coupling = COUPLING_FULL;
    ierr = PetscOptionsEnum("-coupling","\n\tCoupling type","",VFCouplingName,(PetscEnum)ctx->coupling,(PetscEnum*)&ctx->coupling,PETSC_NULL);CHKERRQ(ierr);
    ctx->fileformat = FILEFORMAT_HDF5;
    ierr = PetscOptionsEnum("-format","\n\tFileFormat","",VFFileFormatName,(PetscEnum)ctx->fileformat,(PetscEnum*)&ctx->fileformat,PETSC_NULL);CHKERRQ(ierr);
    ctx->hasCrackPressure = PETSC_TRUE;
    ierr = PetscOptionsTruth("-pressurize","\n\tPressurize cracks","",ctx->hasCrackPressure,&ctx->hasCrackPressure,PETSC_NULL);CHKERRQ(ierr);
    
    ctx->numCracks = 0;
    ierr = PetscOptionsInt("-nc","\n\tNumber of penny-shaped cracks to insert","",ctx->numCracks,&ctx->numCracks,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscMalloc(ctx->numCracks*sizeof(VFPennyCrack),&ctx->crack);CHKERRQ(ierr);
    for (i = 0; i < ctx->numCracks; i++) {
      ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"c%d_",i);CHKERRQ(ierr);
      ierr = VFPennyCrackCreate(&(ctx->crack[i]));CHKERRQ(ierr);
      ierr = VFPennyCrackGet(prefix, &(ctx->crack[i]));CHKERRQ(ierr);
      ierr = VFPennyCrackView(&(ctx->crack[i]),PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }

    ctx->numWells = 0;
    ierr = PetscOptionsInt("-nw","\n\tNumber of wells to insert","",ctx->numWells,&ctx->numWells,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscMalloc(ctx->numWells*sizeof(VFWell),&ctx->well);CHKERRQ(ierr);
    for (i = 0; i < ctx->numWells; i++) {
      ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"w%d_",i);CHKERRQ(ierr);
      ierr = VFWellCreate(&(ctx->well[i]));CHKERRQ(ierr);
      ierr = VFWellGet(prefix, &(ctx->well[i]));CHKERRQ(ierr);
      ierr = VFWellView(&(ctx->well[i]),PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  ctx->timestep = 1;
  ctx->timevalue = 0.;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFPropGet"
/*
  VFPropGet

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFPropGet(VFProp *vfprop)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\n\nVF: variational fracture model specific options:","");CHKERRQ(ierr);
  {
    vfprop->epsilon  = 2.e-1;
    vfprop->eta      = 1.e-5;
    vfprop->irrevtol = 5e-2;
    ierr = PetscOptionsReal("-epsilon","\n\tVariational fracture regularization length (start from 4x cell size)","",vfprop->epsilon,&vfprop->epsilon,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-eta","\n\tArtificial stiffness of cracks ","",vfprop->eta,&vfprop->eta,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-irrevtol","\n\tThreshold on v below which fracture irreversibility is enforced","",vfprop->irrevtol,&vfprop->irrevtol,PETSC_NULL);CHKERRQ(ierr);
    vfprop->permmax = 5.;
    ierr = PetscOptionsReal("-permmax","\n\tPermeability multiplier of cracks (achieved at  v=0.)","",vfprop->permmax,&vfprop->permmax,PETSC_NULL);CHKERRQ(ierr);
    vfprop->atCv = .5;
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFLayerInit"
/*
  VFLayerInit: find the horizontal layer associated to a cell in the most stupid and unoptimzed way

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFLayerInit(VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       ek,l;
  PetscReal      ****coords_array;
  PetscInt       xs,xm,ys,ym,zs,zm;
  
  PetscFunctionBegin;
  ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = PetscMalloc((ctx->ncellz+1) * sizeof(PetscInt),&ctx->layer);CHKERRQ(ierr);

  for (l = 0; l < ctx->ncellz+1; l++) {
    ctx->layer[l] = 0;
  }
  for (ek = zs; ek < zs+zm; ek++) {
    for (l = 0; l < ctx->nlayer; l++){
      if (coords_array[ek][ys][xs][2] > ctx->layersep[l]) {
        ctx->layer[ek] = l;
      }
    }
  }
  if (ctx->verbose > 0) {
    for (ek = 0; ek < ctx->ncellz+1; ek++) {
      if (ek >= zs && ek < zs+zm) {
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ctx->layer[%i]=%i (depth=%g)\n",ek,ctx->layer[ek],coords_array[ek][ys][xs][2]);CHKERRQ(ierr);
      } else {
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"layer %i inactive\n",ek);CHKERRQ(ierr);        
      }
      ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    }
  }
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatPropGet"
/*
  MatPropGet: get all material properties.
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFMatPropGet(MatProp *matprop,PetscInt n)
{
  PetscErrorCode ierr;
  PetscReal      *prop;
  PetscInt       nopt;
  PetscReal      E = 1.;
  PetscReal      nu = 0.;
  PetscReal      alpha = 1.e-5;
  PetscReal      Gc = 1.;
  PetscReal      beta = 1.e-6;
  int            i;
  
  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\n\nVF: material properties:","");CHKERRQ(ierr);
  {  
    ierr = PetscMalloc(n * sizeof(PetscReal),&prop);CHKERRQ(ierr);
    nopt = n;
    for (i = 0; i < n; i++) prop[i] = E;
    ierr = PetscOptionsRealArray("-E","\n\tComma separated list of Young\'s modulii","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) {
      SETERRQ4(PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-E",nopt,__FUNCT__);
    }
    for (i=0; i< n; i++) matprop[i].E = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = nu;
    ierr = PetscOptionsRealArray("-nu","\n\tComma separated list of Poisson\'s ratio","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) {
      SETERRQ4(PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-nu",nopt,__FUNCT__);
    }
    for (i=0; i< n; i++) matprop[i].nu = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = alpha;
    ierr = PetscOptionsRealArray("-alpha","\n\tComma separated list of linear thermal expansion coefficients","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) {
      SETERRQ4(PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-alpha",nopt,__FUNCT__);
    }
    for (i=0; i< n; i++) matprop[i].alpha = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = Gc;
    ierr = PetscOptionsRealArray("-Gc","\n\tComma separated list of fracture toughness (not critical SIF!)","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) {
      SETERRQ4(PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-Gc",nopt,__FUNCT__);
    }
    for (i=0; i< n; i++) matprop[i].Gc = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = beta;
    ierr = PetscOptionsRealArray("-beta","\n\tComma separated list of Biot constants","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) {
      SETERRQ4(PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-beta",nopt,__FUNCT__);
    }
    for (i=0; i< n; i++) matprop[i].beta = prop[i];
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  ierr = PetscFree(prop);CHKERRQ(ierr);

  for (i=0;i < n; i++) {
    matprop[i].lambda = matprop[i].E * matprop[i].nu / (1. + matprop[i].nu) / (1. - 2. * matprop[i].nu);
    matprop[i].mu     = matprop[i].E / (1. + matprop[i].nu) * .5;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFGeometryInitialize"
/*
  VFGeometryInitialize: Creates DA, and coordinates, and other geometric informations

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFGeometryInitialize(VFCtx *ctx,PetscReal *dx,PetscReal *dy,PetscReal *dz)
{
  PetscErrorCode      ierr;
  PetscViewer         viewer,h5viewer;
  char                filename[FILENAME_MAX];  
  PetscReal           *X,*Y,*Z;
  PetscReal           ****coords_array;
  PetscInt            xs,xm,ys,ym,zs,zm;
  int                 i,j,k;
  PetscReal           BBmin[3],BBmax[3];

  PetscFunctionBegin;
  ierr = CartFE_Init();CHKERRQ(ierr);

  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,ctx->ncellx+1,ctx->ncelly+1,ctx->ncellz+1,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&ctx->daVect);CHKERRQ(ierr);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,ctx->ncellx+1,ctx->ncelly+1,ctx->ncellz+1,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&ctx->daScal);CHKERRQ(ierr);
  ierr = CartFE_Element3DCreate(&ctx->e3D);CHKERRQ(ierr);
  /*
    Constructs coordinates Vec
  */
  ierr = DACreateGlobalVector(ctx->daVect,&ctx->coordinates);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->coordinates,"Coordinates");CHKERRQ(ierr);
  
  /*
    Construct coordinates from the arrays of cell sizes
  */
  ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = PetscMalloc3(ctx->ncellx+1,PetscReal,&X,ctx->ncelly+1,PetscReal,&Y,ctx->ncellz+1,PetscReal,&Z);CHKERRQ(ierr);
  Z[0] = 0.;
  Y[0] = 0.;
  X[0] = 0.;
  for (k = 1; k < ctx->ncellz+1; k++) Z[k] = Z[k-1] + dz[k-1];
  for (j = 1; j < ctx->ncelly+1; j++) Y[j] = Y[j-1] + dy[j-1];
  for (i = 1; i < ctx->ncellx+1; i++) X[i] = X[i-1] + dx[i-1];

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        coords_array[k][j][i][2] = Z[k];
        coords_array[k][j][i][1] = Y[j];
        coords_array[k][j][i][0] = X[i]; 
      }
    }
  }
  ierr = PetscFree3(X,Y,Z);CHKERRQ(ierr);
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);    

  ierr = DASetCoordinates(ctx->daVect,ctx->coordinates);CHKERRQ(ierr);
  ierr = DASetCoordinates(ctx->daScal,ctx->coordinates);CHKERRQ(ierr);
  switch (ctx->fileformat) {
    case FILEFORMAT_BIN:
      /*
        As of version 3.1, there is a bug in petsc preventing to save 2 DA with different number of degrees of freedoms
        per node in a single file.
        Even coordinates are part of the DA and do not need to be saved separately, it looks like the coordinate vector
        obtained from DAGetCoordinates cannot be saved properly in an hdf5 file, so we save the coordinate vector anyway
      */
      ierr = PetscLogStagePush(ctx->vflog.VF_IOStage);CHKERRQ(ierr);
      ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.bin",ctx->prefix);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = DAView(ctx->daScal,viewer);CHKERRQ(ierr);
      ierr = VecView(ctx->coordinates,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      ierr = PetscLogStagePop();CHKERRQ(ierr);
      break;
#ifdef PETSC_HAVE_HDF5
    case FILEFORMAT_HDF5:
      /*
        Write headers in multistep XDMF file
      */
      ierr = PetscLogStagePush(ctx->vflog.VF_IOStage);CHKERRQ(ierr);
      ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.xmf",ctx->prefix);CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&ctx->XDMFviewer);
      ierr = XDMFmultistepInitialize(ctx->XDMFviewer);CHKERRQ(ierr);
      
      /* 
        Save cordinate in main hdf5 file
      */
      ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.h5",ctx->prefix);CHKERRQ(ierr);
      ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&h5viewer);CHKERRQ(ierr);
      ierr = VecView(ctx->coordinates,h5viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(h5viewer);CHKERRQ(ierr);
      ierr = PetscLogStagePop();CHKERRQ(ierr);
#endif
  }
  
  /*
    Replace the coordinate vector defined on PETSC_COMM_WORLD with ghosted coordinate vectors on PETSC_COMM_SELF, 
    so that we can query coordinates of ghost points
  */
  ierr = VecDestroy(ctx->coordinates);CHKERRQ(ierr);
  ierr = DAGetGhostedCoordinates(ctx->daScal,&ctx->coordinates);CHKERRQ(ierr);
  
  if (ctx->verbose > 0) {
    ierr = DAGetBoundingBox(ctx->daVect,BBmin,BBmax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reservoir bounding box: (%g,%g) x (%g,%g) x (%g,%g)\n",BBmin[0],BBmax[0],BBmin[1],BBmax[1],BBmin[2],BBmax[2]);CHKERRQ(ierr);
    ierr = DAView(ctx->daVect,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = DAView(ctx->daScal,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = VFLayerInit(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFFieldsInitialize"
/*
  VFFieldsInitialize: Creates fields

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFFieldsInitialize(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode  ierr;
  PetscInt        c; 

  PetscFunctionBegin;
  ierr = DACreateGlobalVector(ctx->daVect,&fields->U);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->U,"Displacement");CHKERRQ(ierr);
  ierr = VecSet(fields->U,0.0);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(ctx->daVect,&fields->BCU);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->BCU,"Displacement");CHKERRQ(ierr);
  ierr = VecSet(fields->BCU,0.0);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(ctx->daScal,&fields->V);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->V,"Fracture");CHKERRQ(ierr);

  ierr = DACreateGlobalVector(ctx->daScal,&fields->VIrrev);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->VIrrev,"VIrrev");CHKERRQ(ierr);
  
  /*
    Create optional penny shaped cracks
  */
  ierr = VecSet(fields->V,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields->VIrrev,1.0);CHKERRQ(ierr);
  for (c = 0; c < ctx->numCracks; c++) {
    ierr = VFPennyCrackBuildVAT2(fields->V,&(ctx->crack[c]),ctx);CHKERRQ(ierr);
    ierr = VecPointwiseMin(fields->VIrrev,fields->V,fields->VIrrev);CHKERRQ(ierr);
  }
  /*
    Create optional wells
  */
  for (c = 0; c < ctx->numWells; c++) {
    ierr = VFWellBuildVAT2(fields->V,&(ctx->well[c]),ctx);CHKERRQ(ierr);
    ierr = VecPointwiseMin(fields->VIrrev,fields->V,fields->VIrrev);CHKERRQ(ierr);
  }
  ierr = VecCopy(fields->VIrrev,fields->V);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(ctx->daScal,&fields->theta);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->theta,"Temperature");CHKERRQ(ierr);
  ierr = VecSet(fields->theta,0.0);CHKERRQ(ierr);
  
  ierr = DACreateGlobalVector(ctx->daScal,&fields->thetaRef);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->thetaRef,"Reference Temperature");CHKERRQ(ierr);
  
  ierr = DACreateGlobalVector(ctx->daScal,&fields->pressure);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->pressure,"Pressure");CHKERRQ(ierr);
  ierr = VecSet(fields->pressure,0.0);CHKERRQ(ierr);
  
  ierr = DACreateGlobalVector(ctx->daScal,&fields->pressureRef);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->pressureRef,"Reference Pressure");CHKERRQ(ierr);
  ierr = VecSet(fields->pressureRef,0.0);CHKERRQ(ierr);
  
  ierr = DACreateGlobalVector(ctx->daScal,&fields->pmult);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->pmult,"Permeability");CHKERRQ(ierr);
  ierr = VecSet(fields->pmult,0.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFBCInitialize"
/*
  VFBCInitialize: Creates and initialize bondary condition data structures

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFBCInitialize(VFCtx *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = BCUInit(&ctx->bcU[0],ctx->preset);CHKERRQ(ierr);
  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"BCU:\n");CHKERRQ(ierr);
    ierr = BCView(&ctx->bcU[0],PETSC_VIEWER_STDOUT_WORLD,3);CHKERRQ(ierr);
  }
  ierr = BCVInit(&ctx->bcV[0],ctx->preset);CHKERRQ(ierr);
  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"BCV:\n");CHKERRQ(ierr);
    ierr = BCView(&ctx->bcV[0],PETSC_VIEWER_STDOUT_WORLD,1);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
  }

#undef __FUNCT__
#define __FUNCT__ "VFSolversInitialize"
/*
  VFSolversInitialize: Creates matrices, RHS, solvers

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFSolversInitialize(VFCtx *ctx)
{
  PetscMPIInt    comm_size;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
  if (comm_size == 1) {
    ierr = DAGetMatrix(ctx->daVect,MATSEQAIJ,&ctx->KU);CHKERRQ(ierr);
  } else {
    ierr = DAGetMatrix(ctx->daVect,MATMPIAIJ,&ctx->KU);CHKERRQ(ierr);
  }
  ierr = MatSetOption(ctx->KU,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(ctx->daVect,&ctx->RHSU);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->RHSU,"RHSU");CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ctx->kspU);CHKERRQ(ierr);
  
  ierr = KSPSetTolerances(ctx->kspU,1.e-8,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetOperators(ctx->kspU,ctx->KU,ctx->KU,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ctx->kspU,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPAppendOptionsPrefix(ctx->kspU,"U_");CHKERRQ(ierr);
  ierr = KSPSetType(ctx->kspU,KSPBCGSL);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ctx->kspU);CHKERRQ(ierr);
  ierr = KSPGetPC(ctx->kspU,&ctx->pcU);CHKERRQ(ierr);
  ierr = PCSetType(ctx->pcU,PCBJACOBI);CHKERRQ(ierr);
  ierr = PCSetFromOptions(ctx->pcU);CHKERRQ(ierr);
  
  if (comm_size == 1) {
    ierr = DAGetMatrix(ctx->daScal,MATSEQAIJ,&ctx->KV);CHKERRQ(ierr);
  } else {
    ierr = DAGetMatrix(ctx->daScal,MATMPIAIJ,&ctx->KV);CHKERRQ(ierr);
  }
  ierr = MatSetOption(ctx->KV,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(ctx->daScal,&ctx->RHSV);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->RHSV,"RHSV");CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ctx->kspV);CHKERRQ(ierr);
  
  ierr = KSPSetTolerances(ctx->kspV,1.e-8,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetOperators(ctx->kspV,ctx->KV,ctx->KV,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ctx->kspV,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPAppendOptionsPrefix(ctx->kspV,"V_");CHKERRQ(ierr);
  ierr = KSPSetType(ctx->kspV,KSPBCGSL);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ctx->kspV);CHKERRQ(ierr);
  ierr = KSPGetPC(ctx->kspV,&ctx->pcV);CHKERRQ(ierr);
  ierr = PCSetType(ctx->pcV,PCBJACOBI);CHKERRQ(ierr);
  ierr = PCSetFromOptions(ctx->pcV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFTimeStepPrepare"
/*
  VFTimeStepPrepare: Prepare for a new time step:
    - Update VIrrev
    - Read boundary displacement from files if necessary
    - Set boundary values of U and V

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFTimeStepPrepare(VFCtx *ctx,VFFields *fields)
{
  char           filename[FILENAME_MAX];
  PetscViewer    viewer;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /*
    Initialize VIrrev with the result of past iteration
  */
  ierr = VecCopy(fields->V,fields->VIrrev);CHKERRQ(ierr);
  
  /* 
    Set boundary values for U and V
  */
  ierr = VecApplyDirichletBC(fields->U,fields->BCU,&ctx->bcU[0]);CHKERRQ(ierr);
  ierr = VecApplyDirichletBC(fields->V,fields->V,&ctx->bcV[0]);CHKERRQ(ierr);

  if (ctx->coupling == COUPLING_NONE) {
    /*
      If we do an uncoupled computation, read Temperature and pressure from file.
    */
    ierr = PetscLogStagePush(ctx->vflog.VF_IOStage);CHKERRQ(ierr);
    switch (ctx->fileformat) {
      case FILEFORMAT_HDF5:
#ifdef PETSC_HAVE_HDF5
        ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
        ierr = VecLoadIntoVector(viewer,fields->theta);CHKERRQ(ierr);
        ierr = VecLoadIntoVector(viewer,fields->pressure);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);    
#endif
        break;
      case FILEFORMAT_BIN:
        SETERRQ(PETSC_ERR_SUP,"Reading from binary files not implemented yet\n");
        break;
      }
    ierr = PetscLogStagePop();CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFElasticityTimeStep"
/*
  VFElasticityTimeStep

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFElasticityTimeStep(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /*
    Assembly, solve
  */
  ierr = VF_StepU(fields,ctx);CHKERRQ(ierr);
  ctx->ElasticEnergy=0;
  ctx->InsituWork=0;
  ierr = VF_UEnergy3D(&ctx->ElasticEnergy,&ctx->InsituWork,&ctx->PressureWork,fields,ctx);CHKERRQ(ierr);
  ctx->TotalEnergy = ctx->ElasticEnergy - ctx->InsituWork;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFFractureTimeStep"
/*
  VFFractureTimeStep

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFFractureTimeStep(VFCtx *ctx,VFFields *fields)
{
  PetscInt        altminit=1;
  Vec             Vold;
  PetscReal       errV=1e+10;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = VecDuplicate(fields->V,&Vold);CHKERRQ(ierr);
  do {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i, alt min step %i\n",ctx->timestep,altminit);CHKERRQ(ierr);
    /*
      Assembly, solve for U
    */
    ierr = VF_StepU(fields,ctx);CHKERRQ(ierr);
    
    ierr = VecCopy(fields->V,Vold);CHKERRQ(ierr);
    /*
      Assembly, solve for V
    */
    ierr = VF_StepV(fields,ctx);CHKERRQ(ierr);

    /* 
      Compute max V change
    */
    ierr = VecAXPY(Vold,-1.,fields->V);CHKERRQ(ierr);
    ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on V: %e\n",errV);CHKERRQ(ierr);
    altminit++;
  } while (errV > ctx->altmintol && altminit <= ctx->altminmaxit);
  
  ctx->ElasticEnergy=0;
  ctx->InsituWork=0;
  ierr = VF_UEnergy3D(&ctx->ElasticEnergy,&ctx->InsituWork,&ctx->PressureWork,fields,ctx);CHKERRQ(ierr);
  ctx->SurfaceEnergy=0;
  ierr = VF_VEnergy3D(&ctx->SurfaceEnergy,fields,ctx);CHKERRQ(ierr);
  ctx->TotalEnergy = ctx->ElasticEnergy + ctx->SurfaceEnergy - ctx->InsituWork;
  ierr = VecDestroy(Vold);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFFinalize"
/*
  VFFinalize

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFFinalize(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;
  char            filename[FILENAME_MAX];
  
  PetscFunctionBegin;
  ierr = PetscFree(ctx->matprop);CHKERRQ(ierr);
  ierr = PetscFree(ctx->layer);CHKERRQ(ierr);
  
  ierr = DADestroy(ctx->daVect);CHKERRQ(ierr);
  ierr = DADestroy(ctx->daScal);CHKERRQ(ierr);

  ierr = KSPDestroy(ctx->kspU);CHKERRQ(ierr);
  ierr = MatDestroy(ctx->KU);CHKERRQ(ierr);
  ierr = VecDestroy(ctx->RHSU);CHKERRQ(ierr); 
  
  ierr = KSPDestroy(ctx->kspV);CHKERRQ(ierr);
  ierr = MatDestroy(ctx->KV);CHKERRQ(ierr);
  ierr = VecDestroy(ctx->RHSV);CHKERRQ(ierr); 
  ierr = VecDestroy(ctx->coordinates);CHKERRQ(ierr);
  
  ierr = VecDestroy(fields->U);CHKERRQ(ierr);
  ierr = VecDestroy(fields->BCU);CHKERRQ(ierr);
  ierr = VecDestroy(fields->V);CHKERRQ(ierr);
  ierr = VecDestroy(fields->VIrrev);CHKERRQ(ierr);
  ierr = VecDestroy(fields->theta);CHKERRQ(ierr);
  ierr = VecDestroy(fields->thetaRef);CHKERRQ(ierr);
  ierr = VecDestroy(fields->pressure);CHKERRQ(ierr);
  ierr = VecDestroy(fields->pressureRef);CHKERRQ(ierr);
  ierr = VecDestroy(fields->pmult);CHKERRQ(ierr);
  
  ierr = PetscViewerDestroy(ctx->energyviewer);CHKERRQ(ierr);
  
  ierr = PetscFree(ctx->crack);CHKERRQ(ierr);
  ierr = PetscFree(ctx->well);CHKERRQ(ierr);
  
  /*
    Close the xdmf multi-step file
  */
#ifdef PETSC_HAVE_HDF5
  if (ctx->fileformat == FILEFORMAT_HDF5) {
    ierr = XDMFuniformgridFinalize(ctx->XDMFviewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(ctx->XDMFviewer);CHKERRQ(ierr);
  }
#endif
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx->prefix);CHKERRQ(ierr);
  ierr = PetscLogPrintSummary(PETSC_COMM_WORLD,filename);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FieldsH5Write"
/*
  FieldsH5Write: Export all fields in HDF5 format using PETSc HDF5 viewer. Also write an XDMF container

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode FieldsH5Write(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;
  char           H5filename[FILENAME_MAX],H5coordfilename[FILENAME_MAX],XDMFfilename[FILENAME_MAX];
  PetscViewer    H5Viewer,XDMFViewer;
  
  PetscFunctionBegin;
  /*
    Write the hdf5 files
  */
#ifdef PETSC_HAVE_HDF5
  ierr = PetscSNPrintf(H5filename,FILENAME_MAX,"%s.%.5i.h5",ctx->prefix,ctx->timestep);CHKERRQ(ierr);
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,H5filename,FILE_MODE_WRITE,&H5Viewer);
  ierr = VecView(fields->U,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->V,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->pmult,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->theta,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->pressure,H5Viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(H5Viewer);CHKERRQ(ierr);
  
  /*
    Write the XDMF Description
  */
  ierr = PetscSNPrintf(XDMFfilename,FILENAME_MAX,"%s.%.5i.xmf",ctx->prefix,ctx->timestep);CHKERRQ(ierr);
  ierr = PetscSNPrintf(H5coordfilename,FILENAME_MAX,"%s.h5",ctx->prefix);CHKERRQ(ierr);
  
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,XDMFfilename,&XDMFViewer);CHKERRQ(ierr);
  ierr = XDMFuniformgridInitialize(XDMFViewer,ctx->timevalue,H5filename);CHKERRQ(ierr); 
  ierr = XDMFtopologyAdd (XDMFViewer,ctx->ncellx+1,ctx->ncelly+1,ctx->ncellz+1,H5coordfilename,"Coordinates");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,ctx->ncellx+1,ctx->ncelly+1,ctx->ncellz+1,3,"Vector","Node",H5filename,"Displacement");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,ctx->ncellx+1,ctx->ncelly+1,ctx->ncellz+1,1,"Scalar","Node",H5filename,"Fracture");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,ctx->ncellx+1,ctx->ncelly+1,ctx->ncellz+1,1,"Scalar","Node",H5filename,"Permeability");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,ctx->ncellx+1,ctx->ncelly+1,ctx->ncellz+1,1,"Scalar","Node",H5filename,"Temperature");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,ctx->ncellx+1,ctx->ncelly+1,ctx->ncellz+1,1,"Scalar","Node",H5filename,"Pressure");CHKERRQ(ierr);
  ierr = XDMFuniformgridFinalize(XDMFViewer);CHKERRQ(ierr); 
  ierr = PetscViewerDestroy(XDMFViewer);CHKERRQ(ierr);
  /*
    Add an entry in the multistep XDMF file
  */
  ierr = XDMFmultistepAddstep(ctx->XDMFviewer,XDMFfilename);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FieldsBinaryWrite"
/*
  FieldsBinaryWrite: Export all fields in PETSc native binary format.

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode FieldsBinaryWrite(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;
  char           filename[FILENAME_MAX];
  PetscViewer    viewer;
  
  PetscFunctionBegin;
  /*
    Write the binary files
  */
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.%.5i.bin",ctx->prefix,ctx->timestep);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
  ierr = VecView(fields->U,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->V,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->pmult,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->theta,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->pressure,viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PermUpdate"
/*
  PermUpdate

  (c) 2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode PermUpdate(Vec V,Vec Pmult,VFProp *vfprop,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscReal      pmult,v;
  PetscInt       xs,xm;
  PetscInt       ys,ym;
  PetscInt       zs,zm;
  PetscInt       i,j,k;
  PetscReal      ***v_array,***pmult_array;

  PetscFunctionBegin;

  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);        
  ierr = DAVecGetArray(ctx->daScal,Pmult,&pmult_array);CHKERRQ(ierr);        
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        v = v_array[k][j][i];
        if (v <= 0) {
          pmult = vfprop->permmax;
        } else if (v <= 1.) {
				  pmult = vfprop->permmax - (vfprop->permmax - 1.0) * v;
        } else {
          pmult = 1.;
        }
        pmult_array[k][j][i] = pmult;
      }
    }
  }
  ierr = DAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);        
  ierr = DAVecRestoreArray(ctx->daScal,Pmult,&pmult_array);CHKERRQ(ierr);        
  PetscFunctionReturn(0);
}
