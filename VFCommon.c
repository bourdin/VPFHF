#include "petsc.h"
#include "CartFE.h"
#include "VFCommon_private.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFHeat.h"
#include "VFWell.h"
#include "VFCracks.h"
#include "VFHeat.h"

#include "xdmf.h"

#undef __FUNCT__
#define __FUNCT__ "VFInitialize"
/*
 VFInitialize: Initialize the VF code. Called by the fortran implementation of VIADAT

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFInitialize(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;
  PetscViewer    optionsviewer;
  char           filename[FILENAME_MAX];

  PetscFunctionBegin;
  ierr = PetscPrintf(PETSC_COMM_WORLD,banner);CHKERRQ(ierr);
#if defined(PETSC_USE_DEBUG)
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      ##########################################################\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                          WARNING!!!                    #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #   For production runs, use a petsc compiled with       #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #   optimization, the performance will be generally      #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #   two or three times faster.                           #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      ##########################################################\n\n\n");CHKERRQ(ierr);
#endif

  ierr = PetscLogBegin();CHKERRQ(ierr);
  ierr = VFCtxGet(ctx);CHKERRQ(ierr);
  ierr = VFGeometryInitialize(ctx);CHKERRQ(ierr);
  ierr = VFPropGet(&ctx->vfprop);CHKERRQ(ierr);

  ierr = PetscMalloc(ctx->nlayer * sizeof(VFMatProp),&ctx->matprop);CHKERRQ(ierr);
  ierr = VFMatPropGet(ctx->matprop,ctx->nlayer);CHKERRQ(ierr);
  ierr = VFResPropGet(&ctx->resprop);CHKERRQ(ierr);
  ierr = VFMatPropFieldsInitialize(ctx, ctx->matprop);

  if (ctx->printhelp) PetscFunctionReturn(0);


  ierr = VFFieldsInitialize(ctx,fields);CHKERRQ(ierr);
  ierr = VFBCInitialize(ctx);CHKERRQ(ierr);
  ierr = VFSolversInitialize(ctx);CHKERRQ(ierr);

  ctx->heatsolver = HEATSOLVER_SNESFEM;
  ierr            = FlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
  ierr            = VF_HeatSolverInitialize(ctx,fields);CHKERRQ(ierr);
  /*
   Save command line options to a file
   */
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.txt",ctx->prefix,0);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&optionsviewer);CHKERRQ(ierr);
  ierr = PetscOptionsView(optionsviewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&optionsviewer);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Option table:\n");CHKERRQ(ierr);
  ierr = PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.ener",ctx->prefix);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&ctx->energyviewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->energyviewer,"#i,Elastic Energy,InsituWork,Surface Energy,Pressure Work,Total Energy\n");CHKERRQ(ierr);
  ierr = PetscViewerFlush(ctx->energyviewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFCtxGet"
/*
 VFCtxGet

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFCtxGet(VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       nopt;
  PetscBool      flg;
  int            i;
  PetscReal      *buffer;
  char           prefix[PETSC_MAX_PATH_LEN+1];

  PetscFunctionBegin;
  /*
   Get options and register help message
   */
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\n\nVF: general options:","");CHKERRQ(ierr);
  {
    ctx->verbose = 0;
    ierr         = PetscOptionsInt("-verbose","\n\tDisplay debug informations about the computation\t","",ctx->verbose,&ctx->verbose,PETSC_NULL);CHKERRQ(ierr);
    /*
      ctx->flowsolver = FLOWSOLVER_KSPMIXEDFEM;
    */
    ierr          = PetscOptionsEnum("-flowsolver","\n\tFlow solver","",VFFlowSolverName,(PetscEnum)ctx->flowsolver,(PetscEnum*)&ctx->flowsolver,PETSC_NULL);CHKERRQ(ierr);
    ctx->units    = UnitaryUnits;
    ierr          = PetscOptionsEnum("-flowunits","\n\tFlow solver","",VFUnitName,(PetscEnum)ctx->units,(PetscEnum*)&ctx->units,PETSC_NULL);CHKERRQ(ierr);
    ctx->flowcase = ALLNORMALFLOWBC;
    ierr          = PetscOptionsEnum("-flowBC","\n\tFlow boundary conditions","",VFFlowBC_Case,(PetscEnum)ctx->flowcase,(PetscEnum*)&ctx->flowcase,PETSC_NULL);CHKERRQ(ierr);

    ctx->fractureflowsolver = FRACTUREFLOWSOLVER_SNESMIXEDFEM;
    ierr                    = PetscOptionsEnum("-fractureflowsolver","\n\tFracture Flow solver","",VFFractureFlowSolverName,(PetscEnum)ctx->fractureflowsolver,(PetscEnum*)&ctx->fractureflowsolver,PETSC_NULL);CHKERRQ(ierr);

    ctx->heatsolver = HEATSOLVER_SNESFEM;
    ierr            = PetscOptionsEnum("-heatsolver","\n\tHeat solver","",VFHeatSolverName,(PetscEnum)ctx->heatsolver,(PetscEnum*)&ctx->heatsolver,PETSC_NULL);CHKERRQ(ierr);
    ctx->Hunits     = UnitaryUnits;
    ierr            = PetscOptionsEnum("-heatunits","\n\tHeat solver units","",VFUnitName,(PetscEnum)ctx->Hunits,(PetscEnum*)&ctx->Hunits,PETSC_NULL);CHKERRQ(ierr);

    ctx->mechsolver = FRACTURE;
    ierr            = PetscOptionsEnum("-mechsolver","\n\tType of simulation","",VFMechSolverName,(PetscEnum)ctx->mechsolver,(PetscEnum*)&ctx->mechsolver,PETSC_NULL);CHKERRQ(ierr);
    ctx->hasInsitu  = PETSC_FALSE;
    nopt            = 6;
    ierr            = PetscMalloc(nopt * sizeof(PetscReal),&buffer);CHKERRQ(ierr);
    for (i = 0; i < nopt; i++) buffer[i]=0.;
    ierr = PetscOptionsRealArray("-insitumin","\n\tIn-situ stresses in the lower z-section.\n\tUse Voigt notations (s11, s22, s33, s23, s13, s12)","",buffer,&nopt,&flg);CHKERRQ(ierr);
    if (nopt > 6 && !ctx->printhelp) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting at most 6 component of the insitu stresses, got %i in %s\n",nopt,__FUNCT__);
    for (i = 0; i < 6; i++) ctx->insitumin[i]=buffer[i];
    ctx->hasInsitu = flg;

    nopt = 6;
    ierr = PetscOptionsRealArray("-insitumax","\n\tIn-situ stresses in the upper z-section, will re-use insitumin if omitted","",buffer,&nopt,&flg);CHKERRQ(ierr);
    if (nopt == 0)
      for (i = 0; i < 6; i++) {
        ctx->insitumax[i]=ctx->insitumin[i];
    } else {
      if (nopt > 6 && !ctx->printhelp) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting at most 6 component of the insitu stresses, got %i in %s\n",nopt,__FUNCT__);
      for (i = 0; i < 6; i++) ctx->insitumax[i]=buffer[i];
    }
    if (ctx->hasInsitu || flg) ctx->hasInsitu = PETSC_TRUE;

    /*
     Move to geometry
     */
    ctx->nlayer      = 1;
    ierr             = PetscOptionsInt("-nlayer","\n\tNumber of layers","",ctx->nlayer,&ctx->nlayer,PETSC_NULL);CHKERRQ(ierr);
    nopt             = ctx->nlayer-1;
    ierr             = PetscMalloc((ctx->nlayer+1) * sizeof(PetscReal),&ctx->layersep);CHKERRQ(ierr);
    ctx->layersep[0] = -1e+30;
    ctx->layersep[1] = 0.;
    ierr             = PetscOptionsRealArray("-layersep","\n\tComma separated list of (nlayer-1) layer interfaces","",&ctx->layersep[1],&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (ctx->nlayer !=1 && nopt != ctx->nlayer-1 && !ctx->printhelp) SETERRQ3(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting %i layers separators, got only %i in %s\n",ctx->nlayer+1,nopt,__FUNCT__);
    /*
     end move to geometry
     */
    ierr = PetscSNPrintf(ctx->prefix,FILENAME_MAX,"TEST");CHKERRQ(ierr);
    ierr = PetscOptionsString("-p","\n\tfile prefix","",ctx->prefix,ctx->prefix,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);

    ctx->altmintol  = 1.e-4;
    ierr            = PetscOptionsReal("-altmintol","\n\tTolerance for alternate minimizations algorithm","",ctx->altmintol,&ctx->altmintol,PETSC_NULL);CHKERRQ(ierr);
    ctx->altminmaxit= 10000;
    ierr            = PetscOptionsInt("-altminmaxit","\n\tMaximum number of alternate minimizations iterations","",ctx->altminmaxit,&ctx->altminmaxit,PETSC_NULL);CHKERRQ(ierr);
    ctx->preset     = SYMXY;
    ierr            = PetscOptionsEnum("-preset","\n\tPreset simulation type","",VFPresetName,(PetscEnum)ctx->preset,(PetscEnum*)&ctx->preset,PETSC_NULL);CHKERRQ(ierr);
    ctx->unilateral = UNILATERAL_NONE;
    ierr            = PetscOptionsEnum("-unilateral","\n\tType of unilateral conditions","",VFUnilateralName,(PetscEnum)ctx->unilateral,(PetscEnum*)&ctx->unilateral,PETSC_NULL);CHKERRQ(ierr);
    ctx->fileformat = FILEFORMAT_HDF5;
    ierr            = PetscOptionsEnum("-format","\n\tFileFormat","",VFFileFormatName,(PetscEnum)ctx->fileformat,(PetscEnum*)&ctx->fileformat,PETSC_NULL);CHKERRQ(ierr);

    ctx->maxtimestep  = 1;
    ierr              = PetscOptionsInt("-maxtimestep","\n\tMaximum number of timestep","",ctx->maxtimestep,&ctx->maxtimestep,PETSC_NULL);CHKERRQ(ierr);
    ctx->maxtimevalue = 1.;
    ierr              = PetscOptionsReal("-maxtimevalue","\n\tMaximum timevalue","",ctx->maxtimevalue,&ctx->maxtimevalue,PETSC_NULL);CHKERRQ(ierr);

    nopt = 6;
    ierr = PetscOptionsRealArray("-BCpres", "\n\tPressure at Boundaries.\n\t (PX0,PX1,PY0,PY1,PZ0,PZ1) negative value if natural BC","",buffer,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt > 6 && !ctx->printhelp) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting at most 6 component of the Pressure BC, got %i in %s\n",nopt,__FUNCT__);
    for (i = 0; i < 6; i++) ctx->BCpres[i]=buffer[i];

    ierr = PetscOptionsRealArray("-BCtheta", "\n\tTemperature at Boundaries.\n\t (TX0,TX1,TY0,TY1,TZ0,TZ1) negative value if natural BC","",buffer,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt > 6 && !ctx->printhelp) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting at most 6 component of the Temperature BC, got %i in %s\n",nopt,__FUNCT__);
    for (i = 0; i < 6; i++) ctx->BCtheta[i]=buffer[i];

    nopt = 3;
    for (i = 0; i < 3; i++) ctx->SrcLoc[i] = 99999;
    ierr = PetscOptionsIntArray("-SrcLoc","\n\t location of source point","",ctx->SrcLoc,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != 3 && nopt != 0) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting at most 3 component of SrcLoc, got %i in %s\n",nopt,__FUNCT__);
    ctx->SrcRate          = 0.0;
    ierr                  = PetscOptionsReal("-SrcRate","\n\tStrength of the source in kg/s","",ctx->SrcRate,&ctx->SrcRate,PETSC_NULL);CHKERRQ(ierr);
    ctx->hasCrackPressure = PETSC_FALSE;
    ierr                  = PetscOptionsBool("-pressurize","\n\tPressurize cracks","",ctx->hasCrackPressure,&ctx->hasCrackPressure,PETSC_NULL);CHKERRQ(ierr);

    ctx->numPennyCracks = 0;
    ierr                = PetscOptionsInt("-npc","\n\tNumber of penny-shaped cracks to insert","",ctx->numPennyCracks,&ctx->numPennyCracks,PETSC_NULL);CHKERRQ(ierr);
    ierr                = PetscMalloc(ctx->numPennyCracks*sizeof(VFPennyCrack),&ctx->pennycrack);CHKERRQ(ierr);
    for (i = 0; i < ctx->numPennyCracks; i++) {
      ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"pc%d_",i);CHKERRQ(ierr);
      ierr = VFPennyCrackCreate(&(ctx->pennycrack[i]));CHKERRQ(ierr);
      ierr = VFPennyCrackGet(prefix, &(ctx->pennycrack[i]));CHKERRQ(ierr);
      if (ctx->verbose > 0) {
        ierr = VFPennyCrackView(&(ctx->pennycrack[i]),PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      }
    }

    ctx->numRectangularCracks = 0;
    ierr                      = PetscOptionsInt("-nrc","\n\tNumber of rectangular cracks to insert","",ctx->numRectangularCracks,&ctx->numRectangularCracks,PETSC_NULL);CHKERRQ(ierr);
    ierr                      = PetscMalloc(ctx->numRectangularCracks*sizeof(VFRectangularCrack),&ctx->rectangularcrack);CHKERRQ(ierr);
    for (i = 0; i < ctx->numRectangularCracks; i++) {
      ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"rc%d_",i);CHKERRQ(ierr);
      ierr = VFRectangularCrackCreate(&(ctx->rectangularcrack[i]));CHKERRQ(ierr);
      ierr = VFRectangularCrackGet(prefix, &(ctx->rectangularcrack[i]));CHKERRQ(ierr);
      if (ctx->verbose > 0) {
        ierr = VFRectangularCrackView(&(ctx->rectangularcrack[i]),PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      }
    }

    ctx->numWells = 0;
    ierr          = PetscOptionsInt("-nw","\n\tNumber of wells to insert","",ctx->numWells,&ctx->numWells,PETSC_NULL);CHKERRQ(ierr);
    ierr          = PetscMalloc(ctx->numWells*sizeof(VFWell),&ctx->well);CHKERRQ(ierr);
    for (i = 0; i < ctx->numWells; i++) {
      ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"w%d_",i);CHKERRQ(ierr);
      ierr = VFWellCreate(&(ctx->well[i]));CHKERRQ(ierr);
      ierr = VFWellGet(prefix, &(ctx->well[i]));CHKERRQ(ierr);
      if (ctx->verbose > 0) {
        ierr = VFWellView(&(ctx->well[i]),PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      }
    }
    ierr = PetscFree(buffer);CHKERRQ(ierr);
  }
  ierr           = PetscOptionsEnd();CHKERRQ(ierr);
  ctx->timestep  = 1;
  ctx->timevalue = 0.;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFGeometryInitialize"
/*
 VFGeometryInitialize: Creates DA, and coordinates, and other geometric informations

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFGeometryInitialize(VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscViewer    viewer,h5viewer;
  char           filename[FILENAME_MAX];
  PetscReal      BBmin[3],BBmax[3];
  PetscReal      *X,*Y,*Z;
  PetscReal      ****coords_array;
  PetscInt       xs,xm,ys,ym,zs,zm;
  int            i,j,k;
  int            nval;
  PetscInt       *n,nx,ny,nz;
  PetscReal      *l,lx,ly,lz;
  const PetscInt *lx1,*ly1,*lz1;
  PetscInt       x_nprocs,y_nprocs,z_nprocs,*olx,*oly,*olz;
  PetscBool      flg;

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\n\nVF-Chevron: geometry options:","");CHKERRQ(ierr);
  {
    ierr = PetscMalloc(3 * sizeof(PetscReal),&l);CHKERRQ(ierr);
    for (i = 0; i < 3; i++) l[i] = 1.;
    nval = 3;
    ierr = PetscOptionsRealArray("-l","\n\tDomain dimensions (default 1.), comma separated","",l,&nval,&flg);CHKERRQ(ierr);
    if (nval != 3 && nval != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting 3 values for option %s, got only %i in %s\n",n,"-l",nval,__FUNCT__);
    lx   = l[0];
    ly   = l[1];
    lz   = l[2];
    ierr = PetscFree(l);CHKERRQ(ierr);

    ierr = PetscMalloc(3 * sizeof(PetscInt),&n);CHKERRQ(ierr);
    for (i = 0; i < 3; i++) n[i] = 11;
    nval = 3;
    ierr = PetscOptionsIntArray("-n","\n\tnumber of grid points (default 11), comma separated","",n,&nval,&flg);CHKERRQ(ierr);
    if (nval != 3 && nval != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting 3 values for option %s, got only %i in %s\n",n,"-n",nval,__FUNCT__);
    nx   = n[0];
    ny   = n[1];
    nz   = n[2];
    ierr = PetscFree(n);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  if (ctx->printhelp) PetscFunctionReturn(0);

  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
                      DMDA_STENCIL_BOX,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,
                      PETSC_NULL,PETSC_NULL,PETSC_NULL,&ctx->daVect);CHKERRQ(ierr);
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
                      DMDA_STENCIL_BOX,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,
                      PETSC_NULL,PETSC_NULL,PETSC_NULL,&ctx->daScal);CHKERRQ(ierr);
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
                      DMDA_STENCIL_BOX,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,4,1,
                      PETSC_NULL,PETSC_NULL,PETSC_NULL,&ctx->daFlow);CHKERRQ(ierr);



  ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,&x_nprocs,&y_nprocs,&z_nprocs,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRanges(ctx->daScal,&lx1,&ly1,&lz1);CHKERRQ(ierr);
  ierr = PetscMalloc(x_nprocs*sizeof(*olx),&olx);CHKERRQ(ierr);
  ierr = PetscMalloc(y_nprocs*sizeof(*oly),&oly);CHKERRQ(ierr);
  ierr = PetscMalloc(z_nprocs*sizeof(*olz),&olz);CHKERRQ(ierr);

  ierr = PetscMemcpy(olx,lx1,x_nprocs*sizeof(*olx));CHKERRQ(ierr);
  ierr = PetscMemcpy(oly,ly1,y_nprocs*sizeof(*oly));CHKERRQ(ierr);
  ierr = PetscMemcpy(olz,lz1,z_nprocs*sizeof(*olz));CHKERRQ(ierr);

  olx[x_nprocs-1]--;
  oly[y_nprocs-1]--;
  olz[z_nprocs-1]--;
  /*
   Finite elements should never use ghost cells, so the stencil for the cell DA is 0
   */
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
                      DMDA_STENCIL_BOX,nx-1,ny-1,nz-1,x_nprocs,y_nprocs,z_nprocs,1,0,
                      olx,oly,olz,&ctx->daScalCell);CHKERRQ(ierr);


  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
                      DMDA_STENCIL_BOX,nx-1,ny-1,nz-1,x_nprocs,y_nprocs,z_nprocs,6,1,
                      olx,oly,olz,&ctx->daVFperm);CHKERRQ(ierr);

  /*
   ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
   DMDA_STENCIL_BOX,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,6,1,
   PETSC_NULL,PETSC_NULL,PETSC_NULL,&ctx->daVFperm);CHKERRQ(ierr);
   */


  ierr = CartFE_Init();CHKERRQ(ierr);
  ierr = CartFE_Element3DCreate(&ctx->e3D);CHKERRQ(ierr);
  /*
   Constructs coordinates Vec
   */
  ierr = DMCreateGlobalVector(ctx->daVect,&ctx->coordinates);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->coordinates,"Coordinates");CHKERRQ(ierr);

  /*
   Construct coordinates from the arrays of cell sizes
   */
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = PetscMalloc3(nx,PetscReal,&X,ny,PetscReal,&Y,nz,PetscReal,&Z);CHKERRQ(ierr);
  Z[0] = 0.;
  Y[0] = 0.;
  X[0] = 0.;
  for (k = 1; k < nz; k++) Z[k] = k * lz / (nz-1.);
  for (j = 1; j < ny; j++) Y[j] = j * ly / (ny-1.);
  for (i = 1; i < nx; i++) X[i] = i * lx / (nx-1.);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++)
      for (i = xs; i < xs + xm; i++) {
        coords_array[k][j][i][2] = Z[k];
        coords_array[k][j][i][1] = Y[j];
        coords_array[k][j][i][0] = X[i];
      }
  }
  ierr = PetscFree3(X,Y,Z);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = DMDASetCoordinates(ctx->daVect,ctx->coordinates);CHKERRQ(ierr);
  ierr = DMDASetCoordinates(ctx->daScal,ctx->coordinates);CHKERRQ(ierr);
  switch (ctx->fileformat) {
  case FILEFORMAT_BIN:
    /*
     As of version 3.1, there is a bug in petsc preventing to save 2 DA with different number of degrees of freedoms
     per node in a single file.
     Even coordinates are part of the DA and do not need to be saved separately, it looks like the coordinate vector
     obtained from DMDAGetCoordinates cannot be saved properly in an hdf5 file, so we save the coordinate vector anyway
     */
    ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.bin",ctx->prefix);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = DMView(ctx->daScal,viewer);CHKERRQ(ierr);
    ierr = VecView(ctx->coordinates,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    break;
#if defined(PETSC_HAVE_HDF5)
  case FILEFORMAT_HDF5:
    /*
     Write headers in multistep XDMF file
     */
    ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.xmf",ctx->prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&ctx->XDMFviewer);
    ierr = XDMFmultistepInitialize(ctx->XDMFviewer);CHKERRQ(ierr);

    /*
     Save cordinate in main hdf5 file
     */
    ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.h5",ctx->prefix);CHKERRQ(ierr);
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&h5viewer);CHKERRQ(ierr);
    ierr = VecView(ctx->coordinates,h5viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&h5viewer);CHKERRQ(ierr);
#endif
  }

  /*
   Replace the coordinate vector defined on PETSC_COMM_WORLD with ghosted coordinate vectors on PETSC_COMM_SELF,
   so that we can query coordinates of ghost points
   */
  ierr = VecDestroy(&ctx->coordinates);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(ctx->daScal,&ctx->coordinates);CHKERRQ(ierr);


  if (ctx->verbose > 0) {
    /*
     Get bounding box from petsc DA
     */
    ierr = DMDAGetBoundingBox(ctx->daVect,BBmin,BBmax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reservoir bounding box: (%g,%g) x (%g,%g) x (%g,%g)\n",BBmin[0],BBmax[0],BBmin[1],BBmax[1],BBmin[2],BBmax[2]);CHKERRQ(ierr);
    ierr = DMView(ctx->daVect,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = DMView(ctx->daScal,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = VFLayerInit(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFPropGet"
/*
 VFPropGet

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
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
    ierr             = PetscOptionsReal("-epsilon","\n\tVariational fracture regularization length (start from 4x cell size)","",vfprop->epsilon,&vfprop->epsilon,PETSC_NULL);CHKERRQ(ierr);
    ierr             = PetscOptionsReal("-eta","\n\tArtificial stiffness of cracks ","",vfprop->eta,&vfprop->eta,PETSC_NULL);CHKERRQ(ierr);
    ierr             = PetscOptionsReal("-irrevtol","\n\tThreshold on v below which fracture irreversibility is enforced","",vfprop->irrevtol,&vfprop->irrevtol,PETSC_NULL);CHKERRQ(ierr);
    vfprop->permmax  = 5.;
    ierr             = PetscOptionsReal("-permmax","\n\tPermeability multiplier of cracks (achieved at  v=0.)","",vfprop->permmax,&vfprop->permmax,PETSC_NULL);CHKERRQ(ierr);
    vfprop->atCv     = .5;
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VFMatPropGet"
/*
 VFMatPropGet: get all material properties.
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFMatPropGet(VFMatProp *matprop,PetscInt n)
{
  PetscErrorCode ierr;
  PetscReal      *prop;
  PetscInt       nopt;
  PetscReal      E     = 1.;
  PetscReal      nu    = 0.;
  PetscReal      alpha = 1.e-5;
  PetscReal      Gc    = 1.;
  PetscReal      beta  = 1.e-4;            /* Normalized by E (MPa) assuming E=10,000 MPa */
  PetscReal      Cp    = 1.;
  PetscReal      rho   = 1.;
  int            i;

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\n\nVF: material properties:","");CHKERRQ(ierr);
  {
    ierr = PetscMalloc(n * sizeof(PetscReal),&prop);CHKERRQ(ierr);
    nopt = n;
    for (i = 0; i < n; i++) prop[i] = E;
    ierr = PetscOptionsRealArray("-E","\n\tComma separated list of Young\'s modulii","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-E",nopt,__FUNCT__);
    for (i=0; i< n; i++) matprop[i].E = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = nu;
    ierr = PetscOptionsRealArray("-nu","\n\tComma separated list of Poisson\'s ratio","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-nu",nopt,__FUNCT__);
    for (i=0; i< n; i++) matprop[i].nu = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = alpha;
    ierr = PetscOptionsRealArray("-alpha","\n\tComma separated list of linear thermal expansion coefficients","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-alpha",nopt,__FUNCT__);
    for (i=0; i< n; i++) matprop[i].alpha = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = Gc;
    ierr = PetscOptionsRealArray("-Gc","\n\tComma separated list of fracture toughness (not critical SIF!)","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-Gc",nopt,__FUNCT__);
    for (i=0; i< n; i++) matprop[i].Gc = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = beta;
    ierr = PetscOptionsRealArray("-beta","\n\tComma separated list of Biot constants","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-beta",nopt,__FUNCT__);
    for (i=0; i< n; i++) matprop[i].beta = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = rho;
    ierr = PetscOptionsRealArray("-rho","\n\tComma separated list of reservoir densities","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-rho",nopt,__FUNCT__);
    for (i=0; i< n; i++) matprop[i].rho = prop[i];

    nopt = n;
    for (i = 0; i < n; i++) prop[i] = Cp;
    ierr = PetscOptionsRealArray("-Cp","\n\tComma separated list of specific heat capacities","",prop,&nopt,PETSC_NULL);CHKERRQ(ierr);
    if (nopt != n && nopt != 0) SETERRQ4(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Expecting %i values for option %s, got only %i in %s\n",n,"-Cp",nopt,__FUNCT__);
    for (i=0; i< n; i++) matprop[i].Cp = prop[i];
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  ierr = PetscFree(prop);CHKERRQ(ierr);

  for (i=0; i < n; i++) {
    matprop[i].lambda = matprop[i].E * matprop[i].nu / (1. + matprop[i].nu) / (1. - 2. * matprop[i].nu);
    matprop[i].mu     = matprop[i].E / (1. + matprop[i].nu) * .5;
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFMatPropFieldsInitialize"
/*
 VFMatPropFieldsInitialize: Creates material property fields
 */
extern PetscErrorCode VFMatPropFieldsInitialize(VFCtx *ctx, VFMatProp *matprop)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = DMCreateGlobalVector(ctx->daScalCell,&matprop->VecGc);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) matprop->VecGc,"Critical Energy Release Rate");CHKERRQ(ierr);
  ierr = VecSet(matprop->VecGc,matprop->Gc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFResPropGet"
/*
 VFResPropGet

 Change it to read properties from input file later
 */
extern PetscErrorCode VFResPropGet(VFResProp *resprop)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"","");CHKERRQ(ierr);
  {
    resprop->perm      = 1.e-1;                    /* Multiply by 1e12 because pressure unit in MPa, viscosity unit in cp, and density is specific density */
    resprop->por       = 0.2;                      /* fraction */
    resprop->Pinit     = 1.;                      /* MPa */
    resprop->Tinit     = 200.;                     /* Celsius */
    resprop->relk      = 1.0;                      /* fraction */
    resprop->visc      = 1.0;                      /* cp */
    resprop->fdens     = 1.0;                      /* specific density */
    resprop->rock_comp = 0.5e-4;                          /* Rock Compressibility */
    resprop->wat_comp  =5e-4;                         /* Water Compressibility */
    resprop->TCond_X   = 0.6;                       /* Water thermal conductivity in W/m-K */
    resprop->TCond_Y   = resprop->TCond_X;
    resprop->TCond_Z   = resprop->TCond_X;
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFFieldsInitialize"
/*
 VFFieldsInitialize: Creates fields

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFFieldsInitialize(VFCtx *ctx,VFFields *fields)
{

  PetscErrorCode ierr;
  PetscInt       c;

  PetscFunctionBegin;
  fields->numfields = 9;
  ierr              = DMCreateGlobalVector(ctx->daVect,&fields->U);CHKERRQ(ierr);
  ierr              = PetscObjectSetName((PetscObject) fields->U,"Displacement");CHKERRQ(ierr);
  ierr              = VecSet(fields->U,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daVect,&fields->BCU);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->BCU,"Displacement");CHKERRQ(ierr);
  ierr = VecSet(fields->BCU,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->V);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->V,"Fracture");CHKERRQ(ierr);
  ierr = VecSet(fields->V,1.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->VIrrev);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->VIrrev,"VIrrev");CHKERRQ(ierr);
  ierr = VecSet(fields->VIrrev,1.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->theta);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->theta,"Temperature");CHKERRQ(ierr);
  ierr = VecSet(fields->theta,ctx->resprop.Tinit);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->thetaRef);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->thetaRef,"Reference Temperature");CHKERRQ(ierr);
  ierr = VecSet(fields->thetaRef,ctx->resprop.Tinit);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->pressure);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->pressure,"Pressure");CHKERRQ(ierr);
  ierr = VecSet(fields->pressure,ctx->resprop.Pinit);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->pressureRef);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->pressureRef,"Reference Pressure");CHKERRQ(ierr);
  ierr = VecSet(fields->pressureRef,ctx->resprop.Pinit);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScalCell,&fields->pmult);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->pmult,"PermeabilityMultiplier");CHKERRQ(ierr);
  ierr = VecSet(fields->pmult,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daFlow,&fields->VelnPress);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->VelnPress,"Velocity and Pressure");CHKERRQ(ierr);
  ierr = VecSet(fields->VelnPress,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daVFperm,&fields->vfperm);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->vfperm,"Permeability_V-field");CHKERRQ(ierr);
  ierr = VecSet(fields->vfperm,1.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daVect,&fields->velocity);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->velocity,"Fluid Velocity");CHKERRQ(ierr);
  ierr = VecSet(fields->velocity,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->VolCrackOpening);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->VolCrackOpening,"Volumetric Crack Opening");CHKERRQ(ierr);
  ierr = VecSet(fields->VolCrackOpening,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->VolLeakOffRate);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->VolLeakOffRate,"Volumetric Leakoff Rate");CHKERRQ(ierr);
  ierr = VecSet(fields->VolLeakOffRate,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daFlow,&fields->FlowBCArray);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->FlowBCArray,"Velocity and Pressure Boundary Values");CHKERRQ(ierr);
  ierr = VecSet(fields->FlowBCArray,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&ctx->Source);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->Source,"Source Term");CHKERRQ(ierr);
  ierr = VecSet(ctx->Source,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->PresBCArray);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->PresBCArray,"Pressure Boundary Values");CHKERRQ(ierr);
  ierr = VecSet(fields->PresBCArray,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&ctx->PresBCArray);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->PresBCArray,"Pressure Boundary Values in ctx");CHKERRQ(ierr);
  ierr = VecSet(ctx->PresBCArray,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daVect,&ctx->VelBCArray);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->VelBCArray,"Velocity boundary Values in ctx");CHKERRQ(ierr);
  ierr = VecSet(ctx->VelBCArray,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&ctx->HeatSource);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->HeatSource,"Heat Source Term");CHKERRQ(ierr);
  ierr = VecSet(ctx->HeatSource,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&ctx->HeatFluxBCArray);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->HeatFluxBCArray,"Heat Flux Boundary Term");CHKERRQ(ierr);
  ierr = VecSet(ctx->HeatFluxBCArray,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daVFperm,&ctx->Cond);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->Cond,"Thermal Conductivity of Rock");CHKERRQ(ierr);
  ierr = VecSet(ctx->Cond,1.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&ctx->TBCArray);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->TBCArray,"Temperature Boundary Term");CHKERRQ(ierr);
  ierr = VecSet(ctx->TBCArray,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daVect,&fields->fracvelocity);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->fracvelocity,"Fracture Fluid Velocity");CHKERRQ(ierr);
  ierr = VecSet(fields->fracvelocity,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daFlow,&fields->fracVelnPress);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->fracVelnPress,"Fracture Velocity and Pressure");CHKERRQ(ierr);
  ierr = VecSet(fields->fracVelnPress,0.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&fields->fracpressure);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fields->fracpressure,"Fracture Pressure");CHKERRQ(ierr);
  ierr = VecSet(fields->fracpressure,0.);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daVect,&ctx->FracVelBCArray);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->FracVelBCArray,"Fracture Velocity boundary Values");CHKERRQ(ierr);
  ierr = VecSet(ctx->FracVelBCArray,0.0);CHKERRQ(ierr);


  /*
   Create optional penny-shaped and rectangular cracks
   */
  ierr = VecSet(fields->V,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields->VIrrev,1.0);CHKERRQ(ierr);
  for (c = 0; c < ctx->numPennyCracks; c++) {
    ierr = VFPennyCrackBuildVAT2(fields->V,&(ctx->pennycrack[c]),ctx);CHKERRQ(ierr);
    ierr = VecPointwiseMin(fields->VIrrev,fields->V,fields->VIrrev);CHKERRQ(ierr);
  }
  for (c = 0; c < ctx->numRectangularCracks; c++) {
    ierr = VFRectangularCrackBuildVAT2(fields->V,&(ctx->rectangularcrack[c]),ctx);CHKERRQ(ierr);
    ierr = VecPointwiseMin(fields->VIrrev,fields->V,fields->VIrrev);CHKERRQ(ierr);
  }
  /*
   Create optional wells
   */
  for (c = 0; c < ctx->numWells; c++) {
    ierr = VFWellBuildVAT2(fields->V,&(ctx->well[c]),ctx);CHKERRQ(ierr);
    ierr = VecPointwiseMin(fields->VIrrev,fields->V,fields->VIrrev);CHKERRQ(ierr);
  }
  ierr        = VecCopy(fields->VIrrev,fields->V);CHKERRQ(ierr);
  ctx->fields = fields;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFBCInitialize"
/*
 VFBCInitialize: Creates and initialize boundary condition data structures

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
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

  ierr = BCPInit(&ctx->bcP[0],ctx);CHKERRQ(ierr);
  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"BCP:\n");CHKERRQ(ierr);
    ierr = BCView(&ctx->bcP[0],PETSC_VIEWER_STDOUT_WORLD,1);CHKERRQ(ierr);
  }

  ierr = BCQInit(&ctx->bcQ[0],ctx);CHKERRQ(ierr);
  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"BCQ:\n");CHKERRQ(ierr);
    ierr = BCView(&ctx->bcQ[0],PETSC_VIEWER_STDOUT_WORLD,1);CHKERRQ(ierr);
  }

  ierr = BCTInit(&ctx->bcT[0],ctx);CHKERRQ(ierr);
  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"BCT:\n");CHKERRQ(ierr);
    ierr = BCView(&ctx->bcT[0],PETSC_VIEWER_STDOUT_WORLD,1);CHKERRQ(ierr);
  }

  ierr = BCTInit(&ctx->bcQT[0],ctx);CHKERRQ(ierr);
  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"BCQT:\n");CHKERRQ(ierr);
    ierr = BCView(&ctx->bcQT[0],PETSC_VIEWER_STDOUT_WORLD,1);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFSolversInitialize"
/*
 VFSolversInitialize: Creates matrices, RHS, solvers

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 Keita Yoshioka yoshk@chevron.com
 */
extern PetscErrorCode VFSolversInitialize(VFCtx *ctx)
{
  PetscErrorCode ierr;
  KSP            kspU,kspV;
  PC             pcU,pcV;
  PetscReal      atol,rtol,stol;
  PetscInt       maxit,maxf;
  Vec            lbV,ubV;
  Mat            JacV;
  Vec            residualV;
  PetscFunctionBegin;
  /*
   U solver initialization
   */
  ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesU);CHKERRQ(ierr);
  ierr = SNESSetDM(ctx->snesU,ctx->daVect);CHKERRQ(ierr);
  ierr = SNESSetOptionsPrefix(ctx->snesU,"U_");CHKERRQ(ierr);
  ierr = SNESSetType(ctx->snesU,SNESKSPONLY);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(ctx->snesU);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daVect,&ctx->UResidual);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->UResidual,"U_Residual");CHKERRQ(ierr);

  /*
   KU and RHSU will go when we update the residual evaluation
   */
  ierr = DMCreateGlobalVector(ctx->daVect,&ctx->RHSU);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->RHSU,"RHSU");CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daVect,PETSC_NULL,&ctx->KU);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daVect,PETSC_NULL,&ctx->JacU);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->KU,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);

  ierr = SNESGetKSP(ctx->snesU,&kspU);CHKERRQ(ierr);
  ierr = KSPSetTolerances(kspU,1.e-8,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetType(kspU,KSPBCGSL);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspU);CHKERRQ(ierr);
  ierr = KSPGetPC(kspU,&pcU);CHKERRQ(ierr);
  ierr = PCSetType(pcU,PCBJACOBI);CHKERRQ(ierr);
  ierr = PCSetFromOptions(pcU);CHKERRQ(ierr);
  /*
   V solver initialization
   */
  ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesV);CHKERRQ(ierr);
  ierr = SNESSetDM(ctx->snesV,ctx->daScal);CHKERRQ(ierr);
  ierr = SNESSetOptionsPrefix(ctx->snesV,"V_");CHKERRQ(ierr);
  if (ctx->verbose > 1) {
    ierr = SNESMonitorSet(ctx->snesV,VF_VSNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = SNESSetFromOptions(ctx->snesV);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&lbV);CHKERRQ(ierr);
  ierr = VecSet(lbV,0.0);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->daScal,&ubV);CHKERRQ(ierr);
  ierr = VecSet(ubV,1.0);CHKERRQ(ierr);
  ierr = SNESVISetVariableBounds(ctx->snesV,lbV,ubV);CHKERRQ(ierr);
  ierr = SNESGetTolerances(ctx->snesV,&atol,&rtol,&stol,&maxit,&maxf);CHKERRQ(ierr);
  ierr = SNESSetTolerances(ctx->snesV,1e-8,rtol,stol,maxit,maxf);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&residualV);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daScal,PETSC_NULL,&JacV);CHKERRQ(ierr);
  ierr = MatSetOption(JacV,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = SNESSetFunction(ctx->snesV,residualV,VF_VResidual,ctx);CHKERRQ(ierr);
  ierr = SNESSetJacobian(ctx->snesV,JacV,JacV,VF_VIJacobian,ctx);CHKERRQ(ierr);

  ierr = SNESGetKSP(ctx->snesV,&kspV);CHKERRQ(ierr);
  ierr = KSPSetTolerances(kspV,1.e-8,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetType(kspV,KSPBCGSL);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspV);CHKERRQ(ierr);
  ierr = KSPGetPC(kspV,&pcV);CHKERRQ(ierr);
  ierr = PCSetType(pcV,PCBJACOBI);CHKERRQ(ierr);
  ierr = PCSetFromOptions(pcV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFLayerInit"
/*
 VFLayerInit: find the horizontal layer associated to a cell in the most stupid and unoptimized way

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu

 */
extern PetscErrorCode VFLayerInit(VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       ek,l;
  PetscReal      ****coords_array;
  PetscInt       xs,xm,ys,ym,zs,zm;
  PetscInt       nz;

  PetscFunctionBegin;
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx->daVect,PETSC_NULL,PETSC_NULL,PETSC_NULL,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscMalloc(nz * sizeof(PetscInt),&ctx->layer);CHKERRQ(ierr);

  for (l = 0; l < nz; l++) ctx->layer[l] = 0;
  for (ek = zs; ek < zs+zm; ek++)
    for (l = 0; l < ctx->nlayer; l++)
      if (coords_array[ek][ys][xs][2] > ctx->layersep[l]) ctx->layer[ek] = l;
  if (ctx->verbose > 0)
    for (ek = 0; ek < nz; ek++) {
      if (ek >= zs && ek < zs+zm) {
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ctx->layer[%i]=%i (depth=%g)\n",ek,ctx->layer[ek],coords_array[ek][ys][xs][2]);CHKERRQ(ierr);
      } else {
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"layer %i inactive\n",ek);CHKERRQ(ierr);
      }
      ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    }
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFTimeStepPrepare"
/*
 VFTimeStepPrepare: Prepare for a new time step:
 - Update VIrrev
 - Read boundary displacement from files if necessary
 - Set boundary values of U and V

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFTimeStepPrepare(VFCtx *ctx,VFFields *fields)
{
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

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFElasticityTimeStep"
/*
 VFElasticityTimeStep

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFElasticityTimeStep(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /*
   Assembly, solve
   */
  ierr              = VF_StepU(fields,ctx);CHKERRQ(ierr);
  ctx->ElasticEnergy=0;
  ctx->InsituWork   =0;
  ctx->PressureWork = 0.;
  ierr              = VF_UEnergy3D(&ctx->ElasticEnergy,&ctx->InsituWork,&ctx->PressureWork,fields,ctx);CHKERRQ(ierr);
  ctx->TotalEnergy  = ctx->ElasticEnergy - ctx->InsituWork - ctx->PressureWork;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFFractureTimeStep"
/*
 VFFractureTimeStep

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu

 */
extern PetscErrorCode VFFractureTimeStep(VFCtx *ctx,VFFields *fields)
{
  PetscInt       altminit=1;
  Vec            Vold;
  PetscReal      errV=1e+10;
  PetscErrorCode ierr;

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

  ctx->ElasticEnergy = 0.;
  ctx->InsituWork    = 0.;
  ierr               = VF_UEnergy3D(&ctx->ElasticEnergy,&ctx->InsituWork,&ctx->PressureWork,fields,ctx);CHKERRQ(ierr);
  ctx->SurfaceEnergy = 0.;
  ierr               = VF_VEnergy3D(&ctx->SurfaceEnergy,fields,ctx);CHKERRQ(ierr);
  ctx->TotalEnergy   = ctx->ElasticEnergy + ctx->SurfaceEnergy - ctx->InsituWork - ctx->PressureWork;
  ierr               = VecDestroy(&Vold);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFFinalize"
/*
 VFFinalize

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFFinalize(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;
  char           filename[FILENAME_MAX];
  PetscViewer    optionsviewer;
  PetscInt       nopts;

  PetscFunctionBegin;

  ierr = PetscFree(ctx->matprop);CHKERRQ(ierr);
  ierr = PetscFree(ctx->layer);CHKERRQ(ierr);
  ierr = PetscFree(ctx->layersep);CHKERRQ(ierr);
  ierr = PetscFree(ctx->pennycrack);CHKERRQ(ierr);
  ierr = PetscFree(ctx->rectangularcrack);CHKERRQ(ierr);
  ierr = PetscFree(ctx->well);CHKERRQ(ierr);

  ierr = DMDestroy(&ctx->daVect);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->daScal);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->daVFperm);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->daFlow);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->daScalCell);CHKERRQ(ierr);

  ierr = MatDestroy(&ctx->KU);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->RHSU);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->UResidual);CHKERRQ(ierr);
  ierr = SNESDestroy(&ctx->snesU);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->JacU);CHKERRQ(ierr);


  ierr = SNESDestroy(&ctx->snesV);CHKERRQ(ierr);

  ierr = VecDestroy(&fields->VelnPress);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->vfperm);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->VolCrackOpening);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->VolLeakOffRate);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->FlowBCArray);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->Source);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->PresBCArray);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->PresBCArray);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->VelBCArray);CHKERRQ(ierr);

  ierr = VecDestroy(&ctx->HeatSource);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->HeatFluxBCArray);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->Cond);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->TBCArray);CHKERRQ(ierr);

  ierr = VecDestroy(&fields->U);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->BCU);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->V);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->VIrrev);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->theta);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->thetaRef);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->pressure);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->velocity);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->pressureRef);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->pmult);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&ctx->energyviewer);CHKERRQ(ierr);

  ierr = VecDestroy(&fields->fracpressure);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->fracvelocity);CHKERRQ(ierr);
  ierr = VecDestroy(&fields->fracVelnPress);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->FracVelBCArray);CHKERRQ(ierr);


  ierr = VF_HeatSolverFinalize(ctx,fields);CHKERRQ(ierr);
  ierr = FlowSolverFinalize(ctx,fields);CHKERRQ(ierr);
  /*
   Close the xdmf multi-step file
   */
#if defined(PETSC_HAVE_HDF5)
  if (ctx->fileformat == FILEFORMAT_HDF5) {
    ierr = XDMFuniformgridFinalize(ctx->XDMFviewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&ctx->XDMFviewer);CHKERRQ(ierr);
  }
#endif
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx->prefix);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&optionsviewer);CHKERRQ(ierr);
  ierr = PetscLogView(optionsviewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&optionsviewer);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.py",ctx->prefix);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&optionsviewer);CHKERRQ(ierr);
  ierr = PetscLogViewPython(optionsviewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&optionsviewer);CHKERRQ(ierr);
  ierr = PetscOptionsAllUsed(&nopts);CHKERRQ(ierr);
  if (nopts > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nWARNING: \nSome options were unused. Check the command line for typos.\n");CHKERRQ(ierr);
    ierr = PetscOptionsLeft();CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FieldsH5Write"
/*
 FieldsH5Write: Export all fields in HDF5 format using PETSc HDF5 viewer. Also write an XDMF container

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode FieldsH5Write(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;
  char           H5filename[FILENAME_MAX],H5coordfilename[FILENAME_MAX],XDMFfilename[FILENAME_MAX];
  PetscViewer    H5Viewer,XDMFViewer;
  PetscInt       nx, ny, nz;

  PetscFunctionBegin;
  /*
   Write the hdf5 files
  */
#if defined(PETSC_HAVE_HDF5)
  ierr = DMDAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscSNPrintf(H5filename,FILENAME_MAX,"%s.%.5i.h5",ctx->prefix,ctx->timestep);CHKERRQ(ierr);
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,H5filename,FILE_MODE_WRITE,&H5Viewer);
  ierr = VecView(fields->U,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->velocity,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->fracvelocity,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->vfperm,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->V,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->pmult,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->theta,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->pressure,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->fracpressure,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->VolCrackOpening,H5Viewer);CHKERRQ(ierr);
  ierr = VecView(fields->VolLeakOffRate,H5Viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&H5Viewer);CHKERRQ(ierr);

  /*
   Write the XDMF Description
   */
  ierr = PetscSNPrintf(XDMFfilename,FILENAME_MAX,"%s.%.5i.xmf",ctx->prefix,ctx->timestep);CHKERRQ(ierr);
  ierr = PetscSNPrintf(H5coordfilename,FILENAME_MAX,"%s.h5",ctx->prefix);CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,XDMFfilename,&XDMFViewer);CHKERRQ(ierr);
  ierr = XDMFuniformgridInitialize(XDMFViewer,ctx->timevalue,H5filename);CHKERRQ(ierr);
  ierr = XDMFtopologyAdd (XDMFViewer,nx,ny,nz,H5coordfilename,"Coordinates");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,3,"Vector","Node",H5filename,"Displacement");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,3,"Vector","Node",H5filename,"Fluid Velocity");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,3,"Vector","Node",H5filename,"Fracture Fluid Velocity");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx-1,ny-1,nz-1,6,"Tensor6","Cell",H5filename,"Permeability_V-field");CHKERRQ(ierr);       /*Cell*/
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,1,"Scalar","Node",H5filename,"Fracture");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx-1,ny-1,nz-1,1,"Scalar","Cell",H5filename,"PermeabilityMultiplier");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,1,"Scalar","Node",H5filename,"Temperature");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,1,"Scalar","Node",H5filename,"Pressure");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,1,"Scalar","Node",H5filename,"Fracture Pressure");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,1,"Scalar","Node",H5filename,"Volumetric Crack Opening");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFViewer,nx,ny,nz,1,"Scalar","Node",H5filename,"Volumetric Leakoff Rate");CHKERRQ(ierr);
  ierr = XDMFuniformgridFinalize(XDMFViewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&XDMFViewer);CHKERRQ(ierr);
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

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
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
  ierr = VecView(fields->velocity,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->V,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->pmult,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->theta,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->pressure,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->VolCrackOpening,viewer);CHKERRQ(ierr);
  ierr = VecView(fields->VolLeakOffRate,viewer);CHKERRQ(ierr);
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
  PetscReal      pmult,  v;
  PetscInt       xs,xm;
  PetscInt       ys,ym;
  PetscInt       zs,zm;
  PetscInt       i,j,k;
  PetscReal      ***v_array,***pmult_array;

  PetscFunctionBegin;

  ierr = DMDAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,Pmult,&pmult_array);CHKERRQ(ierr);
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++)
      for (i = xs; i < xs + xm; i++) {
        v = v_array[k][j][i];
        if (v <= 0) {
          pmult = vfprop->permmax;
        } else {
          if (v <= 1.) {
            pmult = vfprop->permmax - (vfprop->permmax - 1.0) * v;
          } else {
            pmult = 1.;
          }
        }
        pmult_array[k][j][i] = pmult;
      }
  }
  ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,Pmult,&pmult_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
