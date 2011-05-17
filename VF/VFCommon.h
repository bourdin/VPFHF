#ifndef VFCOMMON_H
#define VFCOMMON_H
/*
  VFCommon.h
  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
static const char banner[] = "\n\nVF:\nNumerical implementation of the variational approach to fracture.\n(c) 2010-2011 Blaise Bourdin, Louisiana State University. bourdin@lsu.edu\n\n";

typedef enum { 
  SYMXY,
  SYMXY_DISP,
  SYMX,
  SYMX_DISP,
  SYMY,
  SYMY_DISP,
  NOSYM,
  NOSYM_DISP,
  TEST_CLAMPEDX0,
  TEST_CLAMPEDX1,
  TEST_CLAMPEDX0X1,
  TEST_CLAMPEDY0,
  TEST_CLAMPEDY1,
  TEST_CLAMPEDY0Y1,
  TEST_CLAMPEDZ0,
  TEST_CLAMPEDZ1,
  TEST_CLAMPEDZ0Z1,
  TEST_MANUAL
} VFPreset;
static const char *VFPresetName[] = {"SYMXY","SYMXY_DISP","SYMX","SYMX_DISP","SYMY","SYMY_DISP",
                                     "NOSYM","NOSYM_DISP",
                                     "TEST_CLAMPEDX0","TEST_CLAMPEDX1","TEST_CLAMPEDX0X1",
                                     "TEST_CLAMPEDY0","TEST_CLAMPEDY1","TEST_CLAMPEDY0Y1",
                                     "TEST_CLAMPEDZ0","TEST_CLAMPEDZ1","TEST_CLAMPEDZ0Z1",
                                     "TEST_MANUAL",
                                     "VFPresetName","",0};
typedef enum {
  FRACTURE,
  ELASTICITY
} VFMode;
static const char *VFModeName[] = {
  "FRACTURE",
  "ELASTICITY",
  "VFModeName",
  "",
  0
};

typedef enum {
  UNILATERAL_NONE,
  UNILATERAL_SHEARONLY
} VFUnilateralType;  
static const char *VFUnilateralName[] = {
  "NONE",
  "SHEARONLY",
  "VFUnilateralName",
  "",
  0
};
    
typedef enum {
  COUPLING_NONE,
  COUPLING_GMRSTOVF,
  COUPLING_FULL
} VFCouplingType;
static const char *VFCouplingName[] = {
  "NONE",
  "GMRSTOVF",
  "FULL",
  "VFCouplingName",
  "",
  0
};

typedef enum {
  FILEFORMAT_BIN,
  FILEFORMAT_HDF5
  } VFFileFormatType;
static const char *VFFileFormatName[] = {
  "bin",
  "hdf5",
  "VFFileFormatName",
  "",
  0
};
/* 
  all fields involved in the computations
*/
typedef struct {
  Vec V;
  Vec VIrrev;
  Vec U;
  Vec BCU;
  Vec theta;
  Vec thetaRef;
  Vec pressure;
  Vec pmult;
} VFFields;

typedef struct {
  PetscReal       E,nu;       /* Young modulus and poisson ratio */
  PetscReal       lambda,mu;  /* Lame coefficients               */
  PetscReal       alpha;      /* Linear thermal expansion coef.  */
  PetscReal       Gc;         /* Fracture toughness              */
  PetscReal       beta;       /* Biot's constant                 */
  PetscReal       rho;        /* density */
} MatProp;

typedef struct {
  PetscReal        epsilon;
  PetscReal        eta;
  PetscReal        atCv;
  PetscReal        irrevtol;
  PetscReal        permVcutoff;
} VFProp;

typedef struct {
  PetscLogStage VF_IOStage;
  
  PetscLogStage VF_UAssemblyStage;
  PetscCookie   VF_MatULocalCookie;
  PetscLogEvent VF_MatULocalEvent;
  PetscCookie   VF_VecULocalCookie;
  PetscLogEvent VF_VecULocalEvent;
  
  PetscLogStage VF_USolverStage;

  PetscLogStage VF_VAssemblyStage;
  PetscCookie   VF_MatVLocalCookie;
  PetscLogEvent VF_MatVLocalEvent;
  PetscCookie   VF_VecVLocalCookie;
  PetscLogEvent VF_VecVLocalEvent;

  PetscLogStage VF_VSolverStage;

  PetscLogStage VF_EnergyStage;
  PetscCookie   VF_EnergyLocalCookie;
  PetscLogEvent VF_EnergyLocalEvent;
} VFLog;

typedef struct {
  PetscInt            ncellx,ncelly,ncellz;
  PetscInt            nlayer;
  PetscReal           *layersep;
  PetscInt            *layer;         /* dim=nz+1. gives the layer number of a cell  */
  PetscReal           BoundingBox[6]; /* Reservoir bounding box [Xmin, Xmax, Ymin, Ymax, Zmin, Zmax] */
  BC                  bcU[3];
  BC                  bcV[1];
  DA                  daVect;
  DA                  daScal;
  CartFE_Element3D    e3D;
  char                prefix[PETSC_MAX_PATH_LEN];
  Vec                 coordinates;
  PetscInt            verbose;
  VFPreset            preset;
  Mat                 KU;
  PC                  pcU;
  KSP                 kspU;
  Vec                 RHSU;
  Mat                 KV;
  PC                  pcV;
  KSP                 kspV;
  Vec                 RHSV;
  PetscReal           altmintol;
  PetscInt            altminmaxit;
  MatProp             *matprop;
  VFProp              vfprop;
  VFLog               vflog;
  PetscReal           insitumin[6];
  PetscReal           insitumax[6];
  VFMode              mode;
  PetscInt            saveint;    /* Save partial step every saveint alternate minimizations iterations */
  VFUnilateralType    unilateral;
  VFCouplingType      coupling;
  VFFileFormatType    fileformat;
  PetscViewer         energyviewer;
  PetscViewer         XDMFviewer;
  PetscReal           timevalue;
  PetscInt            timestep;
  PetscReal           ElasticEnergy;
  PetscReal           SurfaceEnergy;
  PetscReal           InsituWork;
  PetscReal           TotalEnergy;
} VFCtx;

extern VFCtx          ctx;
extern VFFields       fields;

extern PetscErrorCode OldVFInitialize(VFCtx *ctx,VFFields *fields);  
extern PetscErrorCode VFCtxGet(VFCtx *ctx);
extern PetscErrorCode VFGeometryInitialize(VFCtx *ctx,PetscReal *dx,PetscReal *dy,PetscReal *dz);
extern PetscErrorCode VFFieldsInitialize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VFBCInitialize(VFCtx *ctx);
extern PetscErrorCode VFSolversInitialize(VFCtx *ctx);
extern PetscErrorCode VFTimeStepPrepare(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VFElasticityTimeStep(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VFFractureTimeStep(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VFFinalize(VFCtx *ctx,VFFields *fields);

extern PetscErrorCode VFPropGet(VFProp *vfprop);
extern PetscErrorCode VFMatPropGet(MatProp *matprop,PetscInt n);
extern PetscErrorCode VFLayerInit(VFCtx *ctx);

extern PetscErrorCode FieldsH5Write(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FieldsBinaryWrite(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VFLogInitialize(VFLog *vflog);

extern PetscErrorCode PermUpdateTruncate(Vec V,Vec Pmult,VFProp *vfprop,VFCtx *ctx);

#endif /* VFCOMMON_H */
