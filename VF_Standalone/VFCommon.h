#ifndef VFCOMMON_H
#define VFCOMMON_H
/*
  VFCommon.h
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
static const char banner[] = "\n\nVF:\nNumerical implementation of the variational approach to fracture.\n(c) 2010-2011 Blaise Bourdin, Louisiana State University. bourdin@lsu.edu\n\n";

typedef enum {
	VELOCITY,
	PRESSURE,
	NOBC
} FlowBCTYPE;

typedef struct {
	FlowBCTYPE     face[6];
	FlowBCTYPE     edge[12];
	FlowBCTYPE     vertex[8];   
} FLOWBC;

static const char *FLOWBCTYPE_NAME[] = {
"NORMALVELOCITY",
"PRESSURE",
"FLOWBCTYPE_NAME",
"",
0
};

typedef struct {
	PetscReal       mu;			/* Fluid viscosity						*/
	PetscReal       rho;        /* Fluid density						*/
	PetscReal		por;
	PetscReal       cf;         /* Fluid compressibility				*/
	PetscReal		beta;		/* Conversion constant					*/
	PetscReal		gamma;		/*Conversion parameter					*/
	PetscReal		alpha;		/*Conversion parameter					*/
	PetscReal		g[3];
} FlowProp; 

typedef enum {
	FieldUnits,			/* Flow computation in field units						*/
	MetricUnits        /* Flow computation in metric units						*/
} FlowUnit; 

static const char *FlowUnitName[] = {
"FieldUnit",
"MetricUnit",
"FlowUnitName",
"",
0
};

typedef enum { 
	ALLNORMALFLOWBC,
	ALLPRESSUREBC
} FlowCases;

static const char *FlowBC_Case[] = {
"ALLNORMALFLOWBC",
"ALLPRESSUREBC",
"FlowBC_Case",
"",
0
};

typedef enum { 
  SYMXY,
  SYMX,
  SYMY,
  NOSYM,
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
static const char *VFPresetName[] = {"SYMXY","SYMX","SYMY","NOSYM",
                                     "TEST_CLAMPEDX0","TEST_CLAMPEDX1","TEST_CLAMPEDX0X1",
                                     "TEST_CLAMPEDY0","TEST_CLAMPEDY1","TEST_CLAMPEDY0Y1",
                                     "TEST_CLAMPEDZ0","TEST_CLAMPEDZ1","TEST_CLAMPEDZ0Z1",
                                     "TEST_MANUAL",
                                     "VFPresetName","",0};
typedef enum {
  FRACTURE,
  ELASTICITY,
  NOMECH
} VFMechSolverType;
static const char *VFMechSolverName[] = {
  "FRACTURE",
  "ELASTICITY",
  "NOMECH",
  "VFMechSolverName",
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
  FLOWSOLVER_FEM,
  FLOWSOLVER_DARCYMIXEDFEMSTEADYSTATE,
  FLOWSOLVER_FAKE,
  FLOWSOLVER_READFROMFILES,
  } VFFlowSolverType;
static const char *VFFlowSolverName[] = {
  "FEM",
  "MixedFEM",
  "FAKE",
  "READFROMFILES",
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
  int numfields;
  Vec V;
  Vec VIrrev;
  Vec U;
  Vec BCU;
  Vec theta;
  Vec thetaRef;
  Vec pressure;
  Vec pressureRef;
  Vec pmult;
	Vec VelnPress;
	Vec vfperm;
	Vec velocity;
} VFFields;

static const char *VFFieldNames[] = {
  "V",
  "VIrrev",
  "Displacement",
  "BCU",
  "theta",
  "thetaref",
  "pressure",
  "pressureRef",
  "pmult",
  "",
  0
};

typedef struct {
  PetscReal       E,nu;       /* Young modulus and poisson ratio */
  PetscReal       lambda,mu;  /* Lame coefficients               */
  PetscReal       alpha;      /* Linear thermal expansion coef.  */
  PetscReal       Gc;         /* Fracture toughness              */
  PetscReal       beta;       /* Biot's constant                 */
  PetscReal       rho;        /* density                         */
} MatProp; //change them to Vec later

typedef struct {
  PetscReal        epsilon;
  PetscReal        eta;
  PetscReal        atCv;
  PetscReal        irrevtol;
  PetscReal        permmax;
  /*
    permmax should be moved to resprop
  */
} VFProp;

typedef struct {
  PetscReal         perm;  /* Permeability in m^2 muliply by 1e12 */
  PetscReal         por;   /* Porosity */
  PetscReal         Pinit; /* Initial Pressure in MPa*/
  PetscReal         Tinit; /* Initial Temperature in C*/ 
  PetscReal         relk;  /* Relative Permeability */
  PetscReal         visc;  /* Viscosity in cp */
  PetscReal         fdens; /* Fluid Density in specific density*/
  PetscReal		    cf;	     /* Rock compressibility in field unit*/
  PetscReal         TCond_X; /* Thermal Conductivity in x-direction */
  PetscReal         TCond_Y; /* Thermal Conductivity in y-direction */
  PetscReal         TCond_Z; /* Thermal COnductivity in z-direction */
  /*
    change them to Vec later.
    Instead, I would suggest keeping the structure this way and add a pointer to a resprop in the main context
    We can read them in a file later, the difficulty is to initialize the files to something reasonable
  */
} ResProp; 

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

  PetscLogStage VF_PAssemblyStage;
  PetscCookie   VF_MatPLocalCookie;
  PetscLogEvent VF_MatPLocalEvent;
  PetscCookie   VF_VecPLocalCookie;
  PetscLogEvent VF_VecPLocalEvent;

  PetscLogStage VF_PSolverStage;
  
  PetscLogStage VF_TAssemblyStage;
  PetscCookie   VF_MatTLocalCookie;
  PetscLogEvent VF_MatTLocalEvent;
  PetscCookie   VF_VecTLocalCookie;
  PetscLogEvent VF_VecTLocalEvent;
  
  PetscLogStage VF_TSolverStage;
} VFLog;

typedef struct {
  PetscInt            nlayer;
  PetscReal           *layersep;
  PetscInt            *layer;         /* dim=nz+1. gives the layer number of a cell  */
  BC                  bcU[3];
  BC                  bcV[1];
  BC                  bcP[1];
  BC                  bcT[1];
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
  Mat                 KP;
  PC                  pcP;
  KSP                 kspP;
  Vec                 RHSP;

  /* 
    Global variables for Mixed Darcy Flow
  */
	Mat					        KVelP;
	PC					        pcVelP;
	KSP					        kspVelP;
	FLOWBC				      bcFlow[4];
	DA					        daFlow;
	DA					        daVFperm;
	FlowProp			      flowprop;
	Vec					        RHSVelP;
	FlowUnit			      units;
	FlowCases			      flowcase;
	PetscReal			      flowrate;
	

  /*
    Global Variables for Heat Transfer
  */	
  Mat                 KT;
  PC                  pcT;
  KSP                 kspT;
  Vec                 RHST;

  PetscReal           altmintol;
  PetscInt            altminmaxit;
  MatProp             *matprop;
  ResProp             resprop;
  VFProp              vfprop;
  VFLog               vflog;
  PetscReal           insitumin[6];
  PetscReal           insitumax[6];
  PetscTruth          hasInsitu;
  PetscTruth          hasCrackPressure;
  PetscReal           BCpres[6];
  PetscReal           BCtheta[6];
  PetscInt            SrcLoc[3];
  PetscReal           SrcRate;
  VFUnilateralType    unilateral;
  VFFlowSolverType    flowsolver;
  VFMechSolverType    mechsolver;
  VFFileFormatType    fileformat;
  PetscViewer         energyviewer;
  PetscViewer         XDMFviewer;
  PetscReal           timevalue;
  PetscInt            timestep;
  PetscReal           maxtimevalue;
  PetscInt            maxtimestep;
  PetscReal           ElasticEnergy;
  PetscReal           SurfaceEnergy;
  PetscReal           InsituWork;
  PetscReal           PressureWork;
  PetscReal           TotalEnergy;
} VFCtx;

extern PetscErrorCode VFCtxGet(VFCtx *ctx);
extern PetscErrorCode VFInitialize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VFGeometryInitialize(VFCtx *ctx);
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
extern PetscErrorCode VFResPropGet(ResProp *resprop);

extern PetscErrorCode FieldsH5Write(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FieldsBinaryWrite(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VFLogInitialize(VFLog *vflog);

extern PetscErrorCode PermUpdate(Vec V,Vec Pmult,VFProp *vfprop,VFCtx *ctx);

#endif /* VFCOMMON_H */
