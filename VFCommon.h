#ifndef VFCOMMON_H
#define VFCOMMON_H
/*
 VFCommon.h
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
static const char banner[] = "\n\nVF:\nNumerical implementation of the variational approach to fracture.\n(c) 2010-2012 Blaise Bourdin, Louisiana State University. bourdin@lsu.edu\n\n";

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
	UnitaryUnits,			/* All variables are unitary, for testing purposes		*/
	FieldUnits,				/* Flow computation in field units						*/
	MetricUnits				/* Flow computation in metric units						*/
} FlowUnit; 

static const char *FlowUnitName[] = {
	"UnitaryUnit",
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
	TEST_MANUAL
} VFPreset;
static const char *VFPresetName[] = {"SYMXY","SYMX","SYMY","NOSYM",
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
    FLOWSOLVER_TS,
	FLOWSOLVER_SNES,
	FLOWSOLVER_FEM,
	FLOWSOLVER_DARCYMIXEDFEMSTEADYSTATE,
	FLOWSOLVER_FAKE,
	FLOWSOLVER_READFROMFILES,
} VFFlowSolverType;
static const char *VFFlowSolverName[] = {
    "TS",
	"SNES",
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
	Vec FVCellndof;
	Vec FVCell;
	Vec	VolCrackOpening;
	Vec FlowBCArray;
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
	char         name[256];
	PetscReal    top[3];
	PetscReal    bottom[3];
	/*
	 PetscReal    rate;
	 BCTYPE       BCV;
	 */
} VFWell;

typedef struct {
	char          name[256];
	PetscReal     center[3];
	PetscReal     r,phi,theta;
} VFPennyCrack;

typedef struct {
	char          name[256];
	PetscReal     corners[9];
	/*
	Corner allocation:
	
	[6,7,8]
	   +---------------------+
	   |                     |
	   |                     |
	   |                     |
	   +---------------------+
	[0,1,2]               [3,4,5]
	*/   
} VFRectangularCrack;

typedef struct {
	PetscReal         perm;     /* Permeability in m^2 muliply by 1e12 */
	PetscReal         por;      /* Porosity */
	PetscReal         Pinit;    /* Initial Pressure in MPa*/
	PetscReal         Tinit;    /* Initial Temperature in C*/ 
	PetscReal         relk;     /* Relative Permeability */
	PetscReal         visc;     /* Viscosity in cp */
	PetscReal         fdens;    /* Fluid Density in specific density*/
	PetscReal         rock_comp; /* Rock compressibility in 1/MPa */
	PetscReal		      wat_comp; /* Water compressibility in 1/MPa */
	PetscReal         TCond_X;  /* Thermal Conductivity in x-direction */
	PetscReal         TCond_Y;  /* Thermal Conductivity in y-direction */
	PetscReal         TCond_Z;  /* Thermal Conductivity in z-direction */
	/*
	 change them to Vec later.
	 Instead, I would suggest keeping the structure this way and add a pointer to a resprop in the main context
	 We can read them in a file later, the difficulty is to initialize the files to something reasonable
	 */
} ResProp; 

typedef struct {
	PetscLogStage VF_IOStage;
	
	PetscLogStage VF_UAssemblyStage;
  //PetscClassId  VF_MatULocalClassId;
  //PetscLogEvent VF_MatULocalEvent;
  //PetscClassId  VF_VecULocalClassId;
  //PetscLogEvent VF_VecULocalEvent;
	
	PetscLogStage VF_USolverStage;
	
	PetscLogStage VF_VAssemblyStage;
	//PetscClassId  VF_MatVLocalClassId;
	//PetscLogEvent VF_MatVLocalEvent;
	//PetscClassId  VF_VecVLocalClassId;
	//PetscLogEvent VF_VecVLocalEvent;
	
	PetscLogStage VF_VSolverStage;
	
	PetscLogStage VF_EnergyStage;
	//PetscClassId  VF_EnergyLocalClassId;
	//PetscLogEvent VF_EnergyLocalEvent;
	
	PetscLogStage VF_PAssemblyStage;
	//PetscClassId  VF_MatPLocalClassId;
	//PetscLogEvent VF_MatPLocalEvent;
	//PetscClassId  VF_VecPLocalClassId;
	//PetscLogEvent VF_VecPLocalEvent;
	
	PetscLogStage VF_PSolverStage;
	
	PetscLogStage VF_TAssemblyStage;
	//PetscClassId  VF_MatTLocalClassId;
	//PetscLogEvent VF_MatTLocalEvent;
	//PetscClassId  VF_VecTLocalClassId;
	//PetscLogEvent VF_VecTLocalEvent;
	
	PetscLogStage VF_TSolverStage;
} VFLog;

typedef struct {
	PetscBool           printhelp;
	PetscInt            nlayer;
	PetscReal          *layersep;
	PetscInt           *layer;         /* dim=nz+1. gives the layer number of a cell  */
	BC                  bcU[3];
	BC                  bcV[1];
	BC                  bcP[1];
	BC                  bcT[1];
	DM                  daVect;
	DM                  daScal;
	CartFE_Element3D    e3D;
	CartFE_Element2D    e2D;
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
	PetscReal			      CrackVolume;
	/* 
	 Global variables for Mixed Darcy Flow
	 */
	Mat						      KVelP;
	PC						      pcVelP;
	KSP						      kspVelP;
	FLOWBC				      bcFlow[4];
	DM						      daFlow;
	DM						      daVFperm;
	FlowProp			      flowprop;
	Vec			    	      RHSVelP;
	FlowUnit			      units;
	FlowCases			      flowcase;
	Vec				  	      Source;

	
	/*
	 Global Variables for Heat Transfer
	 */	
	Mat                 KT;
	PC                  pcT;
	KSP                 kspT;
	Vec                 RHST;
	
	/* 
	 SNES solver for Pressure (or flow - T&P)
	 */
	SNES                snesF;
	
	PetscReal           altmintol;
	PetscInt            altminmaxit;
	MatProp            *matprop;
	ResProp             resprop;
	VFProp              vfprop;
	VFLog               vflog;
	PetscReal           insitumin[6];
	PetscReal           insitumax[6];
	PetscBool           hasInsitu;
	PetscBool           hasCrackPressure;
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
	PetscInt            numWells;
	VFWell             *well;
	PetscInt            numCracks;
	VFPennyCrack       *crack;
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
