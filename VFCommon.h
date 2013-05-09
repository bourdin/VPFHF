#ifndef VFCOMMON_H
#define VFCOMMON_H
/*
 VFCommon.h
 (c) 2010-2013 Blaise Bourdin bourdin@lsu.edu
 */
static const char banner[] = "\n\nVF:\nNumerical implementation of the variational approach to fracture.\n(c) 2010-2012 Blaise Bourdin, Louisiana State University. bourdin@lsu.edu\n\n";

typedef struct {
	PetscReal   mu;             /* Fluid viscosity              */
	PetscReal   rho;            /* Fluid density                */
	PetscReal   por;
	PetscReal   cf;             /* Fluid compressibility        */
	PetscReal   M_inv;
	PetscReal   beta;           /* Conversion constant          */
	PetscReal   gamma;          /*Conversion parameter          */
	PetscReal   alpha;          /*Conversion parameter          */
	PetscReal   theta;          /*Time parameter                */
	PetscReal   timestepsize;   /*Time step size                */
	PetscReal   g[3];
	PetscReal   Cp;			  /* Specific heat capacity */
	PetscReal   Kw;         /* Modulus of liquid               */
} VFFlowProp;

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
	Vec VolCrackOpening;
	Vec VolLeakOffRate;
	Vec FlowBCArray;
	Vec PresBCArray;
	Vec fracpressure;
	Vec fracvelocity;
	Vec fracVelnPress;
} VFFields;

/*
All static const are moved into the files where they are used in order to avoid annoying warnings

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
*/
typedef struct {
	PetscReal       E,nu;       /* Young modulus and poisson ratio */
	PetscReal       lambda,mu;  /* Lame coefficients               */
	PetscReal       alpha;      /* Linear thermal expansion coef.  */
	PetscReal       Gc;         /* Fracture toughness              */
	PetscReal       beta;       /* Biot's constant                 */
	PetscReal       rho;        /* density                         */
	Vec             phi;        /* porosity                        */
	Vec             Ks;         /* Bulk modulus of rock            */
	Vec             Kw;         /* Modulus of liquid               */
	Vec             VecGc;
	PetscReal       Cp;		  /* Specific heat capacity */
} VFMatProp; //change them to Vec later

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


/*
  enums. The matching field names used in PetscOptions are in VFCommon_private.h
*/
/* 
  This needs to be moved to VFWell.h once the headers are reorganized
*/
typedef enum {
	PRESSURE,
	RATE
} WellConstraint;

typedef enum {
	INJECTOR,
	PRODUCER
} WellType;

typedef enum {
	UnitaryUnits,     /* All variables are unitary, for testing purposes    */
	FieldUnits,       /* Flow computation in field units            */
	MetricUnits       /* Flow computation in metric units           */
} VFUnit;

typedef enum {
	ALLNORMALFLOWBC,
	ALLPRESSUREBC
} VFFlowCases;

typedef enum {
	SYMXY,
	SYMX,
	SYMY,
	NOSYM,
	TEST_MANUAL
} VFPreset;

typedef enum {
	FRACTURE,
	ELASTICITY,
	NOMECH
} VFMechSolverType;

typedef enum {
	UNILATERAL_NONE,
	UNILATERAL_SHEARONLY
} VFUnilateralType;

typedef enum {
	FLOWSOLVER_TS,
	FLOWSOLVER_SNES,
	FLOWSOLVER_FEM,
	FLOWSOLVER_KSPMIXEDFEM,
	FLOWSOLVER_SNESMIXEDFEM,
	FLOWSOLVER_TSMIXEDFEM,
	FLOWSOLVER_FAKE,
	FLOWSOLVER_READFROMFILES,
	FLOWSOLVER_NONE,
} VFFlowSolverType;

typedef enum {
	FRACTUREFLOWSOLVER_SNESMIXEDFEM,
	FRACTUREFLOWSOLVER_NONE,
} VFFractureFlowSolverType;

typedef enum {
	HEATSOLVER_SNESFEM,
} VFHeatSolverType;

typedef enum {
	FILEFORMAT_BIN,
	FILEFORMAT_HDF5
} VFFileFormatType;

typedef struct {
	char         name[256];
	PetscReal    top[3];
	PetscReal    bottom[3];
	PetscReal    coords[3];
	PetscReal		Qw;
	PetscReal		Pw;
	PetscReal		rw;
	WellConstraint  condition;
	WellType		type;
	/*
	 PetscReal    rate;
	 BCTYPE       BCV;
	 */
} VFWell;





typedef struct {
	char          name[256];
	PetscReal     center[3];
	PetscReal     r,phi,theta,thickness;
} VFPennyCrack;

typedef struct {
	char          name[256];
	PetscReal     corners[9];
	PetscReal     thickness;
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
	PetscReal         wat_comp; /* Water compressibility in 1/MPa */
	PetscReal         TCond_X;  /* Thermal Conductivity in x-direction */
	PetscReal         TCond_Y;  /* Thermal Conductivity in y-direction */
	PetscReal         TCond_Z;  /* Thermal Conductivity in z-direction */
	/*
	 change them to Vec later.
	 Instead, I would suggest keeping the structure this way and add a pointer to a resprop in the main context
	 We can read them in a file later, the difficulty is to initialize the files to something reasonable
	 */
} VFResProp;

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
	SNES                snesV;
	SNES                snesU;
	Vec                 UResidual;
	Mat                 JacU;
	Mat                 KU;
	Vec                 RHSU;
	/*
	 Global variables for regular FEM Flow
	 */
	Mat                 KP; //stifness
	Mat                 KPlhs; // mass matrix
	Mat                 JacP;
	TS                  tsP;
	PC                  pcP;
	KSP                 kspP;
	SNES                snesP;
	Vec                 RHSP;
	Vec                 RHSPpre;
	Vec                 PrePressure;
	Vec                 PFunct;
	Vec                 PresBC;
	PetscReal           CrackVolume;
	PetscReal           LeakOffRate;
	/*
	 Global variables for Mixed Darcy Flow
	 */
	Mat                 KVelP;
	PC                  pcVelP;
	KSP                 kspVelP;
	DM                  daFlow;
	DM                  daVFperm;
	VFFlowProp          flowprop;
	Vec                 RHSVelP;
	Vec                 RHSVelPpre;
	VFUnit              units;
	VFFlowCases         flowcase;
	Vec                 Source;
	DM                  daScalCell;
	Mat                 KVelPlhs;
	Mat                 JacVelP;
	TS                  tsVelP;
	SNES                snesVelP;
	Vec                 FlowFunct;
	Vec                 PreFlowFields;
	Vec                 Perm;
	Vec                 FlowBC;
	BC                  bcQ[3];
	PetscBool           hasFlowWells;
	PetscBool           hasFluidSources;
	Vec                 VelBCArray;
	Vec				          PresBCArray;
	
	/*
	 Global Variables for Heat Transfer
	 */
	Mat                 KT;
	Mat                 KTlhs;
	
	PC                  pcT;
	KSP                 kspT;
	Mat                 JacT;
	Vec                 PreHeatFields;
	Vec                 RHST;
	Vec                 HeatFunct;
	Vec                 RHSTpre;
//	BC                  bcq[3];
	SNES                snesT;
	Vec                 HeatBC;
	VFUnit              Hunits;
	/*
	 SNES solver for Pressure (or flow - T&P)
	 */
	SNES                snesF;
	
	PetscReal           altmintol;
	PetscInt            altminmaxit;
	VFMatProp          *matprop;
	VFResProp           resprop;
	VFProp              vfprop;
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
	VFHeatSolverType    heatsolver;
	VFMechSolverType    mechsolver;
	VFFileFormatType    fileformat;
	PetscViewer         energyviewer;
	PetscViewer         XDMFviewer;
	PetscReal           timevalue;
	PetscReal           current_time;
	PetscReal           dt;
	PetscInt            timestep;
	PetscReal           maxtimevalue;
	PetscInt            maxtimestep;
	PetscInt            maxiterations;
	PetscReal           ElasticEnergy;
	PetscReal           SurfaceEnergy;
	PetscReal           InsituWork;
	PetscReal           PressureWork;
	PetscReal           TotalEnergy;
	PetscInt            numWells;
	VFWell             *well;
	PetscInt            numPennyCracks;
	VFPennyCrack       *pennycrack;
	PetscInt            numRectangularCracks;
	VFRectangularCrack *rectangularcrack;
	Vec                 HeatSource;
	Vec					        HeatFluxBCArray;
	PetscBool           hasHeatSources;
	Vec                 prevT;
	Vec                 Cond;
	BC                  bcQT[1];        /*heat flux in heat equation*/
	Vec					        TBCArray;
	Mat                 KFracVelP;
	Mat                 JacFracVelP;
	SNES                snesFracVelP;
	Vec                 FracResidual;
	Vec                 RHSFracVelP;
	VFFractureFlowSolverType    fractureflowsolver;
	Vec                 RHSFracVelPpre;
	Mat                 KFracVelPlhs;
	Vec                 PreFracFlowFields;
	VFFields            *fields;
  Vec                 FracVelBCArray;
  BC                  bcFracQ[3];
  Vec                 V;
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
extern PetscErrorCode VFMatPropGet(VFMatProp *matprop,PetscInt n);
extern PetscErrorCode VFLayerInit(VFCtx *ctx);
extern PetscErrorCode VFResPropGet(VFResProp *resprop);

extern PetscErrorCode FieldsH5Write(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FieldsBinaryWrite(VFCtx *ctx,VFFields *fields);

extern PetscErrorCode PermUpdate(Vec V,Vec Pmult,VFProp *vfprop,VFCtx *ctx);
extern PetscErrorCode VFMatPropFieldsInitialize(VFCtx *ctx, VFMatProp *matprop);

#endif /* VFCOMMON_H */
