/*
 VFCommon_private.h
 This has to be the stupidest name in the world...
 We'll have to dela with it until the project is completely reorganized

 (c) 2010-2013 Blaise Bourdin bourdin@lsu.edu
*/
#ifndef VF_Standalone_VFCommon_private_h
#define VF_Standalone_VFCommon_private_h

static const char *VFUnitName[] = {
        "UnitaryUnit",
        "FieldUnit",
        "MetricUnit",
        "VFUnitName",
        "",
        0
};

static const char *VFFlowBC_Case[] = {
	"ALLNORMALFLOWBC",
	"ALLPRESSUREBC",
	"VFFlowBC_Case",
	"",
	0
};

static const char *VFPresetName[] = {"SYMXY",
  "SYMX",
  "SYMY",
  "NOSYM",
	"TEST_MANUAL",
	"VFPresetName",
	"",
	0};

static const char *VFMechSolverName[] = {
	"FRACTURE",
	"ELASTICITY",
	"NOMECH",
	"VFMechSolverName",
	"",
	0
};

static const char *VFUnilateralName[] = {
	"NONE",
	"SHEARONLY",
	"VFUnilateralName",
	"",
	0
};

static const char *VFFlowSolverName[] = {
	"FLOWSOLVER_TS",
	"FLOWSOLVER_SNES",
	"FLOWSOLVER_FEM",
	"FLOWSOLVER_KSPMIXEDFEM",
	"FLOWSOLVER_SNESMIXEDFEM",
	"FLOWSOLVER_TSMIXEDFEM",
	"FAKE",
	"READFROMFILES",
	"FLOWSOLVER_NONE",
	"",
	0
};

static const char *VFFractureFlowSolverName[] = {
	"FRACTUREFLOWSOLVER_SNESMIXEDFEM",
	"FRACTUREFLOWSOLVER_NONE",
	"VFFractureFlowSolverName",
	"",
	0
};

static const char *VFHeatSolverName[] = {
	"HEATSOLVER_SNESFEM",
	"VFHeatSolverName",
	"",
	0
};

static const char *VFFileFormatName[] = {
	"bin",
	"hdf5",
	"VFFileFormatName",
	"",
	0
};






#endif
