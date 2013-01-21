#!/bin/bash      
#PBS -V
#PBS -N test16
#PBS -lnodes=4:ppn=12
#PBS -q normal
#PBS -l walltime=4:00:00


export JOBID=${PBS_JOBID}
export WORKDIR=${PBS_O_WORKDIR}/${JOBID}
mkdir -p ${WORKDIR}
cd ${WORKDIR}

touch 00_INFO.txt
echo  JOBID          ${JOBID} >> 00_INFO.txt

if [ -z "$NX" ]; then
  export NX=354
fi
echo  NX             ${NX} >> 00_INFO.txt 

if [ -z "$NY" ]; then
  export NY=2
fi
echo  NY             ${NY} >> 00_INFO.txt 

if [ -z "$NZ" ]; then
  export NZ=354
fi
echo  NZ             ${NZ} >> 00_INFO.txt 

if [ -z "$LX" ]; then
  export LX=8
fi
echo  LX             ${LX} >> 00_INFO.txt 

if [ -z "$LY" ]; then
  export LY=.01
fi
echo  LY             ${LY} >> 00_INFO.txt 

if [ -z "$LZ" ]; then
  export LZ=8.
fi
echo  LZ             ${LZ} >> 00_INFO.txt 

if [ -z "$C0_CENTER" ]; then
  export C0_CENTER='4.,.005,4.'
fi
echo  C0_CENTER           ${C0_CENTER} >> 00_INFO.txt

if [ -z "$C0_R" ]; then
  export C0_R=.2
fi
echo  C0_R           ${C0_R} >> 00_INFO.txt

if [ -z "$C0_THETA" ]; then
  export C0_THETA=0.
fi
echo  C0_THETA       ${C0_THETA} >> 00_INFO.txt

if [ -z "$C0_PHI" ]; then
  export C0_PHI=0.
fi
echo  C0_PHI         ${C0_PHI} >> 00_INFO.txt

if [ -z "$C0_THICKNESS" ]; then
  export C0_THICKNESS=.045
fi
echo  C0_THICKNESS   ${C0_THICKNESS} >> 00_INFO.txt

if [ -z "$MINVOL" ]; then
  export MINVOL=0.
fi
echo  MINVOL         ${MINVOL} >> 00_INFO.txt

if [ -z "$MAXVOL" ]; then
  export MAXVOL=0.025
fi
echo  MAXVOL         ${MAXVOL} >> 00_INFO.txt

if [ -z "$NUMSTEP" ]; then
  export NUMSTEP=100
fi
echo  NUMSTEP        ${NUMSTEP} >> 00_INFO.txt

if [ -z "$EPSILON" ]; then
  export EPSILON=.045
fi
echo  EPSILON        ${EPSILON} >> 00_INFO.txt

if [ -z "$GC" ]; then
  export GC=.8
fi
echo  GC             ${GC} >> 00_INFO.txt

if [ -z "$ETA" ]; then
  export ETA=1e-8
fi
echo  ETA            ${ETA} >> 00_INFO.txt

if [ -z "$MODE" ]; then
  export MODE=1
fi
echo  MODE           ${MODE} >> 00_INFO.txt

if [ -z "$INSITUMIN" ]; then
  export INSITUMIN='0,0,0,0,0,0'
fi
echo  INSITUMIN      ${INSITUMIN} >> 00_INFO.txt

if [ -z "$INSITUMAX" ]; then
  export INSITUMAX='0,0,0,0,0,0'
fi
echo  INSITUMAX      ${INSITUMAX} >> 00_INFO.txt

if [ -z "$KSPUOPTS" ]; then
  export KSPUOPTS='-U_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 -U_mg_levels_ksp_type chebyshev -U_mg_levels_pc_type sor -U_pc_gamg_agg_nsmooths 1 -U_pc_gamg_threshold 0 -U_pc_gamg_type agg -U_ksp_type cg -U_pc_type gamg'
fi
echo  KSPUOPTS       ${KSPUOPTS} >> 00_INFO.txt

if [ -z "$KSPVOPTS" ]; then
  export KSPVOPTS='-V_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 -V_mg_levels_ksp_type chebyshev -V_mg_levels_pc_type sor -V_pc_gamg_agg_nsmooths 1 -V_pc_gamg_threshold 0 -V_pc_gamg_type agg -V_ksp_type cg -V_pc_type gamg'
fi
echo  KSPVOPTS       ${KSPVOPTS} >> 00_INFO.txt

echo  EXTRAOPTS      ${EXTRAOPTS} >> 00_INFO.txt


mpiexec ${VFDIR}/ValidationTests/test16 -p ${JOBID} -n ${NX},${NY},${NZ} -l ${LX},${LY},${LZ} -epsilon ${EPSILON} -eta ${ETA} -gc ${GC} -nu 0. \
      -nc 1 -c0_center ${C0_CENTER} -c0_R ${C0_R} -c0_theta ${C0_THETA} -c0_phi ${C0_PHI} -c0_thickness ${C0_THICKNESS} \
      -mode ${MODE} -minvol ${MINVOL} -maxvol ${MAXVOL} -maxtimestep ${NUMSTEP} -insitumin ${INSITUMIN} -insitumax ${INSITUMAX} \
      ${KSPUOPTS} ${KSPVOPTS} ${EXTRAOPTS}

