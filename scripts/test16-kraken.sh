#!/bin/bash      
#PBS -V                   # Inherit the submission environment
#PBS -N HS2D              # Job Name
#PBS -lnodes=4:ppn=12
#PBS -q normal            # Queue name "normal"
#PBS -l walltime=12:00:00 # Run time (hh:mm:ss) - 1.5 hours
#PBS -M bourdin@lsu.edu   # Use email notification address
#PBS -m be                # Email at Begin and End of job


export JOBID=${PBS_JOBID}
export WORKDIR=${PBS_O_WORKDIR}/${JOBID}
mkdir -p ${WORKDIR}
cd ${WORKDIR}

if [ -z "$NX" ]; then
  export NX=354
fi
if [ -z "$NY" ]; then
  export NY=2
fi
if [ -z "$NZ" ]; then
  export NZ=354
fi
if [ -z "$LX" ]; then
  export LX=8
fi
if [ -z "$LY" ]; then
  export LY=.01
fi
if [ -z "$LZ" ]; then
  export LZ=8.
fi
if [ -z "$C0_X" ]; then
  export C0_X=4.
fi
if [ -z "$C0_Y" ]; then
  export C0_Y=.005
fi
if [ -z "$C0_Z" ]; then
  export C0_Z=4.
fi
if [ -z "$C0_R" ]; then
  export C0_R=.2
fi
if [ -z "$C0_THETA" ]; then
  export C0_THETA=0.
fi
if [ -z "$C0_PHI" ]; then
  export C0_PHI=0.
fi
if [ -z "$C0_THICKNESS" ]; then
  export C0_THICKNESS=.045
fi
if [ -z "$MINVOL" ]; then
  export MINVOL=0.
fi
if [ -z "$MAXVOL" ]; then
  export MAXVOL=0.025
fi
if [ -z "$NUMSTEP" ]; then
  export NUMSTEP=100
fi
if [ -z "$EPSILON" ]; then
  export EPSILON=.045
fi
if [ -z "$GC" ]; then
  export GC=.8
fi
if [ -z "$ETA" ]; then
  export ETA=1e-8
fi
if [ -z "$MODE" ]; then
  export MODE=1
fi
if [ -z "$INSITUMIN" ]; then
  export INSITUMIN='0,0,0,0,0,0'
fi
if [ -z "$INSITUMAX" ]; then
  export INSITUMAX='0,0,0,0,0,0'
fi

echo  JOBID          ${JOBID}         >> 00_INFO.txt
echo  NX             ${NX}            >> 00_INFO.txt 
echo  NY             ${NY}            >> 00_INFO.txt 
echo  NZ             ${NZ}            >> 00_INFO.txt 
echo  LX             ${LX}            >> 00_INFO.txt 
echo  LY             ${LY}            >> 00_INFO.txt 
echo  LZ             ${LZ}            >> 00_INFO.txt 
echo  GC             ${GC}            >> 00_INFO.txt
echo  C0_R           ${C0_R}          >> 00_INFO.txt
echo  C0_THETA       ${C0_THETA}      >> 00_INFO.txt
echo  C0_PHI         ${C0_PHI}        >> 00_INFO.txt
echo  C0_THICKNESS   ${C0_THICKNESS}  >> 00_INFO.txt
echo  INSITUMIN      ${INSITUMIN}     >> 00_INFO.txt
echo  INSITUMAX      ${INSITUMAX}     >> 00_INFO.txt
echo  EPSILON        ${EPSILON}       >> 00_INFO.txt
echo  ETA            ${ETA}           >> 00_INFO.txt
echo  MODE           ${MODE}          >> 00_INFO.txt
echo  EXTRAOPTS      ${EXTRAOPTS}     >> 00_INFO.txt


mpiexec ${VFDIR}/ValidationTests/test16 -p ${JOBID} -n ${NX},${NY},${NZ} -l ${LX},${LY},${LZ} -epsilon ${EPSILON} -eta ${ETA} -gc ${GC} -nu 0. \
      -nc 1 -c0_center ${C0_X},${C0_Y},${C0_Z} -c0_R ${C0_R} -c0_theta ${C0_THETA} -c0_phi ${C0_PHI} -c0_thickness ${C0_THICKNESS} \
      -mode ${MODE} -minvol ${MINVOL} -maxvol ${MAXVOL} -maxtimestep ${NUMSTEP} -insitumin ${INSITUMIN} -insitumax ${INSITUMAX} ${EXTRAOPTS} \
      -U_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 -U_mg_levels_ksp_type chebyshev -U_mg_levels_pc_type sor \
      -U_pc_gamg_agg_nsmooths 1 -U_pc_gamg_threshold 0 -U_pc_gamg_type agg -U_pc_type gamg \
      -V_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 -V_mg_levels_ksp_type chebyshev -V_mg_levels_pc_type sor \
      -V_pc_gamg_agg_nsmooths 1 -V_pc_gamg_threshold 0 -V_pc_gamg_type agg -V_pc_type gamg \
      -altmintol 1e-4 -altminmaxit 100


