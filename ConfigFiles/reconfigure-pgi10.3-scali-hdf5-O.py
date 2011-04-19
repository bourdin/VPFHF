#!/usr/bin/python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--with-debugging=1',
    '--with-clib-autodetect=0',
    '--with-mpi-dir=/opt/scali/',
    '--with-shared=0',
    '--with-blas-lapack-dir=/vend/pgi/10.3/linux86-64/10.3/lib',
    '--with-hdf5-lib=[/data/etc_rpe_ucr/VF/hdf5-1.8.4-pgi10.3-scali/lib/libhdf5.a,/data/etc_rpe_ucr/VF/hdf5-1.8.4-pgi10.3-scali/lib/libhdf5_hl.a,/opt/scali/lib64/libmpio.so,/usr/lib64/librt.so,/usr/lib64/libz.so]',
    'LIBS=-L/vend/pgi/10.3/linux86-64/10.3/lib -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lnspgc -lpgc -lrt -lpthread -lm -lgcc -lc -lgcc',
    'PETSC_ARCH=pgi10.3-scali-hdf5-O',
    '--with-fortranlib-autodetect=0',
    'CFLAGS=-tp k8-64 -mcmodel=medium -fast',
    'FFLAGS=-tp k8-64 -mcmodel=medium -fast',
    '--with-hdf5-include=/data/etc_rpe_ucr/VF/hdf5-1.8.4-pgi10.3-scali/include',
    '--with-vendor-compilers=pgi',
  ]
  configure.petsc_configure(configure_options)
