#!/bin/sh

if [ -z "${VFDIR}" ]
then
	echo "ERROR: VFDIR not set"
	exit
fi
if [ ! -d $VFDIR ] 
then
	echo "ERROR: $VFDIR directory does not exist"
	exit
fi

if [ -z "${GMRSDIR}" ]
then
	echo "ERROR GMRSDIR not set"
	exit
fi
if [ ! -d $GMRSDIR ]
then
  echo "ERROR: $GMRSDIR directory does not exist"
	exit
fi

for srcdir in make sample workscali
do 
  cd $GMRSDIR/$srcdir
  echo "Creating symlinks in $GMRSDIR/$srcdir" 
  ln -sf $VFDIR/GMRS_Files/$srcdir/* .
done

