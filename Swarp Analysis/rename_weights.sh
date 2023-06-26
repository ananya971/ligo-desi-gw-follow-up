#!/bin/bash

for i in V*fits
do
   error="e.$i"
   echo ln $error `basename $i .fits`.error.fits
   ln $error `basename $i .fits`.error.fits
done
