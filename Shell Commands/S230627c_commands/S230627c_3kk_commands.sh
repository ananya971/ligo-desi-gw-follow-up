
WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_103631_451645  20230627 i
cd ligodesi_S230627c_103631_451645/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=159.132741406714  DEC=45.2793855713682  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_104309_523323  20230627 i
cd ligodesi_S230627c_104309_523323/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=160.788605223521  DEC=52.5566485559508  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_104912_575014  20230627 i
cd ligodesi_S230627c_104912_575014/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=162.300420730946  DEC=57.8372551359261  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_104524_515616  20230627 i
cd ligodesi_S230627c_104524_515616/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=161.351290633544  DEC=51.9378457634401  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_104229_494334  20230627 i
cd ligodesi_S230627c_104229_494334/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=160.623914804067  DEC=49.7263326582477  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_105127_574236  20230627 i
cd ligodesi_S230627c_105127_574236/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=162.862873413447  DEC=57.7100968822253  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_104303_512905  20230627 i
cd ligodesi_S230627c_104303_512905/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=160.763119607896  DEC=51.4848207404392  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_104002_453739  20230627 i
cd ligodesi_S230627c_104002_453739/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=160.009499232963  DEC=45.6275428108155  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_102937_382422  20230627 i
cd ligodesi_S230627c_102937_382422/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=157.407258495435  DEC=38.4061125626019  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits

WORKDIR=/data/wst/u/ps1ipp/arri/3kk
mkdir $WORKDIR/bias
mkdir $WORKDIR/dark
mkdir $WORKDIR/flat
cd $WORKDIR
/data/wst/u/ps1ipp/arri/3kk/script_3kk_ligodesi.sh  ligodesi_S230627c_103814_470157  20230627 i
cd ligodesi_S230627c_103814_470157/20230627_i
/data/wst/u/ps1ipp/arri/3kk/script_3kk_header_for_astrometry.sh  RA=159.558488055666  DEC=47.0325902417434  CCD=1  mepdnybo-n2023*.fits
sethead  FILTER=i mepdnybo-n*.fits
for i in mepdnybo-*.fits; do  script.astrometry.sh $i 1;  done
FIRST=`echo vme*.fits | awk '{print $1}'`
divvy  -e          -2  -3  -wcsok  -ref ${FIRST}    vme*.fits 
sumfits -n tvvmepdnybo-n2023*.fits
