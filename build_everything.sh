set -e
export PATH=/usr/local/bin:/usr/bin:/bin:$PATH
export NCDF_PREFIX=`nc-config --prefix`
export EXTRANCDF=`nc-config --flibs`
export MY_ARCH=gfortran
export BIN_OS=LinuxGF

cd ALL
echo "Now Cleaning..."
csh Make_clean.com &> /dev/null

echo "Setting root info..."
sed -i "s/^setenv ARCH amd64/setenv ARCH $MY_ARCH/" Make_all.src
sed -i "s[^setenv NCDF .*[setenv NCDF $NCDF_PREFIX[" Make_ncdf.src
sed -i "s[^setenv EXTRANCDF .*[setenv EXTRANCDF \"$EXTRANCDF\"[" Make_ncdf.src

echo "Building bin..."
cd bin
sed -i "s/^# setenv OS $BIN_OS/setenv OS $BIN_OS/" Make_all.com
# Cause script to exit with non-zero value on compile error.
sed -i 's[\($FC .*\)[\1\n\tif ($status) exit 1[' Make_all.com
# Fix some code so that build works.
sed -i "s[EXTERNAL IARGC[[" hycom_skill.F
# Remove some targets that do not build.
sed -i "s[\bhycom_profile_argo\b[[" Make_all.com
sed -i 's/\bhycom_profile_hybgen.*\b//' Make_all.com
csh Make_clean.com &> /dev/null
csh -x Make_all.com
csh -x Make_ncdf.com
cd ..

echo "Building topo..."
cd topo/src
sed -i "s/^setenv ARCH amd64/setenv ARCH $MY_ARCH/" Make_all.src
# Cause stdout and stderr to be sent to output instead of log files.
sed -i "s[>&! Make_${m}[[" Make_all.com
sed -i "s[>&! Make_${m}[[" Make_ncdf.com
# Cause the make scripts to return non-zero on error.
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_ncdf.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_all.com
# Fix build errors as directed by Allen.
sed -i 's/a,i\//a,i10\//' bathy_01min.f
sed -i 's/a,i\//a,i10\//' landsea_01min.f
sed -i 's/a,i\//a,i10\//' bathy_30sec.f
sed -i 's/a,i\//a,i10\//' landsea_30sec.f
csh Make_all.com
csh Make_ncdf.com
cd ../..

echo "Building archive..."
cd archive/src
sed -i "s/^setenv ARCH amd64/setenv ARCH $MY_ARCH/" Make_all.src
sed -i "s[>&! Make_${m}[[" Make_all.com
sed -i "s[>&! Make_${m}[[" Make_ncdf.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_all.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_ncdf.com
csh -x Make_all.com
csh -x Make_ncdf.com
cd ../..

echo "Building force..."
cd force/src
sed -i "s/^setenv ARCH amd64/setenv ARCH $MY_ARCH/" Make_all.src
sed -i "s[>&! Make_${m}[[" Make_all.com
sed -i "s[>&! Make_${m}[[" Make_ncdf.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_all.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_ncdf.com
csh Make_all.com
csh Make_ncdf.com
cd ../..

echo "Building meanstd..."
# As confirmed by Allen there are no netCDF tools in this directory.
cd meanstd/src
sed -i "s/^setenv ARCH amd64/setenv ARCH $MY_ARCH/" Make_all.src
sed -i "s[>&! Make_${m}[[" Make_all.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_all.com
csh Make_all.com
cd ../..

#echo "Building plot..."
#cd plot/src
#sed -i 's/^setenv ARCH amd64/setenv ARCH gfortran/' Make_all.src
#sed -i "s[>&! Make_${m}[[" Make_all.com
#sed -i "s[>&! Make_${m}[[" Make_ncdf.com
#sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_all.com
#sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_ncdf.com
#csh Make_all.com
#csh Make_ncdf.com
#cd ../..

echo "Building relax..."
cd relax/src
sed -i "s/^setenv ARCH amd64/setenv ARCH $MY_ARCH/" Make_all.src
sed -i "s[>&! Make_${m}[[" Make_all.com
sed -i "s[>&! Make_${m}[[" Make_ncdf.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_all.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_ncdf.com
sed -i 's/^sst_modas:       $(MODS) sst_modas.o              $(LIBN)/sst_modas:       $(MODS) sst_modas.o/' Makefile
csh Make_all.com
csh Make_ncdf.com
cd ../..

echo "Building sample..."
cd sample/src
sed -i "s/^setenv ARCH amd64/setenv ARCH $MY_ARCH/" Make_all.src
sed -i "s[>&! Make_${m}[[" Make_all.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_all.com
csh Make_all.com
cd ../..

echo "Building subregion..."
cd subregion/src
sed -i "s/^setenv ARCH amd64/setenv ARCH $MY_ARCH/" Make_all.src
sed -i "s[>&! Make_${m}[[" Make_all.com
sed -i "s[>&! Make_${m}[[" Make_ncdf.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_all.com
sed -i 's/echo "Make failed:" \${m}/echo "Make failed:" \${m} \n\texit 99/' Make_ncdf.com
# Remove some targets that don't build because no matching code file. According to Allen 
# these are old and should no longer be needed.
sed -i "s[\bsub_topog\b[[" Make_all.com
sed -i "s[\bsub_topog_panam\b[[" Make_all.com
sed -i "s[\bsubregion\b[[" Make_all.com
sed -i "s[\bisub_topog\b[[" Make_all.com
# Remove this link which should be removed by Make_all.com.
rm netcdf.inc
csh Make_all.com
csh Make_ncdf.com
cd ../..

echo "Building model..."
cd ../ATLb2.00
# Need to create a config file for the model for gfortran.
rm -f config/intelGF_one
cat > config/intelGF_one <<End-of-message
FC            =	gfortran
FCFFLAGS      =	-fPIC -fno-second-underscore -O2 -march=native -m64 -mcmodel=medium -fdefault-real-8 -fdefault-double-8
CC            =	gcc
CCFLAGS       =	-O -m64 -mcmodel=medium
CPP           =	cpp -P
CPPFLAGS      =	-DIA32 -DREAL8 -DENDIAN_IO -DTIMER
LD            =	\$(FC)
LDFLAGS       =	-v \$(FCFFLAGS)
EXTRALIBS     = 
SHELL         = /bin/sh
RM            = \rm -f
.c.o:
	\$(CC) \$(CPPFLAGS) \$(CCFLAGS)  -c \$*.c

.f.o:
	\$(FC)             \$(FCFFLAGS) -c \$*.f

.F.o:
	\$(FC) \$(CPPFLAGS) \$(FCFFLAGS) -c \$*.F
End-of-message
# Change the makefile to use this new config file.
cd src_2.2.18_22_one
sed -i 's[include ../config/$(ARCH)_$(TYPE)[include ../config/intelGF_one[' Makefile
# As directed by Allen, fix some code issues.
sed -i "s[iabs[abs[" mod_floats.F
# Build hycom
make
