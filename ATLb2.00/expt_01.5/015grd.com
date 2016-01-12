#!/bin/csh
#$ -S /bin/csh
#$ -pe mpi_4hr_ibm_dnode   16
#$ -v  GRD_TOTAL_MPI_TASKS=16
#$ -cwd
#
set echo
set timestamp
C
C --- Preamble.
C
setenv OS `uname`
switch ($OS)
case 'AIX':
    setenv TMPDIR /usr/var/tmp/${user}
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    echo 'configured for AIX only'
    exit (1)
endsw
C
if ($?JOB_NAME) then
    setenv GRD_JOB_NAME ${REQNAME}
    setenv GRD_STEP_ID  $$
endif
C
C --- Automatic Run Script.
C --- Submit with msub, or msub_csh, or msub_grdcommand.
C --- Multiple segment version, set no. of segments on foreach below.
C
cd /usr/people/${user}/hycom/ATLb2.00/expt_01.5
C
C --- E is expt, P is permenant directory, S is /tmp directory.
C
setenv E 015
setenv P $cwd
switch ($OS)
case 'AIX':
    setenv S `echo $cwd | awk '{print "/usr/var/tmp" substr($0,12,length)}'`
    breaksw
endsw
C
ls -laFq
C
C --- check the RUNNING flag.
C
if ( -e RUNNING && ! -e RUNNING_$GRD_STEP_ID) then
C
C --- MODEL IS ALREADY RUNNING - EXIT.
C
  exit
endif
touch RUNNING
touch RUNNING_$GRD_STEP_ID
C
C --- Number of segments specified by '( SEG1 SEG2 ... )' in foreach.
C --- Use '( SEG1 )' for a single segment per run.
C
foreach seg (SEG1)
echo "Starting SEGMENT " $seg
C
C --- Check the model list.
C
touch LIST
if ( -z LIST ) then
C
C --- MODEL RUN IS COMPLETE - EXIT.
C
  /bin/rm -f RUNNING
  /bin/rm -f RUNNING_$GRD_STEP_ID
  exit
endif
C
C --- Take the next run from LIST.
C
setenv Y01 `head -1 LIST | awk '{print $1}'`
C
if ( `echo $Y01 | egrep 'com' | wc -l` != 0 ) then
C
C --- Precalculated model script.
C
  setenv SCRIPT ${Y01}
  if (! -e $SCRIPT ) then
C
C --- PREDEFINED SCRIPT DOES NOT EXIST - EXIT.
C
    ls -laFq
    /bin/rm -f RUNNING
    /bin/rm -f RUNNING_$GRD_STEP_ID
    exit
  endif
else
C
C --- Generate the next model script.
C
  setenv AB `head -1 LIST | awk '{print $2}'`
  setenv SCRIPT ${E}y${Y01}${AB}.com
  /bin/rm -f ${SCRIPT}
  awk -f ${E}.awk y01=${Y01} ab=${AB} ${E}.com > ${SCRIPT}
endif
C
C --- Run the Script.
C
set script = $SCRIPT
set reqname = `echo ${GRD_JOB_NAME}`
ln -fs ${reqname}.log $script:r.log
C
mkdir -p $S
cp  ${SCRIPT} $S/${SCRIPT}
cd  $S
C
C ---------------------------------------------------------------------------
csh ${SCRIPT}
C ---------------------------------------------------------------------------
C
cd  $P
C
C --- Clean Up.
C
ls -laFq
C
mv LIST LISTOLD
head -1 LISTOLD | comm -23 LISTOLD - > LIST
diff LIST LISTOLD
/bin/rm -f LISTOLD
C
C --- end of SEGMENT foreach loop
C
end
C
C --- Final Clean Up.
C
/bin/rm -f RUNNING
/bin/rm -f RUNNING_$GRD_STEP_ID
C
touch LIST
if ( -z LIST ) then
C
C --- MODEL RUN IS COMPLETE - EXIT.
C
  exit
endif
C
C --- Submit the next job.
C
cat - > $TMPDIR/$$.awk <<'E-o-A'
BEGIN { for ( i=65;i <= 89; i++)
		c[sprintf("%c",i)] = sprintf("%c",i+1)
}
/[0-9]$/ { printf("%sA\n",$0) }
/[A-Y]$/ { printf("%s%s\n",substr($0,1,length($0)-1),c[substr($0,length($0))]) }
'E-o-A'
#
set newname = `echo ${reqname} | awk -f $TMPDIR/$$.awk`
#
/bin/rm -f $TMPDIR/$$.awk
if ( `echo $newname | wc -l` != 0 ) then
    if ($?JOBNAME) then
        setenv JOBNAME ${newname}
#       /bin/nice -7 csh ${E}grd.com >& ${newname}.log &
        csh ${E}grd.com >& ${newname}.log &
    else
        qsub -N      ${newname} \
                   -o $cwd/${newname}.log -j y ${E}grd.com
#       qsub -o ${newname}.log ${E}grd.com
#        sed -e "s/^# *@ *job_name.*=.*/# @ job_name = ${newname}/" ${E}grd.com | qsub -
        if ($status != 0) then
          uname -a
          llclass -l priority
        endif
    endif
else
C
C --- newname failed
C
    echo $newname
    echo $newname | wc -l
endif
C
C --- Exit.
C
exit
