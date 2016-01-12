#!/bin/csh
#
#@ job_name         = XXXrll
#@ output           = $(job_name).log
#@ error            = $(job_name).log
#@ restart          = yes
#@ job_type         = parallel
#@ network.MPI      = csss,not_shared,US
#@ environment      = MP_EUILIB=us     
#@ node             = 1
#@ total_tasks      = 8
#@ node_usage       = not_shared
#@ resources        = ConsumableCpus(1) ConsumableMemory(700mb)
#@ wall_clock_limit = 2:00:00
#@ account_no       = NRLSS018
#@ class            = standard
#@ queue
#
set echo
set timestamp
C
C --- Preamble.
C
setenv OS `uname`
switch ($OS)
case 'AIX':
    hostname
    setenv TMPDIR /scr/${user}
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    echo 'configured for AIX only'
    exit (1)
endsw
C
if ($?JOBNAME) then
    setenv LOADL_JOB_NAME ${JOBNAME}
    setenv LOADL_STEP_ID  $$
endif
C
C --- Automatic Run Script.
C --- Submit with msub, or msub_csh, or msub_ll command.
C --- Multiple segment version, set no. of segments on foreach below.
C
cd ~/hycom/ATLb2.00/expt_02.0
C
C --- E is expt, P is permenant directory, S is /tmp directory.
C
setenv E 020
setenv P $cwd
switch ($OS)
case 'AIX':
#                           substitute /scr for /u/home
    setenv S `echo $cwd | awk '{print "/scr" substr($0,8,length)}'`
    breaksw
endsw
C
ls -laFq
C
C --- check the RUNNING flag.
C
if ( -e RUNNING && ! -e RUNNING_$LOADL_STEP_ID) then
C
C --- MODEL IS ALREADY RUNNING - EXIT.
C
  exit
endif
touch RUNNING
touch RUNNING_$LOADL_STEP_ID
C
C --- Number of segments specified by '( SEG1 SEG2 ... )' in foreach.
C --- Use '( SEG1 )' for a single segment per run.
C
#foreach seg ( SEG1 SEG2 SEG3 SEG4 SEG5 SEG6 SEG7 SEG8 SEG9 SEG10 SEG11 SEG12 SEG13 SEG14 SEG15 SEG16 SEG17 SEG18 SEG19 SEG20 SEG21 SEG22 SEG23 SEG24 )
foreach seg (SEG1 SEG2 SEG3 SEG4 SEG5)
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
  /bin/rm -f RUNNING_$LOADL_STEP_ID
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
    /bin/rm -f RUNNING_$LOADL_STEP_ID
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
set reqname = ${LOADL_JOB_NAME}
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
/bin/rm -f RUNNING_$LOADL_STEP_ID
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
/[0-9]$/  { printf("%s%s\n",$0,"A")}
/[A-Y]$/  { printf("%s%s\n",substr($0,1,length($0)-1),c[substr($0,length($0))])}
/[1-9]Z$/ { printf("%s%s%s\n",substr($0,1,length($0)-2),substr($0,length($0)-1,1)+1,"A")}
/0Z$/     { next }
'E-o-A'
#
set newname = `echo ${reqname} | awk -f $TMPDIR/$$.awk`
#
/bin/rm -f $TMPDIR/$$.awk
if ( `echo $newname | wc -l` != 0 ) then
    if ($?JOBNAME) then
        setenv JOBNAME ${newname}
#       /bin/nice -7 csh ${E}rll.com >& ${newname}.log &
        csh ${E}rll.com >& ${newname}.log &
    else
        sed -e "s/^# *@ *job_name.*=.*/# @ job_name = ${newname}/" ${E}rll.com | llsubmit -
        if ($status != 0) then
          uname -a
          llclass -l batch
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
