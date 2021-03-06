#! /bin/csh -x
#BSUB -n 4
#BSUB -W 1:00
#BSUB -q primary
#BSUB -N
#
set echo
C
C --- Preamble.
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
#   assumes /usr/5bin is before /bin and /usr/bin in PATH.
    setenv TMPDIR /tmp
    breaksw
case 'OSF1':
#   setenv TMPDIR /tmp
    breaksw
case 'IRIX64':
    setenv TMPDIR /tmp
    breaksw
case 'unicos':
    setenv ACCT "NO2031"
#   newacct $ACCT
#   source ~/.login
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
if ($?JOBNAME) then
    setenv LSB_JOBNAME ${JOBNAME}
    setenv LSB_JOBID   $$
endif
C
C --- Automatic Run Script.
C --- Submit with msub, or msub_csh, or msub_lsf command.
C --- Multiple segment version, set no. of segments on foreach below.
C
cd ~/hycom/ATLb2.00/expt_01.5
C
C --- E is expt, P is permenant directory, S is /tmp directory.
C
setenv E 015
setenv P $cwd
switch ($OS)
case 'SunOS':
#                           substitute /scr for /u/home
    setenv S `echo $cwd | awk '{print "/scr" substr($0,8,length)}'`
#                           substitute /net/ajax/scrf for /home/ajax
#   setenv S `echo $cwd | awk '{print "/net/ajax/scrf" substr($0,11,length)}'`
    breaksw
case 'OSF1':
#                           substitute /work for /u
    setenv S `echo $cwd | awk '{print "/work" substr($0,3,length)}'`
#                           substitute /tmp for /???????
#   setenv S `echo $cwd | awk '{print "/tmp" substr($0,9,length)}'`
    breaksw
case 'IRIX64':
#                           substitute /scr for /u/home
    setenv S `echo $cwd | awk '{print "/scr" substr($0,8,length)}'`
    breaksw
case 'unicos':
#                           substitute /tmp for /u/b
    setenv S `echo $cwd | awk '{print "/tmp" substr($0,5,length)}'`
    breaksw
endsw
C
ls -laFq
C
C --- check the RUNNING flag.
C
if ( -e RUNNING && ! -e RUNNING_$LSB_JOBID) then
C
C --- MODEL IS ALREADY RUNNING - EXIT.
C
  exit
endif
touch RUNNING
touch RUNNING_$LSB_JOBID
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
  /bin/rm -f RUNNING_$LSB_JOBID
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
    /bin/rm -f RUNNING_$LSB_JOBID
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
set reqname = ${LSB_JOBNAME}
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
/bin/rm -f RUNNING_$LSB_JOBID
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
#       /bin/nice -7 csh ${E}lsf.com >& ${newname}.log &
        csh ${E}lsf.com >& ${newname}.log &
    else
        kinit -R
        bsub -J ${newname} \
             -o ${newname}.log < ${E}lsf.com
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
