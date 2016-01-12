#
#@$-lT 36:00:00
#@$-lt 35:55:00
#@$-lm 17Mw
#@$-lM 17Mw
#@$-lf 5Gb
#@$-s  /bin/csh
#@$-q  mtask
#@$
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
case 'Linux':
    setenv TMPDIR /tmp
    breaksw
case 'OSF1':
    setenv TMPDIR /tmp
    /bin/rm /tmp/$user
    mkdir ~/scratch
    /bin/ln -s ~/scratch /tmp/$user
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
    setenv QSUB_REQNAME R${JOBNAME}
    setenv QSUB_REQID   $$
endif
C
C --- Automatic Run Script.
C --- Submit with msub, or msub_csh, or msub_nqs command.
C --- Multiple segment version, set no. of segments on foreach below.
C
C --- SR is the scratch disk root.
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
      setenv SR /net/hermes/scrb
    else if (-e /scr) then
      setenv SR /scr
    else
      setenv SR $TMPDIR
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
      setenv SR /external/fast
    else if (-e /scr) then
      setenv SR /scr
    else
      setenv SR $TMPDIR
    endif
    breaksw
case 'OSF1':
    if (-e /scr) then
      setenv SR /scr
    else
      setenv SR $TMPDIR
    endif
    breaksw
case 'IRIX64':
    if (-e /scr) then
      setenv SR /scr
    else
      setenv SR $TMPDIR
    endif
    breaksw
case 'unicos':
    if (-e /scr) then
      setenv SR /scr
    else
      setenv SR $TMPDIR
    endif
    breaksw
endsw
printenv SR
C
C --- E is expt, P is permenant directory, S is scratch directory.
C
cd ~/hycom/ATLb2.00/expt_02.0
setenv P $cwd
printenv P
setenv S `echo $P | sed -e "s?^.*$user?$SR/$user?"`
printenv S
setenv E `echo $P | awk '{print substr($0,length-3,2) substr($0,length,1)}'`
printenv E
C
/bin/ls -laFq
C
C --- check the RUNNING flag.
C
if ( -e RUNNING && ! -e RUNNING_$QSUB_REQID) then
C
C --- MODEL IS ALREADY RUNNING - EXIT.
C
  exit
endif
touch RUNNING
touch RUNNING_$QSUB_REQID
C
C --- Number of segments specified by '( SEG1 SEG2 ... )' in foreach.
C --- Use '( SEG1 )' for a single segment per run.
C
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
  /bin/rm -f RUNNING_$QSUB_REQID
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
    /bin/ls -laFq
    /bin/rm -f RUNNING
    /bin/rm -f RUNNING_$QSUB_REQID
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
set reqname = `echo ${QSUB_REQNAME} | awk '{print substr($0,2)}'`
/bin/ln -fs ${reqname}.log $script:r.log
C
mkdir -p $S
cp  ${SCRIPT} $S/${SCRIPT}
cd  $S
C
C ---------------------------------------------------------------------------
/bin/csh ${SCRIPT}
C ---------------------------------------------------------------------------
C
cd  $P
C
C --- Clean Up.
C
/bin/ls -laFq
C
/bin/mv LIST LISTOLD
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
/bin/rm -f RUNNING_$QSUB_REQID
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
        if ($OS == "Linux") then
#         fix a Linux NFS/tcsh bug
          (/bin/csh ${E}nqs.com |& cat - >& ${newname}.log) &
        else
#         /bin/nice -7 /bin/csh ${E}nqs.com >& ${newname}.log &
          /bin/csh ${E}nqs.com >& ${newname}.log &
        endif
    else
        qsub -mb -me -p 63 -r ${newname} \
                       -eo -o ${newname}.log ${E}nqs.com
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
