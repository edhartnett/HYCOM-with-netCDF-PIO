#!/bin/csh
#
#@ job_name         = 024y001T_pipe
#@ output           = $(job_name).log
#@ error            = $(job_name).log
#@ restart          = yes
#@ job_type         = parallel
#@ network.MPI      = css0,not_shared,US
#@ environment      = MP_EUILIB=us     
#@ node             = 1
#@ total_tasks      = 4
#@ node_usage       = not_shared
#@ wall_clock_limit =  1:00:00
#@ account_no       = NRLSS018
#@ class            = batch
#@ queue
#
set echo
set time = 1
set timestamp
C
C --- setup a master/slave pair of twin HYCOM simulations.
C
cd /scr/${user}/hycom/ATLa2.00/expt_02.4
/bin/rm -f hycom_pipe
mknod      hycom_pipe p
C
C --- master
C
/bin/rm -f              dataT01/PIPE_base.out
echo "../hycom_pipe" >! dataT01/PIPE_MASTER
C
C --- slave
C
/bin/rm -f              dataT03/PIPE_test.out
echo "../hycom_pipe" >! dataT03/PIPE_SLAVE
C
C --- start both simulations, master first.
C
cd ~/hycom/ATLa2.00/expt_02.4
mkdir dataT01
mkdir dataT03
C
csh 024y001T01.com >&! 024y001T01.log &
C
csh 024y001T03.com >&! 024y001T03.log &
C
C --- wait for them to end.
C
wait
cd /scr/${user}/hycom/ATLa2.00/expt_02.4
/bin/rm -f hycom_pipe
/bin/rm -f dataT01/PIPE_MASTER
/bin/rm -f dataT03/PIPE_SLAVE
