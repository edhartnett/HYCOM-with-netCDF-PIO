#!/bin/csh
#
#@ job_name         = 016y001T_pipe
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
#@ wall_clock_limit = 1:00:00
#@ account_no       = NRLSS018
#@ class            = standard
#@ queue
#
set echo
set time = 1
set timestamp
C
C --- setup a master/slave pair of twin HYCOM simulations.
C
cd /scr/${user}/hycom/ATLb2.00/expt_01.6
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
cd ~/hycom/ATLb2.00/expt_01.6
mkdir dataT01
mkdir dataT03
C
csh PIPE/016y001T01.com >&! PIPE/016y001T01.log &
C
csh PIPE/016y001T03.com >&! PIPE/016y001T03.log &
C
C --- wait for them to end.
C
wait
cd /scr/${user}/hycom/ATLb2.00/expt_01.6
/bin/rm -f hycom_pipe
/bin/rm -f dataT01/PIPE_MASTER
/bin/rm -f dataT03/PIPE_SLAVE
