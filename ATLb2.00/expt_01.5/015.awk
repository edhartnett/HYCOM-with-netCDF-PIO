#
# --- awk script that reads a template model run shell script and produces
# --- the actual shell script for the part year run identified by
# --- y01 and ab.
#
# --- Usage:  awk -f 999.awk y01=002 ab=c 999.com > 999y002c.com
#
# --- Alan J. Wallcraft,  NRL,  August 1999.
#

#
# --- For a few hour/day run
# ---     Set np < 0, with the the number of hours per run = -np. 
# ---     Set ny = 1
# ---     Note that ab should be dXXXhYY where XXX is Julian day and YY is hour
# --- For part year runs:
# ---     Indicate the number of parts, np, to split the year into here.
# ---     Set ny = 1
# ---     Note that ab should be between a and the np-th letter alphabetically
# --- For whole year runs:
# ---     Indicate the number of years, ny, per run here.
# ---     Set np = 1
# ---     Note that ab is not used
# --- For actual years from 1901:
# ---     Set nd = 365
# --- This pattern might require changes for a new expt.
#
# --- Also set nd to the number of days per year (360 or 365 or 366).
# --- Note that nd == 365 indicates actual years from 1901
#

BEGIN { np = 1
        ny = 1
        nd = 360
	if      (np < 0)
		ti = -np/24.0
	else if (nd == 365 )
		ti = ny*366/np
	else
        	ti = ny*nd/np
	if      (np < 0)
		nq = 1
	else
		nq = np
	for ( i=0; i < nq; i++) {
		cb[sprintf("%c",i+97)] = sprintf("%c",i+98)
		ia[sprintf("%c",i+97)] = i
		}
	cb[sprintf("%c",nq+96)] = "a"
}

#
# --- supply a single line input "LIMITS"  to generate the standard limits file
# --- supply a single line input "LIMITI"  to generate the first    limits file
# --- supply a single line input "LIMITS?" to generate a std. ?-day limits file
# --- supply a single line input "LIMITI?" to generate a 1st  ?-day limits file
#
/^LIMIT[IS][0-9]*$/ {
		if ( np < 0 ) {
			ta  = substr(ab, 2, 3)
			tb  = substr(ab, 6, 2)
			}
		Y01 = int(y01)
		if (nd == 365) {
#			model day = wind day = days since Jan 1st 1901
			ty = nd*(Y01-1) + (Y01-1 - (Y01-1)%4)/4 + 1
			}
		else {
			ty = nd*(Y01-1)
			}
		if      ( np == 1) 
			ts = ty
		else if ( np < 0 )
			ts = ty + ta - 1 + tb/24.0
		else
			ts = ty + ia[ab]*ti

		td = substr($0, 7, length)
		if ( td == "" ) {
#			LIMITS or LIMITI
			if (nd == 365 && ia[ab] == np-1 && Y01%4 != 0)
				tm = ty + 365
			else
				tm = ts + ti
			}
		else {
#			LIMITS? or LIMITI?
			tm = ts + td
			}
		if ( substr($1,6,1) == "S")
			printf( " %14.5f %14.5f\n",  ts, tm )
		else
			printf( " %14.5f %14.5f\n", -ts, tm )
		next
}

#
# --- ab is a parameter passed in the call to this script,
# ---  it is a single lowercase letter indicating which part to run
# ---  or of the form dXXXhYY for few hour/day runs.
# --- Not used for whole year runs.
#
/^setenv A /  { 
		if      ( np == 1) 
			printf("setenv A \"\"\n")
		else if ( np < 0 )
			printf("setenv A \"%7s\"\n",ab)
		else
			printf("setenv A \"%1s\"\n",ab)
		next
}

/^setenv B /  { 
		if     ( np == 1)  {
			printf("setenv B \"\"\n")
			YXX = int(y01) + ny
			}
		else if ( np < 0 ) {
			ta  = substr(ab, 2, 3)
			tb  = substr(ab, 6, 2)
			tc  = tb - np
			if     (tc < 24) {
				tb = tc
				}
			else {
				ta = ta + tc/24
				tb = tc % 24
				}
			ylen = 365 + (4 - (int(y01) % 4))/4
			if      (ta > ylen) {
				ta  = ta - ylen
				YXX = int(y01) + 1
				}
			else {
				YXX = int(y01)
				}
			printf("setenv B \"d%3.3dh%2.2d\"\n",ta,tb)
			}
		else {
			printf("setenv B \"%1s\"\n",cb[ab])
			if (cb[ab] == ab) {
				YXX = int(y01) + 1
				}
			else {
				YXX = int(y01)
				}
			}
		next
}


#
# --- y01 is a parameter passed in the call to this script,
# ---  it represents the initial year of the model run.
#
/^setenv Y01/ { Y01 = int(y01)
		printf("setenv Y01 \"%03d\"\n",Y01)
		next
}

#
# --- this line must come after setenv B (see above)
#
/^setenv YXX/ { printf("setenv YXX \"%03d\"\n",YXX)
		next
}


#
# --- Detect bad use of the "C" comment command.
#
/^C .*[^\\][()><{}[\]]/ { print "C --- BAD COMMENT LINE DELETED."
		next
}


#
# --- All the above actions end with next, so we only get here if none
# --- of the above are true.
#
	      { print
}
