#
# --- awk script to provide the interface depth at a particular latitude
# --- usage: awk -f sig_lat.awk lat=20.0 sig_21.0
#

/^#/	{
	first = 0
	end = 0
	next
	}

	{
	if (first == 0) {
		first = 1
		r = $1 + 0.25
		l = -90.0
		d = $3
		}
	if ($2 - lat >= 0.0) {
		s = (lat - l) / ($2 - l)
		d = s * $3 + (1 - s) * d
		printf("%5.1f %6.1f %6.1f\n",r,lat,d)
		end = 1
		exit
		}
	l = $2
	d = $3
	}

END	{
	if (end == 0) {
		printf("%5.1f %6.1f %6.1f\n",r,lat,d)
		}
	}
	
