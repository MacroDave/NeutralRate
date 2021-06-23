subroutine MUE_Stage2(matrix x_mue2)

	smpl @all

	!T = @rows(x_mue2)
	!N = @columns(x_mue2)

	matrix x = @subextract(x_mue2,1,2,!T,!N)
	matrix y = @subextract(x_mue2,1,1,!T,1)

	matrix(!T-2*4+1) stat = 0
	matrix(!T,1) xbind = 1

    for !i = 4 to !T-4
		vector(!i) tempxr = 0
		matplace(xbind,tempxr,1,1)
		matrix xr = @hcat(x, xbind)
		!NN = @columns(x_mue2)
      	matrix xi = @inverse(@transpose(xr)*xr)
      	matrix b  = @inverse(@transpose(xr)*xr)*@transpose(xr)*y
      	scalar s3 = @csum(@epow(y-xr*b,2))/(!T-@columns(xr))
		vector tempv = b(!NN)/@sqrt(s3*xi(!NN,!NN))
		matplace(stat,tempv,!i+1-4,1)
    next

    scalar ew = 0
    for !i = 1 to @rows(stat)
        ew = ew+exp(@epow(stat(!i),2)/2)
    next
    scalar ew  = log(ew/!lengthstat)
    scalar mw  = @csum(@epow(stat,2)) / @rows(stat)
    scalar qlr = @max(@epow(stat,2))

    'Values are from Table 3 in Stock and Watson (1998)
    'Test Statistic: Exponential Wald (EW)
	vector(31) valew
	valew.fill 0.426, 0.476, 0.516, 0.661, 0.826, 1.111, _
               1.419, 1.762, 2.355, 2.91,  3.413, 3.868, 4.925, _
               5.684, 6.670, 7.690, 8.477, 9.191, 10.693, 12.024, _
               13.089, 14.440, 16.191, 17.332, 18.699, 20.464, _
               21.667, 23.851, 25.538, 26.762, 27.874
    
	'Test Statistic: Mean Wald (MW)
	vector(31) valmw
    valmw.fill 0.689, 0.757, 0.806, 1.015, 1.234, 1.632, _
               2.018, 2.390, 3.081, 3.699, 4.222, 4.776, 5.767, _
               6.586, 7.703, 8.683, 9.467, 10.101, 11.639, 13.039, _
               13.900, 15.214, 16.806, 18.330, 19.020, 20.562, _
               21.837, 24.350, 26.248, 27.089, 27.758

	'Test Statistic: QLR
	vector(31) valql
    valql.fill 3.198, 3.416, 3.594, 4.106, 4.848, 5.689, _
               6.682, 7.626, 9.16,  10.66, 11.841, 13.098, 15.451, _
               17.094, 19.423, 21.682, 23.342, 24.920, 28.174, 30.736, _
               33.313, 36.109, 39.673, 41.955, 45.056, 48.647, 50.983, _
               55.514, 59.278, 61.311, 64.016
    
    'Median-unbiased estimator of lambda_g for given values of the test
    'statistics are obtained using the procedure described in the 
    'footnote to Stock and Watson (1998) Table 3.
	scalar lame=NA
	scalar lamm=NA
	scalar lamq=NA

    if (ew <= valew(1)) then
        scalar lame = 0
    else
        for !i=1 to 31-1
            if ew > valew(!i) then
			if ew<=valew(!i+1) then
				scalar lame = !i-1+(ew-valew(!i))/(valew(!i+1)-valew(!i))
           	endif 
		 endif
        next
    endif

    if (mw <= valmw(1)) then
        lamm = 0
    else
        for !i=1 to 31-1
            if mw > valmw(!i) then
			if mw <= valmw(!i+1) then
                	scalar lamm = !i-1+(mw-valmw(!i))/(valmw(!i+1)-valmw(!i))
            	endif
        	   endif
    		next
	endif

    if qlr <= valql(1) then
        lamq = 0
    else
        for !i=1 to 31-1
            if qlr > valql(!i) then
			if qlr <= valql(!i+1) then
	                scalar lamq = !i-1+(qlr-valql(!i))/(valql(!i+1)-valql(!i))
     		endif
		  endif
        next
    endif

    if (@isna(lame) or @isna(lamm) or @isna(lamq)) then
		@uiprompt("At least one statistic has an NA value. Check to see if your EW, MW, and/or QLR value is outside of Table 3.")
    endif
    
	vector(3) stats
		stats(1) = ew
		stats(2) = mw
		stats(3) = qlr
	vector(3) lams
		lams(1) = lame
		lams(2) = lamm
		lams(3) = lamq

	scalar lambda_z = lame/!T

	'Cleanup
	delete(noerr) lame lamm lamq mw qlr s3 stat tempv tempxr valew valmw valql xi xr y x xbind

endsub


