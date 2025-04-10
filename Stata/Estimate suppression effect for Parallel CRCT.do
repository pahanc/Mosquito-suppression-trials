*Perform linear regressions on mean counts per compound, month and cluster to estimate the suppression effect

*For each value of the number of clusters per arm 
foreach c_per_arm in 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 {

*For each replicate
forval i=1/100 {
	display `i'
	local pathname = "Sim_dat_Apr_para_2yr_"
	local pathname2 = "\counts_df_"
	local filename = ".csv"  
	local pathfilename= "`pathname'"+"`c_per_arm'"+"`pathname2'"+"`c_per_arm'"+"_"+"`i'"+"`filename'"
	local outfilename="`pathname'"+"`c_per_arm'"+"\"+"p_values_"+"`c_per_arm'"+".txt"
	import delimited "`pathfilename'", clear
	collapse (mean) treatment count, by(cluster period)
	gen istrial = 1 if period>1
	replace istrial=0 if period==1
	collapse (mean) treatment count, by(cluster istrial)
	gen count_bline=.
	sum cluster
	forval k = 1/`r(N)' {
		local clust = cluster[`k']
		if (istrial[`k']==0)	{
			local val=count[`k'] 
		}
		replace count_bline=`val' if cluster==`clust' & istrial==1
	}
	drop if istrial==0
	gen log_count=log(count)
	gen log_count_bline=log(count_bline)
	regress log_count log_count_bline i.treatment
	file open myfile using "`outfilename'", write append
	local t = _b[1.treatment]/_se[1.treatment]
	local p = 2*ttail(e(df_r),abs(`t'))
	set more off
	file write myfile ("`p'") _n 
	file close myfile
	set more off
}

}
