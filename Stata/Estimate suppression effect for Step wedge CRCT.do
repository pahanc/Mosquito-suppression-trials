ssc install swpermute

*Perform linear regressions on mean counts per compound, month and cluster to estimate the suppression effect

*For each value of the number of clusters per sequence
foreach c_per_seq in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 {
	
*For each replicate
forval i=1/100 {
	display `i'
	local pathname =  "Sim_dat_Apr_"
	local pathname2 = "\counts_df_"
	local filename = ".csv"  
	local pathfilename= "`pathname'"+"`c_per_seq'"+"`pathname2'"+"`c_per_seq'"+"_"+"`i'"+"`filename'"
	local outfilename="Sim_dat_Apr_"+"`c_per_seq'"+"\"+"p_values_"+"`c_per_seq'"+".txt"
	import delimited "`pathfilename'", clear
	collapse (mean) count (count) N=count, by(cluster period treatment)
	gen istrial = 1 if period>1
	replace istrial=0 if period==1
	gen count_bline=.
	sum cluster
	forval k = 1/`r(N)' {
		*display `k'
		local clust = cluster[`k']
		if (istrial[`k']==0)	{
			local val=count[`k'] 
		}
		*display `val'
		replace count_bline=`val' if cluster==`clust' & istrial==1
	}
	drop if istrial==0
	gen log_count = log(count)
	gen log_count_bline=log(count_bline)
	set seed 72906
	swpermute _b[1.treatment], cluster(cluster) period(period) intervention(treatment) withinperiod weightperiod(variance _se[1.treatment]^2) reps(500) nodots:    regress log_count log_count_bline i.treatment
	file open myfile using "`outfilename'", write append
	matrix test = r(p)
	local p=el(test,1,4)
	set more off
	file write myfile ("`p'") _n 
	file close myfile
	set more off
}

}