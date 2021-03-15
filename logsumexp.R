logsumexp = function(x,dimension=c(1,2)){
# % Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
# %   By default dim = 1 (row).
	if (is.null(dimension)) dimension=1 ;
# % subtract the largest in each row
	y = apply(x,dimension,max)
	x = x-y#sweep(x,dimension,y,"-")
	#s = y+log(rowSums(exp(x)))
	s = y+ log(apply(exp(x), dimension, sum))
	i = is.infinite(y)
	if (sum(i) > 0) s[i]=y[i]
	return(s)
}