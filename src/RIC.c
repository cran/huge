void RIC(double *X, int *dd, int *nn, double *lambda_opt)
{
	int d,n;
	int i,j,k,m;
	
	d = dd[0];
	n = nn[0];
	
	double lambda_min,lambda_max,tmp;
	
	lambda_min = 99999999;
	
	for(j=0;j<d;j++)
	{	
		lambda_max = 0;
		for(i=1;i<n;i++)
		{
			for(k=0;k<j;k++)
			{
				tmp = 0;
				for(m=0;m<(n-i);m++)
					tmp = tmp + X[j*n+m+i]*X[k*n+m];
				for(m=(n-i);m<n;m++)
					tmp = tmp + X[j*n+m-(n-i)]*X[k*n+m];
				if(tmp<0)
					tmp = -tmp;
				if(tmp>lambda_max)
					lambda_max = tmp;
			}
			for(k=j+1;k<d;k++)
			{
				tmp = 0;
				for(m=0;m<(n-i);m++)
					tmp = tmp + X[j*n+m+i]*X[k*n+m];
				for(m=(n-i);m<n;m++)
					tmp = tmp + X[j*n+m-(n-i)]*X[k*n+m];
				if(tmp<0)
					tmp = -tmp;
				if(tmp>lambda_max)
					lambda_max = tmp;
			}
			
		}
		if(lambda_max<lambda_min)
			lambda_min = lambda_max;
	}
	lambda_opt[0] = lambda_min;
}
