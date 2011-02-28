void SPMBscr(double *S, int *idx_scr, double *lambda, int *nnlambda, int *dd, int *nnscr, double *x, int *col_cnz, int *row_idx)
{
	int d,nlambda,nscr;
    int m,i,j,k;
	
    int junk_i,resize_i,tmp5,tmp6,tmp7,cnz;
    int w_idx,rss_idx,size_a,size_i;

    double ilambda;
    double tmp1,tmp2,tmp3,tmp4;
            
    int iter_ext,iter_int;
    double move_ext,move_int;
    double gap_ext,gap_int;
    double thol;
    
    thol = 1e-4;
	
	d = dd[0];
	nlambda = nnlambda[0];
    nscr = nnscr[0];
	
	double w0[d],w1[d];
	int idx_a[nscr],idx_i[nscr]; 
	
	cnz = 0;
	
    for(m=0;m<d;m++)
    {
        size_i = nscr;
        size_a = 0;
        for(j=0;j<d;j++)
            w0[j] = 0;
        
        for(j=0;j<nscr;j++)
        {
            idx_i[j] = idx_scr[m*nscr+j];
            idx_a[j] = -1;
        }
		
		for(i=0;i<nlambda;i++)
        {
            ilambda = lambda[i];
            gap_ext = 9999;
            iter_ext = 0;
			tmp7 = 1;
            while(iter_ext<10000 && tmp7 >0)
            {
                iter_ext++;
                move_ext = 0;  
				tmp6 = size_a;
                for(j=0;j<size_i;j++)
                {
                    w_idx = idx_i[j];
                    tmp1 = S[w_idx*d+m] + w0[w_idx];
                    
                    tmp5 = w_idx*d;
                    for(k=0;k<size_a;k++)
                    {                   
                        rss_idx = idx_a[k];
                        tmp1 = tmp1 - S[tmp5+rss_idx]*w0[rss_idx];
                    }  
                    
                    if(tmp1 > ilambda)
                    {
                        w1[w_idx] = tmp1 - ilambda;
                        *(idx_a+size_a) = w_idx;
                        size_a++;
                        *(idx_i+j) = -1;
                    }
                    
                    else if(tmp1 <-ilambda)
                    {
                        w1[w_idx] = tmp1 + ilambda;
                        *(idx_a+size_a) = w_idx;
                        size_a++;
                        *(idx_i+j) = -1;
                    }
                    
                    else w1[w_idx] = 0;
                    
                    tmp2 = w1[w_idx] - w0[w_idx];
                    if(tmp2>0)
                        move_ext = move_ext + tmp2;
                    else move_ext = move_ext - tmp2;
                    
                    w0[w_idx] = w1[w_idx];
                }
                
                tmp7 = size_a - tmp6;
                
                junk_i = 0;
                resize_i = size_i;
                
                for(k=0;k<size_i;k++)
                {
                    if(idx_i[k]==-1)
                    {
                        junk_i = junk_i + 1;
                        resize_i--;
                    }
                    else *(idx_i+k-junk_i) = *(idx_i+k);
                }
                size_i = resize_i;
                gap_int = 9999;
                iter_int = 0;
                while(gap_int>thol && iter_int <10000)
                {
                    iter_int = iter_int + 1;
                    move_int = 0;
                    for(j=0;j<size_a;j++)
                    {
                        w_idx = idx_a[j];
						tmp1 = S[w_idx*d+m] + w0[w_idx];
                        
                        tmp5 = w_idx*d;
                        for(k=0;k<size_a;k++)
                        {                   
                            rss_idx = idx_a[k];
                            tmp1 = tmp1 - S[tmp5+rss_idx]*w0[rss_idx];
                        }
                        
                        if(tmp1 > ilambda)
                        {
                            w1[w_idx] = tmp1 - ilambda;
                        }
                        
                        else if(tmp1 <-ilambda)
                        {
                            w1[w_idx] = tmp1 + ilambda;
                        }                    
                        
                        else w1[w_idx] = 0;
                        
                        tmp2 = w1[w_idx] - w0[w_idx];
                        
                        if(tmp2>0)
                            move_int = move_int + tmp2;
                        else move_int = move_int - tmp2;
						
                        w0[w_idx] = w1[w_idx];
                    }
                    tmp3 = 0;
                    for(k=0;k<size_a;k++)
                    {                   
                        tmp4 = w0[idx_a[k]];
                        if(tmp4>0)
                            tmp3 = tmp3 + tmp4;
                        else tmp3 = tmp3 - tmp4;
                    }
                    gap_int = move_int/tmp3;
                }
			}
			for(j=0;j<size_a;j++)
			{
				w_idx = idx_a[j];
				x[cnz] = w1[w_idx];
				row_idx[cnz] = i*d+w_idx;
				cnz++;
			}
        }
		col_cnz[m+1]=cnz;
    }
}
