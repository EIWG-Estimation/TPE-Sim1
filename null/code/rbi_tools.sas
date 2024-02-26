/*******************************************************************************
| Name       : rbi_tools.sas
| Purpose    : FCMP Functions for imputing in the datastep
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 21JUN21 
|-------------------------------------------------------------------------------
| FCMP Call Routine List:
|--------------------------------------------------------------------------------
|
| Name     : cond_impute_mvn(y,m,v,i) 
| Purpose  : Conditional Imputation of missing y vector values using the MVN  
| Arguments: y [REQ] = 1 by n defined array with observed and missing data 
|            m [REQ] = 1 by n defined array with the means to use 
|            v [REQ] = n by n defined array with the vc matrix to use 
|            i [REQ] = 1 by n defined array to store the imputed flags
|  
***********************************************************************************/;

proc fcmp outlib = work.functions.rbi_tools;

  subroutine cond_impute_mvn( y[*], m[*], v[*,*], i[*]);   
    outargs y, i;
    length msg1 msg2 $200;

	*** CHECK DIMENSIONS MATCH ***;
	ydim  = dim1(y);
	mdim  = dim1(m);
	vdim1 = dim1(v);
	vdim2 = dim2(v);
	idim  = dim1(i);
    if not ( ydim = mdim = vdim1 = vdim2 = idim ) then do;
       msg1 = "ER"||upcase("ror:(FCMP):")||"The Function COND_IMPUTE_MVN does not have matching array sizes.";
       msg2 = "ER"||upcase("ror:(FCMP): Dimensions:");
       put msg1;
       put msg2 ydim1= mdim1= vdim1= vdim2=;
	end;
	else do;
	
	  *** WORK OUT THE SIZE OF OBSERVED (N1) AND MISSING (N2) PARTITIONS ***;
      do ii = 1 to dim(y);
	    if y[ ii ] ne . then n1 + 1;
	  end;
      n2 = dim(y) - n1;
	
      if n1 = dim(y) then do;

        msg1 = "WA"||upcase("rning:(FCMP):")||"The Function COND_IMPUTE_MVN does not have missing values in Y array. No imputation performed.";
        put msg1;

      end;
      else if n1 = 0 then do;

        msg1 = "NO"||upcase("te:(FCMP):")||"Values supplied to Y array in COND_IMPUTE_MVN are all missing. MVN data simulated based on M and V arrays.";
        put msg1;
     
  	    *** SIMULATE VALUE FROM MVN ***;
        array zvec[1,1] / nosymbols;
        array smat[1,1] / nosymbols;
        array sxz [1,1] / nosymbols;
        array imp [1,1] / nosymbols;
        call dynamic_array(zvec,n2,1);
        call dynamic_array(smat,n2,n2);
        call dynamic_array(sxz,n2,1);
        call dynamic_array(imp,n2,1);

        do ii = 1 to dim(y);
          zvec[ ii,1 ] = rand("NORMAL",0,1);
	    end;

        call chol(v,smat);
        call mult(smat,zvec,sxz);
        call addmatrix(m,sxz,imp);
 
        *** CREATE ZEROS FOR ALL IMPUTATION VARIABLES ***;
        do ii = 1 to dim(y);
          i[ ii ] = 0;
	    end;

        *** UPDATE Y VECTOR WITH IMPUATIONS AND FLAG THE IMPUTATION VARIABLE ***;
        do ii = 1 to dim(y);
          y[ ii ] = imp[ ii,1 ];
          i[ ii ] = 1;
	    end;

      end;
      else do;

	    *** CREATE LOCATION ARRAYS FOR OBSERVED (L1) AND MISSING (L2) VALUES ***;
        array l1[1] / nosymbols;
        array l2[1] / nosymbols;
        call dynamic_array( l1, n1 );
        call dynamic_array( l2, n2 );

	    *** FILL LOCATION ARRAYS WITH THE INDEX FOR EACH OBSERVED AND MISSING VALUE ***;
        ncomp = 0;
	    nmiss = 0;
        do ii = 1 to dim(y);
	      if y[ ii ] ne missing then do;
            ncomp + 1;
		    l1[ ncomp ] = ii;
          end;
		  else do;
            nmiss + 1;
		    l2[ nmiss ] = ii; 
		  end;
	    end;

        *** CREATE PARTITION ARRAYS AND FILL USING LOCATION ARRAYS ***;
        array y1  [1,1] / nosymbols;
        array y2  [1,1] / nosymbols;
        array m1  [1,1] / nosymbols;
        array m2  [1,1] / nosymbols;
        array v11 [1,1] / nosymbols;
        array v12 [1,1] / nosymbols;
        array v21 [1,1] / nosymbols;
        array v22 [1,1] / nosymbols;

        call dynamic_array( y1, n1, 1 );
        call dynamic_array( y2, n2, 1 );
        call dynamic_array( m1, n1, 1 );
        call dynamic_array( m2, n2, 1 );
        call dynamic_array( v11, n1, n1 );
        call dynamic_array( v12, n1, n2 );
        call dynamic_array( v21, n2, n1 );
        call dynamic_array( v22, n2, n2 );

	    do ii = 1 to n1;
          y1[ ii ] = y[ l1[ii] ];
          m1[ ii ] = m[ l1[ii] ];
        end;
        do ii = 1 to n2;
          y2[ ii ] = y[ l2[ii] ];
          m2[ ii ] = m[ l2[ii] ];
        end;
	    do ii = 1 to n1;
 	      do jj = 1 to n1;
            v11[ ii,jj ] = v[ l1[ii],l1[jj] ];
	      end;
          do jj = 1 to n2;
            v12[ ii,jj ] = v[ l1[ii],l2[jj] ];
	      end;
	    end;
	    do ii = 1 to n2;
 	      do jj = 1 to n1;
            v21[ ii,jj ] = v[ l2[ii],l1[jj] ];
	      end;
          do jj = 1 to n2;
            v22[ ii,jj ] = v[ l2[ii],l2[jj] ];
	      end;
	    end;

        *** CREATE CONDITIONAL MEANS USING FORMULA M_2|1 = M2 + V21*INV(V11)*(Y-M1) ***;
        array residual[1,1] / nosymbols;
	    array invsigma[1,1] / nosymbols;
        array stdresid[1,1] / nosymbols;
        array condinfo[1,1] / nosymbols;
        array cond_mvn_mu[1,1] / nosymbols;
        call dynamic_array(residual,n1,1);
        call dynamic_array(invsigma,n1,n1);
        call dynamic_array(stdresid,n1,1);
        call dynamic_array(condinfo,n2,1);
        call dynamic_array(cond_mvn_mu,n2,1);

        call subtractmatrix(y1,m1,residual);
        call inv(v11,invsigma);
	    call mult(invsigma,residual,stdresid);
	    call mult(v21,stdresid,condinfo);
	    call addmatrix(m2,condinfo,cond_mvn_mu);

        *** CREATE CONDITIONAL VC MATRIX USING FORMULA V_2|1 = V22 - V21*INV(V11)*V12 ***;
        array matrix1 [1,1] / nosymbols;
        array matrix2 [1,1] / nosymbols;
        array cond_mvn_vc[1,1] / nosymbols;
        call dynamic_array(matrix1,n1,n2);
        call dynamic_array(matrix2,n2,n2);
        call dynamic_array(cond_mvn_vc,n2,n2);
  
        *** SOMETIMES ROUNDING ERRORS SO MAY NEED TO AVERAGE THE OFF DIAGONALS ***;
        call mult(invsigma,v12,matrix1);
        call mult(v21,matrix1,matrix2);
	    call subtractmatrix(v22,matrix2,cond_mvn_vc);

	    *** SIMULATE VALUE FROM THE CONDITIONAL MVN ***;
        array zvec[1,1] / nosymbols;
        array smat[1,1] / nosymbols;
        array sxz [1,1] / nosymbols;
        array imp [1,1] / nosymbols;
        call dynamic_array(zvec,n2,1);
        call dynamic_array(smat,n2,n2);
        call dynamic_array(sxz,n2,1);
        call dynamic_array(imp,n2,1);

        do ii = 1 to n2;
          zvec[ ii,1 ] = rand("NORMAL",0,1);
	    end;

        call chol(cond_mvn_vc,smat);
        call mult(smat,zvec,sxz);
        call addmatrix(cond_mvn_mu,sxz,imp);
 
        *** CREATE ZEROS FOR ALL IMPUTATION VARIABLES ***;
        do ii = 1 to dim(y);
          i[ ii ] = 0;
	    end;

        *** UPDATE Y VECTOR WITH IMPUATIONS AND FLAG THE IMPUTATION VARIABLE ***;
        do ii = 1 to n2;
          y[ l2[ ii ] ] = imp[ ii,1 ];
		  i[ l2[ ii ] ] = 1;
	    end;

      end;
    end;
  endsub;

quit;

options cmplib = work.functions; 

