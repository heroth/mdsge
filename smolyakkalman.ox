nldge::GaussLogLikeli()
{
	decl it, ip, loglikeli, mdata,
		 zUp, zLeft, Q, H, dH, miH,
		 mxm=<>, mxsu=<>, mxsl=<>, mym=<>, mysu=<>, mysl=<>,
		 vL, mPxnew, dif,
		 vWt, vWm,
		 vMx11, mPx11,
		 	vMxe11, mPxC11, mPxeC11, mNxe11, mNx11, mNe11, mNp11,
		 		mNxt1, vMxt1, mDxt1, mPxt1,
			vMxmt1, mPxCt1, mPxmCt1, mNxmt1, mNmt1, mNpt1, // mPxt1
		 		mNyt1,	vMyt1, mDyt1, //mDxt1
  		 dPyt1,
		 	mPyt1, mPxyt1,
		 miPyt1, mk, vMxtt, mPxtt, dLPy, dLPys,
		 mZ, mZc, vpc;
	if(!m_bFilterGaussNonLin)
	{	LinearMeasure(&mZ, &mZc);
		vpc = m_vSss - (m_mP * m_vSss') ';			/* constant policy term */  

	}
		 
	mdata = GetData(m_asY);

	vWt  = m_oIntTime.GetWeights();					// prior weights
	vWm  = m_oIntMeas.GetWeights();					// prior weights
	
									/* fill up states shock covariance matrix for state without shocks */
	zUp   = zeros(m_cS - m_cSS	, m_cSS			); 	/* matrix on top of given shock covariance */
	zLeft = zeros(m_cS 			, m_cS - m_cSS 	); 	/* matrix to the left of both */												
	Q 	  = ( zLeft ~ ( zUp | m_mSSbE * m_mSSbE' ));		/* state shock covariance matrix */

	H   = m_mMSbE .^2;						// measurement shocks cov
	dH  = determinant(H);
	miH = invert(H);
	
	mPx11  = unit(m_cS);				 	   	// state cov prior
//	do											
//	{	mPxnew = m_mP*mPx11*m_mP'+Q;			
//    	dif    = norm(mPxnew-mPx11,0);			
//    	mPx11  = mPxnew;
//	}while(dif > 10e-8);						
	mPx11    *= 1e-8;
	vMx11	 = m_vSss;						// state prior

	loglikeli = 0;										
	for(it = 0; it < sizer(mdata); it++)
	{
									// time update
										// moments of nodes
		vMxe11  = vMx11 ~ zeros(1, m_cSS);						// xe _t-1|t-1	joint mean of states & shocks
		mPxC11  = choleski(mPx11)';							// Px _t-1|t-1	cholesky of state cov
		if(mPxC11==0) return -.Inf;							// 				stop if decomposition failed
		mPxeC11 = diagcat(mPxC11, m_mSSbE);						// Pxe_t-1|t-1  joint cov of state & shock
		mNxe11 = m_oIntTime.GetNodes(vMxe11,mPxeC11);			// generate x,e,w ~ N( xe^_t-1|t-1, Pxe_t-1|t-1
		mNx11   = mNxe11[][0:m_cS-1];							// xi _t-1|t-1 state		
		mNe11   = mNxe11[][m_cS:];			    				// e  _t-1|t-1 shocks																																			
		if(!m_bFilterGaussNonLin)
		{	//mNp11   = m_vXss + (m_mC*(mNx11-m_vSss)') ' ;		// vxf	= m_vXss + ( m_mC*(vsf- m_vSss)' )'; 
			mNxt1   = vpc + ( m_mP * mNx11') ';			// vsf = m_vSss + ( m_mP*(vs - m_vSss)' )' + vss;
			vMxt1 = vWt' * mNxt1;							// x^ _t  |t-1 mean
			mDxt1 = mNxt1 - vMxt1;							// Px _t  |t-1 cov transition
			for(ip=mPxt1=0; ip<rows(vWt); ip++) mPxt1 += vWt[ip]*(mDxt1[ip][]' * mDxt1[ip][]);
			mPxt1 += Q;
		}
		else
		{	mNp11   = m_oApprox.Interpolate(mNx11);			// pi _t-1|t-1 policy
										// time update
			fg(&mNxt1, mNx11, mNp11, mNe11);				// x  _t  |t-1 transition
			vMxt1 = vWt' * mNxt1;						// x^ _t  |t-1 mean
			mDxt1 = mNxt1 - vMxt1;						// Px _t  |t-1 cov transition
			for(ip=mPxt1=0; ip<rows(vWt); ip++) mPxt1 += vWt[ip]*(mDxt1[ip][]' * mDxt1[ip][]);
		}
									// measurement update
											// moments of nodes
		vMxmt1  = vMxt1 ~ zeros(1, m_cMS);					// xm _t|t-1	joint mean of states & shocks
		mPxCt1  = choleski(mPxt1)';						// Px _t|t-1	cholesky of state cov
		if(mPxCt1==0) return -.Inf;						// 				stop if decomposition failed
		mPxmCt1 = diagcat(mPxCt1, m_mMSbE);					// Pxm_t|t-1  joint cov of state & shock
		mNxmt1  = m_oIntMeas.GetNodes(vMxmt1,mPxmCt1);		// generate x,e,w ~ N( xe^_t-1|t-1, Pxe_t-1|t-1
		mNxt1   = mNxmt1[][0:m_cS-1];						// xi _t-1|t-1 state		
		mNmt1   = mNxmt1[][m_cS:];			    			// e  _t-1|t-1 shocks
		if(!m_bFilterGaussNonLin)
		{
			//mNpt1   = m_vXss + (m_mC*(mNxt1-m_vSss)') ' ;	// vxf	= m_vXss + ( m_mC*(vsf- m_vSss)' )'; 
			mNyt1   = (mZc   + (mZ * mNxt1') )'	 + mNmt1 ;		// vyf = ( mZc + mZ * vsf' )' + vms;						
		}
		else
		{
			mNpt1   = m_oApprox.Interpolate(mNxt1);	// pi _t-1|t-1 policy
									// measurement update
			fy(&mNyt1, mNxt1, mNpt1, mNmt1);		   	// yi  _t|t-1 measurement update
		}
		vMyt1   = vWm'  * mNyt1;					// y^  _t|t-1 mean of observation
		mDyt1   = mNyt1 - vMyt1;					// 			  deviation observation
		mDxt1	= mNxt1 - vMxt1;					// 			  deviation state
		for(ip=mPyt1=mPxyt1=0; ip<rows(vWm); ip++)
		{	mPyt1  += vWm[ip]*(mDyt1[ip][]' * mDyt1[ip][]);// Py  _t|t-1 covariance of observation
			mPxyt1 += vWm[ip]*(mDxt1[ip][]' * mDyt1[ip][]);// Pxy _t|t-1 cross covariance
		}
		miPyt1 = invertgen(mPyt1,1);


		mDyt1  = mdata[it][] - vMyt1;					// inovation  == mdata[it][] - (mZ*vMxt1' + mZc)' for linear estimation
		dLPy   = logdet(mPyt1, &dLPys);
		if(dLPys==0 || dLPys==-2 || dLPys==2) 	
			return -.Inf;
		else									dLPy *= dLPys;
		    loglikeli -= 0.5*(m_cY * log(2*M_PI)+dLPy + mDyt1 * miPyt1 * mDyt1');

		    mk     = mPxyt1 * miPyt1;					// K  _t 	kalman gain
		    vMxtt  = vMxt1  + (mdata[it][] - vMyt1)*mk';		// mu _t|t  mean of importance density
		    mPxtt  = mPxt1  - mk*mPyt1*mk';				// S  _t|t  cov  of importance density
	
		    vMx11 = vMxtt;
		    mPx11 = mPxtt;

		    mxm |= vMxtt; 						// remember for graph
		    mxsu|= vMxtt + sqrt(diagonal(mPxtt));
		    mxsl|= vMxtt - sqrt(diagonal(mPxtt));
			
		    mNpt1   = m_oApprox.Interpolate(vMx11);						
		    fy(&mNyt1, vMx11, mNpt1, zeros(1, m_cMS));						
		    mym  |= mNyt1;
		    mysu |= mNyt1 + diagonal(choleski(mPyt1));
		    mysl |= mNyt1 - diagonal(choleski(mPyt1)); 

	}

	return loglikeli;
}
