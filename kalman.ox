nldge::LinearLogLikeli()
{	decl zUp, zLeft,						/* to fill up state shock covariance matrix */
		 Z, Zc, 					   	/* measurement coefficients and constant term */
		 Tc, Q, H,						/* state equation coefficients and constant termn */
		 xt, Pt,						/* prior state & error covariance */
		 dif, Pnew,						/* to iterate to get start value for Pt */
		 loglikeli,
		 t, iss,						/* loop over observations */
		 xtt, Ptt,						/* posterior state, error covariance */
		 Ft, iFt, vt,						/* how to call that ??? */
		 mData,
		 ip,mxm=<>, mxsu=<>, mxsl=<>, mym=<>, mysu=<>, mysl=<>; 
									/* to get gives steady state, policy & transitions matrix */

	LinearMeasure(&Z, &Zc);					  	/* linearize measurement equations */
	mData = GetData(m_asY);
	
	Tc = m_vSss' - m_mP * m_vSss';					/* constant policy term */  
							  		/* fill up states shock covariance matrix for state without shocks */
	zUp   = zeros(m_cS - m_cSS	, m_cSS			); 	/* matrix on top of given shock covariance */
	zLeft = zeros(m_cS 			, m_cS - m_cSS 	); 	/* matrix to the left of both */												

	Q = ( zLeft ~ ( zUp | m_mSSbE * m_mSSbE' ));			/* state shock covariance matrix */
  	H = m_mMSbE * m_mMSbE';						/* measurement shocks covariance matrix */

	xt  = m_vSss';							/* start of posterior is steady state state */
	Pt  = unit(m_cS);				 	   	/* start for prediction variance at zeroth period */
//	do							  	/* iterate */
//	{	Pnew = m_mP*Pt*m_mP'+Q;					/* over this expression */
//    	dif  = norm(Pnew-Pt,0);						/* as long as infinity norm of change */
//    	Pt   = Pnew;
//	}while(dif > 10e-8);						/* is large */
	Pt   *= 1e-8;

	loglikeli = 0;							/* initialize loglikelihood */
	for(t = 0; t < sizer(mData); ++t)				/* iterate over all observations */
	{	xtt = m_mP*xt+Tc;					/* prior state */
    	Ptt = m_mP * Pt * m_mP' + Q;    
		Ft  = Z*Ptt*Z'+ H;					/* prior state covariance */ 

//		if(determinant(Ft)==0) return -.Inf;		

		iFt = invertgen(Ft,1);
		if(iFt==0) return -.Inf;
		vt  = mData[t][]' - Z*xtt - Zc;				/* inovation */
    	loglikeli -= 0.5*(m_cY * log(2*M_PI)+log(determinant(Ft))+vt'*iFt*vt);

		xt  = xtt + Ptt*Z'*iFt*vt;				/* posterior state */
    	Pt  = Ptt - Ptt*Z'*iFt*Z*Ptt;					/* posterior state error covariance */

	}
	return loglikeli;
}
