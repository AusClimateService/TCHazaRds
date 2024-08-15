#include <Rcpp.h>
using namespace Rcpp;

//' @title TC Track Distance and Direction From Output Grid Point
//' @description Grid point time series of TC distance and direction.
//' @param Gridlon single Grid point longitude
//' @param Gridlat single Grid point latitude
//' @param TClon vector of TC longitudes
//' @param TClat vector of TC latitudes
//' @return two columns for distance in km and cartesian direction in degrees, counterclockwise from the x axis.
//' //@example RdistPi(142,-14,c(144,145),c(-11,-12))
// [[Rcpp::export]]
NumericMatrix RdistPi(float Gridlon, float Gridlat, NumericVector TClon, NumericVector TClat)
{

  //Rlam is a combination of ...
  //R: radius in km from the centre of the TC to the gird point
  //lam: Direction(geographic bearing, positive clockwise) !!!isn't this the Cartesian direction, counter clockwise from the x axis!!!
  // Haversine formula constants
  const float Rearth = 6372797.560856f;
  const float pi = 3.141592f;
  const float piOn180 = pi / 180;
  // Convert latitude and longitude to radians
  float lat2 = Gridlat * piOn180;
  float lon2 = Gridlon * piOn180;

  // Vector to store results
  int n = TClon.size();
  NumericMatrix Rlam(n, 2);

  // Precompute constants used in the loop
  float cos_lat2 = cosf(lat2);
  float sin_lat2 = sinf(lat2);

  for (int i = 0; i < n; i++) {
    float lat1 = TClat[i] * piOn180;
    float lon1 = TClon[i] * piOn180;

    // Calculate differences in latitude and longitude
    float dlat = lat2 - lat1;
    float dlon = lon2 - lon1;

    // Haversine formula
    float a = sinf(dlat / 2.0f) * sinf(dlat / 2.0f) + cosf(lat1) * cos_lat2 * sinf(dlon / 2.0f) * sinf(dlon / 2.0f);
    float c = 2.0f * atan2f(sqrtf(a), sqrtf(1.0f - a));
    Rlam(i, 0) = c * Rearth / 1000.0f; // Convert to km

    // Calculate the bearing (direction)
    float x = sinf(dlon) * cos_lat2;
    float y = cosf(lat1) * sin_lat2 - sinf(lat1) * cos_lat2 * cosf(dlon);
    Rlam(i, 1) = atan2f(y, x) / piOn180;
  }

  return Rlam;
}

//' @title Jelesnianski Wind Profile Time Series
//' @description wind profile time series at a grid point
//' @param f single coriolis parameter at the centre of TC in hz
//' @param vMax maximum wind velocity calculation in m/s
//' @param rMax radius of maximum winds in km
//' @param R vector of distances from grid points to TC centre in km
//' @return array with two columns for velocity and then vorticity.
//' //@example JelesnianskiWindProfilePi(-1e-4,20,20,50)
// [[Rcpp::export]]
NumericMatrix JelesnianskiWindProfilePi(NumericVector f, NumericVector vMax, NumericVector rMax, NumericVector R)
{
	//
	int n = R.size();
	NumericMatrix VZ(n,2);

 	float Vi, Ri,fi, vMaxi,rMaxi,Zi, sf;


	for(int i = 0;i < n; i++){
		//
		fi = f[i];
		sf = (fi / fabs(fi));
		Ri = R[i];
	  vMaxi = vMax[i];
	  rMaxi = rMax[i];
		Vi = 2.0f * vMaxi * rMaxi * Ri / (rMaxi *rMaxi + Ri * Ri) * sf;
		VZ(i,0) = Vi;
		Zi = (sf * 2.0f * vMaxi * rMaxi / (rMaxi *rMaxi + Ri * Ri) + sf * 2.0f * vMaxi * rMaxi * (rMaxi *rMaxi - Ri * Ri) /((rMaxi *rMaxi + Ri * Ri) * (rMaxi *rMaxi + Ri * Ri)));
		VZ(i,1) = Zi;
	}
	return VZ;
}

//' @title Rankine Wind Profile Time Series
//' @description wind profile time series at a grid point
//' @param f single coriolis parameter at the centre of TC in hz
//' @param vMax maximum wind velocity calculation in m/s
//' @param rMax radius of maximum winds in km
//' @param R vector of distances from grid points to TC centre in km
//' @return array with two columns for velocity and then vorticity.
//' //@example RankineWindProfilePi(-1e-4,20,20,50)
// [[Rcpp::export]]
NumericMatrix RankineWindProfilePi(NumericVector f, NumericVector vMax, NumericVector rMax, NumericVector R)
{
  //
  int n = R.size();
  NumericMatrix VZ(n,2);

  float Vi, Ri,fi, vMaxi,rMaxi,Zi, sf;
  float alpha = 0.5;

  for(int i = 0;i < n; i++){
    //
    fi = f[i];
    sf = (fi / fabs(fi));
    Ri = R[i];
    vMaxi = vMax[i];
    rMaxi = rMax[i];
    if(Ri <= rMaxi){
      Vi = vMaxi*Ri/rMaxi;
      Zi = sf*((vMaxi*(Ri/rMaxi)+vMaxi/rMaxi));
    }
    else{
      Vi = vMaxi * powf(rMaxi / Ri,alpha)*sf;
      Zi = Vi/Ri-alpha*vMaxi*rMaxi/(powf(Ri,alpha+1.0f));
    }
    VZ(i,0) = Vi;
    VZ(i,1) = Zi;
  }
  return VZ;
}


//' @title Holland Wind Profile Time Series
//' @description wind profile time series at a grid point
//' @param f single coriolis parameter at the centre of TC in hz
//' @param vMax maximum wind velocity calculation in m/s
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param rho density of air in Kg/m3
//' @param beta exponential term for Holland vortex
//' @param R vector of distances from grid points to TC centre in km
//' @return array with two columns for velocity and then vorticity.
//' //@example HollandWindProfilePi(-1e-4,20,20,10,1.15,1.2,50)
// [[Rcpp::export]]
NumericMatrix HollandWindProfilePi(NumericVector f, NumericVector vMax, NumericVector rMax, NumericVector dP, float rho, NumericVector beta, NumericVector R)
{
	//Holland profile. For `r < rMax`, we reset the wind field to a
	//cubic profile to avoid the barotropic instability mentioned in
	//Kepert & Wang (2001).

	int n = R.size();
  NumericMatrix VZ(n,2);
  float Vi, Ri,fi,vMaxi,rMaxi,betai,dPi, Zi;
	float E, d2Vm,dVm,aa,bb,cc;
	float delta, edelta;
  E = expf(1.0f);

	for(int i = 0;i < n; i++){
		//
		Ri = R[i];
	  fi = f[i];
	  vMaxi = vMax[i];
	  rMaxi = rMax[i];
	  betai = beta[i];
	  dPi = dP[i]*100.0f;


		d2Vm = ((betai * dPi * (-4.0f * betai *betai *betai * dPi / rho - (-2.0f + betai *betai) * E * (fi * rMaxi) *(fi * rMaxi))) / (E * rho * sqrtf((4.0f * betai * dPi) / (E * rho) + (fi * rMaxi) *(fi * rMaxi)) * (4.0f * betai * dPi * rMaxi *rMaxi / rho + E * (fi * rMaxi *rMaxi) *(fi * rMaxi *rMaxi))));

    dVm = (-fabs(fi)/2.0f + (E * (fi*fi) * rMaxi * sqrtf((4.0f * betai * dPi / rho) / E + (fi * rMaxi) * (fi * rMaxi))) / (2.0f * (4.0f * betai * dPi / rho + E * (fi * rMaxi) * (fi * rMaxi))));

		aa = ((d2Vm / 2.0f - (dVm - vMaxi / rMaxi) / rMaxi) / rMaxi);

		bb = (d2Vm - 6.0f * aa * rMaxi) / 2.0f;

		cc = dVm -3.0f * aa * rMaxi * rMaxi - 2.0f * bb * rMaxi;

		delta = powf(rMaxi / Ri, betai);
		edelta = expf(-delta);
		if (Ri <= rMaxi)
		{
			Vi = (Ri * (Ri * (Ri * aa + bb) + cc));
			Zi = Ri * (Ri * 4.0f * aa + 3.0f * bb) + 2.0f * cc;
		}
		else
		{
		  //dPi = dPi*100.0f; // !need to multiply by dPi 100 to get to Pa not in TCRM  https://github.com/GeoscienceAustralia/tcrm/blob/1916233f7dfdecf6a1b5f5b0d89f3eb1d164bd3e/wind/windmodels.py#L381
			Vi =  (sqrtf((dPi * betai / rho) * delta * edelta + (Ri * fi / 2.0f)*(Ri * fi / 2.0f)) - Ri *fabs(fi) / 2.0f);
			//Zi = ((sqrtf((dPi * betai / rho) * delta * edelta + (Ri * fi / 2.0f)*(Ri * fi / 2.0f))) / Ri - fabs(fi) + edelta * (2.0f * (betai * betai) * dPi * (delta - 1.0f) * delta + rho * edelta * (fi * Ri) *(fi * Ri)) / (2.0f * rho * Ri * sqrtf(4.0f * (betai * dPi / rho) * delta * edelta + (fi * Ri) *(fi * Ri))));
		  Zi = fabs(fi) + (betai*betai * dPi * (delta * delta) * edelta / (2.0f * rho * Ri) - betai*betai * dPi * delta * edelta / (2.0f * rho * Ri) + Ri * fi * fi / 4.0f) / sqrtf(betai * dPi * delta * edelta / rho + (Ri * fi / 2)*(Ri * fi / 2)) + (sqrtf(betai * dPi * delta * edelta / rho + (Ri * fi / 2)*(Ri * fi / 2))) / Ri;
		}
		VZ(i,0) = Vi*fi / fabs(fi);
		VZ(i,1) = Zi*fi / fabs(fi);
	}
  return VZ;
}

//' @title Holland Pressure Profile Time Series
//' @description Pressure profile time series at a grid point.
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param cP TC central pressure in hPa
//' @param beta exponential term for Holland vortex
//' @param R vector of distances from grid points to TC centre in km
//' @return vector of pressures.
//' //@example HollandPressureProfilePi(20,20,980,1.2,50)
// [[Rcpp::export]]
NumericVector HollandPressureProfilePi(NumericVector rMax, NumericVector dP, NumericVector cP, NumericVector beta, NumericVector R)
{
	//Holland pressure profile
	int n = R.size();
  NumericVector P(n);


	float Ri,rMaxi,betai,cPi,dPi;

	for(int i = 0;i < n; i++){
		//
		Ri = R[i];
	  rMaxi = rMax[i];
	  betai = beta[i];
	  cPi = cP[i];
	  dPi = dP[i];

		P[i] = cPi + dPi*exp(-1.0f*pow(rMaxi / Ri, betai));
	}
	return P;

}

//' @title New Holland Wind Profile Time Series
//' @description Wind profile time series at a grid point. Holland et al. 2010.  In this version, the exponent is allowed to vary linearly outside the radius of maximum wind. I.e. rather than take the square root, the exponent varies around 0.5.Currently this version does not have a corresponding vorticity profile set up in wind Vorticity, so it cannot be applied in some wind field modelling.
//' @param f single coriolis parameter at the centre of TC in hz
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param rho density of air in Kg/m3
//' @param R vector of distances from grid points to TC centre in km
//' @param vMax maximum wind velocity calculation in m/s
//' @param beta exponential term for Holland vortex
//' @return array with two columns for velocity and then vorticity.
//' //@example NewHollandWindProfilePi(-1e-4,20,20,1.15,-14,50,1.3)
// [[Rcpp::export]]
NumericMatrix NewHollandWindProfilePi(NumericVector f, NumericVector rMax, NumericVector dP, float rho, NumericVector R, NumericVector vMax, NumericVector beta)
{
	//Holland et al. 2010.  In this version, the exponent is allowed to
	//vary linearly outside the radius of maximum wind.i.e.rather than
	//	take the sqare root, the exponent varies around 0.5.Currently
	//	this version does not have a corresponding vorticity profile set up
	//	in windVorticity, so it cannot be applied in some wind field modelling.
	int n = R.size();
  NumericMatrix VZ(n,2);
	float Ri,fi,rMaxi,vMaxi,dPi;
	float delta, edelta;

	float Bs, deltag, edeltag, rgterm, xn, xx;


	float rGale = 250.0; // Radius for gale force wind. This should be user defined

	for(int i = 0;i < n; i++){
		//
		Ri = R[i];
	  fi = f[i];
	  rMaxi = rMax[i];
	  vMaxi = vMax[i];
	  dPi = dP[i];
	  //TClati = TClat[i];

		Bs = beta[i];//(-0.000044f * powf(dPi / 100.0f, 2.0f) + 0.01 * (dPi / 100.0f) - 0.014f * fabs(TClati) + 1.0);
		deltag = powf(rMaxi / rGale, Bs);
		edeltag = exp(-1.0f * deltag);
		rgterm = Bs * 100.0f * dPi * deltag * edeltag / rho;
		xn = log(17.0f) / log(rgterm);
		xx = 0.5;

		if (Ri > rMaxi)
		{
			xx = (0.5 + (Ri - rMaxi) * (xn - 0.5) / (rGale - rMaxi));
		}

		delta = powf(rMaxi / Ri, Bs);
		edelta = exp(-delta);

		VZ(i,0) = (fi / fabs(fi)) * vMaxi * pow( delta * edelta, xx);
		VZ(i,1) = 0.0f;// Warning dummy value

	}
	return(VZ);
}

//' @title Double Holland Wind Profile Time Series
//' @description Wind profile time series at a grid point. McConochie *et al*'s double Holland vortex model based on Cardone *et al*, 1994.This application is the Coral Sea adaptation of the double vortex model and it can also be used for concentric eye - wall configurations.
//' @param f single coriolis parameter at the centre of TC in hz
//' @param vMax maximum wind velocity calculation in m/s
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param cP TC central pressure in hPa
//' @param rho density of air in Kg/m3
//' @param beta exponential term for Holland vortex
//' @param R vector of distances from grid points to TC centre in km
//' @return array with two columns for velocity and then vorticity.
//' //@example DoubleHollandWindProfilePi(-1e-4,20,20,10,980,1.15,1.2,50)
// [[Rcpp::export]]
NumericMatrix DoubleHollandWindProfilePi(NumericVector f, NumericVector vMax, NumericVector rMax, NumericVector dP, NumericVector cP, float rho, NumericVector beta, NumericVector R)
{
	//McConochie *et al*'s double Holland vortex model (based on Cardone *et
	//al*, 1994).This application is the Coral Sea adaptation of the
	//double vortex model(it can also be used for concentric eye - wall
	//configurations).
	//

	int n = R.size();
  NumericMatrix VZ(n,2);
	float Vi, Ri,fi,dPi,rMax1,vMaxi,cPi,betai;
	float E, d2Vm, aa, bb, cc;

	float cubic = 0.0f; //cubic profile (cubic == 0.1f) to avoid the barotropic instability mentioned in Kepert 2001

	float beta1, beta2;

	float rMax2 = 150.0f;
  //float rMax1 = rMax;
	float gradientV1, gradientV2;

	float chi, psi;

	float dp2,dp1,nu,mu,enu,emu;

	for(int i = 0;i < n; i++){
	//
	Ri = R[i];
	fi = f[i];
	dPi = dP[i]*100.0f; //pa not hPa
	rMax1 = rMax[i];
	vMaxi = vMax[i];
	cPi = cP[i]*100.0f;
	betai = beta[i];

	if (dPi < 1500.0f)
	{
		dp2 = ((dPi / 1500.0f) * (800.0f + (dPi - 800.0f) / 2000.0f));
	}
	else
	{
		dp2 = 800.0f + (dPi - 800.0f) / 2000.0f;
	}

	dp1 = dPi - dp2;


	//Second derivative of the profile
	beta1 = betai;
	beta2 = 7.2f - cPi / 16000.0f;

	beta2 = betai-0.1f; //7.2f - cPi / 16000.0f;

	E = exp(1.0f);

	nu = pow((rMax2 / rMax1), beta2);
  //missing powf on second line
	d2Vm = (-1.0f / (8.0f * pow(4.0f * beta1 * dp1 / (rho * E) + (4.0f * beta2 * dp2 / rho) * nu * exp(-nu) + powf(rMax1 * fi, 2.0f), 1.5f))*
		powf(-(4.0f * (beta1 *beta1) * dp1 / (rho * rMax1 * E)) + (4.0f * (beta1 * beta1) * dp1 / (rho * rMax1 * E)) - (4 * (beta2 *beta2) * dp2 / rho) *
		(nu / rMax1) * exp(-nu) + (4.0f * (beta2 *beta2) * dp2 / rho) *((nu *nu) / rMax1) * exp(-nu) + 2.0f * rMax1 * fi *fi, 2.0f) +
		1.0f / (4.0f * sqrt((4 * beta1 * dp1 / (rho * E)) +
		(4.0f * beta2 * dp2 / rho) * nu * 2.0f +
		exp(-nu) + pow(rMax1 * fi,2.0f)))*
		((4.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))+
		(4.0f * (beta1 *beta1) * dp1 / (rho * (rMax1 *rMax1) * E))-
		(12.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))-
		(4.0f * (beta1 *beta1) * dp1 / (rho * (rMax1 *rMax1) * E))+
		(4.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))+
		(4.0f * (beta2 *beta2*beta2) * dp2 / rho) *
		(nu / (rMax1 *rMax1)) * exp(-nu)+
		(4.0f * (beta2 *beta2) * dp2 / rho) *
		(nu / (rMax1 *rMax1)) * exp(-nu)-
		(12.0f * (beta2 *beta2*beta2) * dp2 / rho) *
		(nu *nu) / (rMax1 *rMax1) * exp(-nu)-
		(4.0f * (beta2 *beta2) * dp2 / rho) *
		(nu *nu) / (rMax1 *rMax1) * exp(-nu)+
		(4.0f * (beta2 *beta2*beta2) * dp2 / rho) *
		(nu *nu*nu) / (rMax1 *rMax1) * exp(-nu)+
		2.0f * fi *fi));



		mu = powf(rMax1 / Ri, beta1);
		nu = powf(rMax2 / Ri, beta2);
		emu = exp(-mu);
		enu = exp(-nu);

		chi = beta1 * dp1 / rho;
		psi = beta2 * dp2 / rho;

		gradientV1 = (chi) * mu * emu;
		gradientV2 = (psi) * nu * enu;


		aa = (d2Vm / 2.0f - (-vMaxi / rMax1) / rMax1) / rMax1;
		bb = (d2Vm - 6.0f * aa * rMax1) / 2.0f;
		cc = -3.0f * aa * rMax1 * rMax1 - 2.0f * bb * rMax1;

		Vi = (fi / fabs(fi) * sqrt(gradientV1 + gradientV2 + (Ri *fi / 2.0f) *(Ri *fi / 2.0f)) - Ri * fabs(fi) / 2.0f);
		//cubic profile to avoid the barotropic instability
		if(cubic == 1.0f){
		if (dPi >= 1500.0f && Ri <= rMax1)
		  {
		  	Vi = (fi / fabs(fi) * Ri * (Ri * (Ri * aa + bb) + cc));
		  }
		}
		VZ(i,0) = Vi;
		VZ(i,1) = 0.0f;
			//(f / fabs(f) * sqrtf(chi * delta * edelta + psi * nu * enu + (f * Ri / 2.0f) *(f * Ri / 2.0f)) / Ri -
			//fabs(f) + (0.5f) *
			//(chi * ddelta * edelta * (1 - delta) +
			//psi * dgamma * egamma * (1 - gamma) +
			//R * self.f ** 2) /
			//np.sqrt(chi * delta * edelta + psi * gamma *
			//egamma + (self.f * R / 2) ** 2))


	}
	return VZ;
}

//' @title Double Holland Pressure Profile Time Series
//' @description Pressure profile time series at a grid point
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param cP TC central pressure in hPa
//' @param beta exponential term for Holland vortex
//' @param R vector of distances from grid points to TC centre in km
//' @return vector of pressures.
//' //@example DoubleHollandPressureProfilePi(20,20,980,1.2,50)
// [[Rcpp::export]]
NumericVector DoubleHollandPressureProfilePi(NumericVector rMax, NumericVector dP, NumericVector cP,  NumericVector beta, NumericVector R)
{
	//Holland pressure profile
	int n = R.size();
  NumericVector P(n);

	float Ri,dPi,rMaxi,cPi;
	float dp1,dp2;
	float beta1, beta2;
	float nu, mu, enu, emu;

	float rMax2 = 150.0f;


	for(int i = 0;i < n; i++){
	dPi = dP[i]*100;
	rMaxi = rMax[i];
	cPi = cP[i]*100;

	if (dPi < 1500.0f)
	{
		dp2 = (dPi / 1500.0f)*(800.0f + (dPi - 800.0f) / 2000.0f);
	}
	else
	{
		dp2 = 800.0f + (dPi - 800.0f) / 2000.0f;
	}

	dp1 = dPi - dp2;

	beta1 = beta[i];
	//beta1 = 7.3f - cP / 16000.0f;
	//beta2 = 7.2f - cPi / 16000.0f;
	beta2 = beta1-0.1f;


		//
		Ri = R[i];
		mu = powf(rMaxi / Ri,beta1);
		nu = powf(rMax2 / Ri,beta2);
		emu = exp(-mu);
		enu = exp(-nu);
		P[i] = (cPi + dp1*emu + dp2*enu)/100.0f;
	}
	return P;
}

//' @title Hubbert Wind Field Time Series
//' @description Time series vortex Wind, wind vectors. Hubbert, G.D., G.J.Holland, L.M.Leslie and M.J.Manton, 1991: A Real - Time System for Forecasting Tropical Cyclone Storm Surges. *Weather and Forecasting*, **6 * *, 86 - 97
//' @param f single coriolis parameter at the centre of TC in hz
//' @param rMax radius of maximum winds in km
//' @param vFm input forward velocity of TC
//' @param thetaFm input forward direction of TC
//' @param Rlam two columns for distances and direction from grid points to TC centre in km
//' @param V velocity profile
//' @param surface equals one if winds are reduced from the gradient level to the surface, otherwise gradient winds.
//' @return array with two columns for zonal and meridional wind speed vector-components.
//' //@example HubbertWindFieldPi(-1e-4,20,2,10,rbind(c(50,35),c(45,40)),c(20,20))
// [[Rcpp::export]]
NumericMatrix HubbertWindFieldPi(NumericVector f, NumericVector rMax, NumericVector vFm, NumericVector thetaFm, NumericMatrix Rlam, NumericVector V,float surface)
{
	//
	//Hubbert, G.D., G.J.Holland, L.M.Leslie and M.J.Manton, 1991:
	//A Real - Time System for Forecasting Tropical Cyclone Storm Surges.
	//	*Weather and Forecasting*, **6 * *, 86 - 97

	int n = V.size();
  NumericMatrix UwVw(n,2);

	float Km = 0.70f; //reduction factor from boundary layer to water surface
  if(surface < 1.0f) Km = 1;
	float inflow;
	float Ri,fi,vFmi,Vi,rMaxi;
	float lami;
	float thetaMax = 70.0f; //70.0f;
	float pi = 3.141592f;
	float piOn180 = pi / 180.0f;
	float thetaFmRAD;
  float sf;
	float thetaMaxAbsolute, asym, Vsf, phi;

	for(int i = 0;i < n; i++){
	  fi = f[i];
	  sf = fi/fabs(fi);
	  vFmi = vFm[i];
	  rMaxi = rMax[i];
	  thetaFmRAD = thetaFm[i]*piOn180;
		//V = self.velocity(R)
		Ri = Rlam(i,0);
		lami = Rlam(i,1)*pi/180.0f;
		Vi = V[i];
		//CB didn't have -1.0f*sign(f)
    inflow = -sf*25.0f;
		if (Ri < rMaxi)
		{
			inflow = 0;
		}

		inflow *= piOn180;
    //CB didn't have * pi / 180.0f
		thetaMaxAbsolute = thetaFmRAD + thetaMax*-sf* piOn180;
		asym = vFmi * cosf(thetaMaxAbsolute - lami + pi);
		Vsf =  Km * (Vi + asym);
		phi = inflow - lami;
    //TCRM has an additional factor 1.069 to convert to 1-munite sustained wind speed
		UwVw(i,0) = Vsf *sinf(phi);
		UwVw(i,1) = Vsf *cosf(phi);

	}
	return UwVw;
}

//' @title McConochie Wind Field Time Series
//' @description Time series vortex Wind, wind vectors. McConochie, J.D., T.A.Hardy and L.B.Mason, 2004: Modelling tropical cyclone over - water wind and pressure fields. Ocean Engineering, 31, 1757 - 1782.
//' @param rMax radius of maximum winds in km
//' @param vMax maximum wind velocity calculation in m/s
//' @param vFm input forward velocity of TC
//' @param thetaFm input forward direction of TC
//' @param Rlam two columns for distances and direction from grid points to TC centre in km
//' @param V velocity profile
//' @param f coriolis parameter at the centre of TC in hz
//' @param surface equals one if winds are reduced from the gradient level to the surface, otherwise gradient winds.
//' @return array with two columns for zonal and meridional wind speed vector-components.
//' //@example McConochieWindFieldPi(-1e-4,20,2,10,rbind(c(50,35),c(45,40)),c(20,20))
// [[Rcpp::export]]
NumericMatrix McConochieWindFieldPi(NumericVector rMax, NumericVector vMax, NumericVector vFm, NumericVector thetaFm, NumericMatrix Rlam, NumericVector V,float f,float surface)
{
	//
	//McConochie, J.D., T.A.Hardy and L.B.Mason, 2004:
	//Modelling tropical cyclone over - water wind and pressure fields.
	//	Ocean Engineering, 31, 1757 - 1782.
	int n = V.size();
  NumericMatrix UwVw(n,2);

	//float Km = 0.70;
	float inflow;
	float Ri,rMaxi,vMaxi,vFmi, Vi;
	float lami;
	float thetaMax = 70.0f;
  float pi = 3.141592f;
	float thetaMaxAbsolute, asym, Vsf, phi,swrf;
	float thetaFmRAD;
  float sf;
  sf = f/fabs(f);
  float piOn180 = pi / 180.0f;

	for(int i = 0;i < n; i++){
		//
		//V = self.velocity(R)
		thetaFmRAD = thetaFm[i]*piOn180;
	  rMaxi = rMax[i];
	  vMaxi = vMax[i];
	  vFmi = vFm[i];
		inflow = 25.0f;
		Ri = Rlam(i,0);
		lami = Rlam(i,1) * piOn180;
		Vi = V[i];
		if (Ri < 1.2f*rMaxi)
		{
			//
			//inflow = 10.0f + 75.0f * (Ri / rMaxi - 1.0f);
			inflow = 75.0f * (Ri / rMaxi ) - 65.0f; //actual eq from McConochie
		}
		if (Ri < rMaxi)
		{
			//
			inflow = 10.0f * Ri / rMaxi;
		}
		inflow = inflow*piOn180;  //missing -sign(f) i.e. f/fabs(f)

		thetaMaxAbsolute = thetaFmRAD + thetaMax*-1.0f*sf*pi/180;
		phi = inflow - lami;

		asym = (0.5f * (1.0f + cosf(thetaMaxAbsolute - lami)) * vFmi * (Vi / vMaxi));
		Vsf = Vi + asym;

		swrf = 0.81f;
		// had an extra ;
		if (abs(Vsf) >= 6.0f)
		{
			swrf = 0.81f - (2.93f * (fabs(Vsf) - 6.0f) / 1000.0f);
		}
		if (abs(Vsf) >= 19.5f)
		{
			swrf = 0.77f - (4.31f * (fabs(Vsf) - 19.5f) / 1000.0f);
		}
		if (abs(Vsf) >= 45.0f)
		{
			swrf = 0.66f;
		}
		if(surface < 1.0f) swrf = 1.0f;
    // TCRM has an additional factor to convert to 1-munite sustained wind speed:
		UwVw(i,0) = swrf * Vsf * sinf(phi);
		UwVw(i,1) = swrf * Vsf * cosf(phi);

	}
	return UwVw;
}

//' @title Kepert Wind Field
//' @description Time series vortex Wind, wind vectors. Kepert, J., 2001: The Dynamics of Boundary Layer Jets within the Tropical Cyclone Core.Part I : Linear Theory.J.Atmos.Sci., 58, 2469 - 2484
//' @param rMax radius of maximum winds in km
//' @param vMax maximum wind velocity calculation in m/s
//' @param vFm input forward velocity of TC
//' @param thetaFm input forward direction of TC
//' @param f single coriolis parameter at the centre of TC in hz
//' @param Rlam two columns for distances and direction from grid points to TC centre in km
//' @param VZ array two columns velocity then vorticity
//' @param surface equals one if winds are reduced from the gradient level to the surface, otherwise gradient winds.
//' @return array with two columns for zonal and meridional wind speed vector-components.
//' //@example KepertWindField(20,20,2,10,-1e-4,rbind(c(50,35),c(45,40)),rbind(c(20,2),c(22,3)))
// [[Rcpp::export]]
NumericMatrix KepertWindFieldPi(NumericVector rMax, NumericVector vMax, NumericVector vFm, NumericVector thetaFm, NumericVector f, NumericMatrix Rlam, NumericMatrix VZ,float surface)
{
	// Kepert, J., 2001: The Dynamics of Boundary Layer Jets within the
	//Tropical Cyclone Core.Part I : Linear Theory.J.Atmos.Sci., 58,
	//	2469 - 2484
	//Orginal code { Written Jeff Kepert, Bureau of Meteorology, 1998-2000.
	//Copyright the Bureau of Meteorology.
	//Please do not distribute orignal work without my knowledge.
	
	//The model is, so far as I know, robust, except if the storm is
	//close to inertially neutral (e.g. b too big). Note that it was written
	//to understand the dynamics, not to make accurate predictions -  
	//the constants (C, K, etc) have not been tuned to observations.
	//Note also that because of a linearisation in the derivation, the
	//model does not produce the correct limit in the limit r -> infinity.
	//This code uses a Holland (1980) parametric profile, but should also
	//work with other reasonable parametric profiles.
	//Created on Fri Oct  2 09:53:29 2015. Port from matlab code.
	//@author: Jeff }
	
	NumericVector V = VZ( _ , 0 );
	int n = V.size();
  NumericMatrix UwVw(n,2);

	float Ri,fi,rMaxi,vMaxi, Vi, Zi,fs,Ks;
	float lami,lami2;
	float K = 50.0f; //diffusivity
	float Cd = 0.002f; // Constant drag coeff


	float Vt,al,be,gam, albe;
	float chi, eta, psi;

	float A0r, A0i,u0s,v0s,Amr,Ami,ums,vms,Apr,Api,ups,vps,vFmi,Umod;
	float us, vs, usf, vsf,phi;
	float pi = 3.141592f;
	float thetaFmi;
  float piOn180 = pi / 180.0f;
  //Vm = fabs(max(V));

	for(int i = 0;i < n; i++){
		lami = Rlam(i,1);
	  lami = lami * piOn180;
	  thetaFmi = thetaFm(i);
		thetaFmi = thetaFmi * piOn180;
	  rMaxi = rMax[i]*1000.0f;
		Ri = Rlam(i,0)*1000.0f;
		fi = f(i);
		fs = fi/fabs(fi);

		Vi = VZ(i,0);
		Zi = VZ(i,1);
    vFmi = vFm(i);
    Umod = vFmi;
    vMaxi = vMax(i);
    if((vFmi > 0) and ((vMaxi/vFmi) < 5.0f)){
      Umod = vFmi*fabs(1.25f * (1.0f - (vFmi / vMaxi)));
    }

		Vt = Umod;
		if (Ri > (2.0f * rMaxi))  //this is not the core https://github.com/GeoscienceAustralia/tcrm/blob/1916233f7dfdecf6a1b5f5b0d89f3eb1d164bd3e/wind/windmodels.py#L1063
		{
			Vt = Umod * expf(-powf((Ri / (2.0f * rMaxi)) - 1.0f, 2.0f));
		}


		al = ((2.0f * Vi / Ri) + fi) / (2.0f * K);
		be = (fi + Zi) / (2.0f * K);
		gam = Vi / (2.0f * K * Ri);


		albe = sqrtf(al / be);

		chi = fabs((Cd / K) * Vi / sqrtf(sqrt(al * be)));
		eta = fabs((Cd / K) * Vi / sqrtf(sqrt(al * be) + fabs(gam)));
		psi = fabs((Cd / K) * Vi / sqrtf(fabs(sqrt(al * be) - fabs(gam))));

		Ks = ((chi * chi) + 2.0f * chi + 2.0f)/(2.0f * chi * chi + 3.0f * chi + 2.0f); //gradient to surface wind reduction factor Eq 30 in Kepert 2001
    if(surface < 1.0f) Ks = 1.0f;
		// converted from complex number formula to this
		A0r = -(chi  * (1.0f + 0.0f * (1.0f + chi)) * Vi) / (2.0f * chi * chi + 3.0f * chi + 2.0f);
		A0i = -(chi  * (1.0f + 1.0f * (1.0f + chi)) * Vi) / (2.0f * chi * chi + 3.0f * chi + 2.0f);

		//Symmetric surface wind component
		u0s = A0r *albe * fs;
		v0s = A0i;

		// converted from complex number formula to this
		Amr = -(psi * (1.0f + 2.0f * albe + (1.0f+0.0f)*(1.0f + albe) * eta) * Vt) / (albe * ((2.0f + 2.0f * 0.0f) * (1.0f + eta * psi) + 3.0f * psi + 3.0f * 0.0f * eta));
		Ami = -(psi * (1.0f + 2.0f * albe + (1.0f+1.0f)*(1.0f + albe) * eta) * Vt) / (albe * ((2.0f + 2.0f * 1.0f) * (1.0f + eta * psi) + 3.0f * psi + 3.0f * 1.0f * eta));

		if (fabs(gam) > sqrtf(al*be))
		{
		  Amr =-(psi * (1.0f + 2.0f * albe + (1.0f+0.0f)*(1.0f + albe) * eta) * Vt) / (albe * ((2.0f - 2.0f * 0.0f + 3.0f * (eta + psi) + (2.0f + 2.0f * 0.0f) * eta * psi)));
		  Ami =-(psi * (1.0f + 2.0f * albe + (1.0f+1.0f)*(1.0f + albe) * eta) * Vt) / (albe * ((2.0f - 2.0f * 1.0f + 3.0f * (eta + psi) + (2.0f + 2.0f * 1.0f) * eta * psi)));

		}

		//First asymmetric surface component
		//ums = (Amr * expf(-0.0f * lami * fs)) * albe;
		//vms = (Ami * expf(-1.0f * lami * fs)) * fs;
		//DEV branch is different to MAIN branch of TCRM for the following lines
		//ums = (Amr * (cosf(-(lami-thetaFmi) * fs) + 0.0f * sinf(-(lami-thetaFmi) * fs))) * albe;  //https://en.wikipedia.org/wiki/Euler%27s_formula
		//vms = (Ami * (cosf(-(lami-thetaFmi) * fs) + 1.0f * sinf(-(lami-thetaFmi) * fs))) * fs;
		lami2 = lami - thetaFmi;
		ums = (Amr * (cosf(lami2 * fs) )) * albe;  //https://en.wikipedia.org/wiki/Euler%27s_formula
		vms = (Ami * (- 1.0f * sinf(lami2 * fs))) * fs;


		Apr = -(eta * (1.0f - 2.0f * albe + (1.0f + 0.0f) * (1.0f - albe) * psi) * Vt) / (albe * ((2.0f + 2.0f * 0.0f) * (1.0f + eta * psi) + 3.0f * eta + 3 * 0.0f * psi));
		Api = -(eta * (1.0f - 2.0f * albe + (1.0f + 1.0f) * (1.0f - albe) * psi) * Vt) / (albe * ((2.0f + 2.0f * 1.0f) * (1.0f + eta * psi) + 3.0f * eta + 3 * 1.0f * psi));

		if (fabs(gam) > sqrtf(al*be))
		{
			Apr = -(eta * (1.0f - 2.0f * albe + (1.0f - 0.0f) * (1.0f - albe) * psi) * Vt) / (albe * (2.0f + 2.0f * 0.0f + 3.0f * (eta + psi) + (2.0f - 2.0f * 0.0f) * eta * psi));
			Api = -(eta * (1.0f - 2.0f * albe + (1.0f - 1.0f) * (1.0f - albe) * psi) * Vt) / (albe * (2.0f + 2.0f * 1.0f + 3.0f * (eta + psi) + (2.0f - 2.0f * 1.0f) * eta * psi));
		}

		//Second asymmetric surface component
		//ups = (Apr * expf(0.0f * lami * fs)) * albe;
		//vps = (Api * expf(1.0f * lami * fs)) * fs;
		//ups = (Apr * (cosf((lami-thetaFmi) * fs) + 0.0f * sinf((lami-thetaFmi) * fs))) * albe;  //https://en.wikipedia.org/wiki/Euler%27s_formula
		//vps = (Api * (cosf((lami-thetaFmi) * fs) + 1.0f * sinf((lami-thetaFmi) * fs))) * albe;

		ups = (Apr * (cosf(lami2 * fs))) * albe;  //https://en.wikipedia.org/wiki/Euler%27s_formula
		vps = (Api * (1.0f * sinf(lami2 * fs))) * fs;//albe;


		//Total surface wind in (moving coordinate system)
		us = u0s + ups + ums;
		vs = v0s + vps + vms + Vi;

		usf = us + Vt * cosf(lami2);
		vsf = vs - Vt * sinf(lami2);
		phi = atan2f(usf, vsf);
	  //if(phi > pi){
	  //  phi = pi - phi;
	  //}
	  float square_term = sqrtf(usf * usf + vsf * vsf);
		UwVw(i,0) = square_term * sinf(phi - lami) * Ks;
		UwVw(i,1) = square_term * cosf(phi - lami) * Ks;

	}
	return UwVw;
}


//' @title TC Distance and Direction From Output Grid Points
//' @description Grid points distance and direction to TC.
//' @param Gridlon vector of Grid point longitudes
//' @param Gridlat vector of Grid point latitudes
//' @param TClon single TC longitude
//' @param TClat single TC latitude
//' @return two columns for distance in km and cartesian direction in degrees, counter clockwise from the x axis.
//' //@example Rdist(c(144,145),c(-11,-12),142,-14)
// [[Rcpp::export]]
NumericMatrix Rdist(NumericVector Gridlon, NumericVector Gridlat, float TClon, float TClat)
{

  //haversine formula
  int n = Gridlon.size();

  NumericMatrix Rlam(n,2);
  //Rlam is a combination of ...
  //R: radius in km from the centre of the TC to the gird point
  //lam: Direction(geographic bearing, positive clockwise) !!!isn't this the Cartesian direction, counter clockwise from the x axis!!!
  float Rearth = 6372797.560856f;
  float pi = 3.141592f;
  float piOn180 = pi / 180.0f;

  float lat1 = TClat * piOn180;
  float lon1 = TClon * piOn180;

  float cos_lat1 = cosf(lat1);
  float sin_lat1 = sinf(lat1);


  for(int i = 0;i < n; i++){
    float lat2 = Gridlat[i] * piOn180;
    float lon2 = Gridlon[i] * piOn180;

    float dlat = lat2 - lat1;
    float dlon = lon2 - lon1;

    float cos_lat2 = cosf(lat2);

    float a = sinf(dlat / 2.0f) * sinf(dlat / 2.0f) + cos_lat1 * cos_lat2 * sinf(dlon / 2.0f) * sinf(dlon / 2.0f);
    float c = 2.0f * atan2f(sqrtf(a), sqrtf(1.0f - a));
    Rlam(i,0) = c * Rearth/1000.0f;//convert to km

    float x = sinf(dlon) * cos_lat2;
    float y = cos_lat1 * sinf(lat2) - sin_lat1 * cos_lat2 * cosf(dlon);
    Rlam(i,1) = atan2f(y, x) / piOn180;
  }
  return Rlam;
}

//' @title Jelesnianski Wind Profile
//' @description wind profile at grid points
//' @param f single coriolis parameter at the centre of TC in hz
//' @param vMax maximum wind velocity calculation in m/s
//' @param rMax radius of maximum winds in km
//' @param R vector of distances from grid points to TC centre in km
//' @return array with two columns for velocity and then vorticity.
//' //@example JelesnianskiWindProfile(-1e-4,20,20,50)
// [[Rcpp::export]]// [[Rcpp::export]]
NumericMatrix JelesnianskiWindProfile(float f, float vMax, float rMax, NumericVector R)
{
  //
  int n = R.size();
  NumericMatrix VZ(n,2);

  float Vi, Ri,Zi, sf;

  sf = (f / fabs(f));
  for(int i = 0;i < n; i++){
    //
    Ri = R[i];
    Vi = 2.0f * vMax * rMax * Ri / (rMax *rMax + Ri * Ri) * sf;
    VZ(i,0) = Vi;
    Zi = (sf * 2.0f * vMax * rMax / (rMax *rMax + Ri * Ri) + sf * 2.0f * vMax * rMax * (rMax *rMax - Ri * Ri) /((rMax *rMax + Ri * Ri) * (rMax *rMax + Ri * Ri)));
    VZ(i,1) = Zi;
  }
  return VZ;
}

//' @title Holland Wind Profile
//' @description wind profile at grid points
//' @param f single coriolis parameter at the centre of TC in hz
//' @param vMax maximum wind velocity calculation in m/s
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param rho density of air in Kg/m3
//' @param beta exponential term for Holland vortex
//' @param R vector of distances from grid points to TC centre in km
//' @return array with two columns for velocity and then vorticity.
//' //@example HollandWindProfile(-1e-4,20,20,10,1.15,1.2,50)
// [[Rcpp::export]]
NumericMatrix HollandWindProfile(float f, float vMax, float rMax, float dP, float rho, float beta, NumericVector R)
{
  //Holland profile. For `r < rMax`, we reset the wind feld to a
  //cubic profile to avoid the barotropic instability mentioned in
  //Kepert & Wang (2001).

  int n = R.size();
  NumericMatrix VZ(n,2);
  float Vi, Ri, Zi,fs;
  float E, d2Vm,dVm,aa,bb,cc;
  float delta, edelta;
  fs = f / fabs(f);
  dP *= 100.0f;
  E = expf(1.0f);

  for(int i = 0;i < n; i++){
    //
    Ri = R[i];


    if (Ri <= rMax){

      d2Vm = ((beta * dP * (-4.0f * beta *beta *beta * dP / rho - (-2.0f + beta *beta) * E * (f * rMax) *(f * rMax))) / (E * rho * sqrtf((4.0f * beta * dP) / (E * rho) + (f * rMax) *(f * rMax)) * (4.0f * beta * dP * rMax *rMax / rho + E * (f * rMax *rMax) *(f * rMax *rMax))));

      dVm = (-fabs(f)/2.0f + (E * (f*f) * rMax * sqrtf((4.0f * beta * dP / rho) / E + (f * rMax) * (f * rMax))) / (2.0f * (4.0f * beta * dP / rho + E * (f * rMax) * (f * rMax))));

      aa = ((d2Vm / 2.0f - (dVm - vMax / rMax) / rMax) / rMax);

      bb = (d2Vm - 6.0f * aa * rMax) / 2.0f;

      cc = dVm -3.0f * aa * rMax * rMax - 2.0f * bb * rMax;

      Vi = (Ri * (Ri * (Ri * aa + bb) + cc));
      Zi = Ri * (Ri * 4.0f * aa + 3.0f * bb) + 2.0f * cc;
    }
    else
    {
      delta = powf(rMax / Ri, beta);
      edelta = expf(-delta);
      //dP = dP*100.0f; // !need to multiply by dP 100 to get to Pa not in TCRM  https://github.com/GeoscienceAustralia/tcrm/blob/1916233f7dfdecf6a1b5f5b0d89f3eb1d164bd3e/wind/windmodels.py#L381
      Vi =  (sqrtf((dP * beta / rho) * delta * edelta + (Ri * f / 2.0f)*(Ri * f / 2.0f)) - Ri *fabs(f) / 2.0f);
      //Zi = ((sqrtf((dP * beta / rho) * delta * edelta + (Ri * f / 2.0f)*(Ri * f / 2.0f))) / Ri - fabs(f) + edelta * (2.0f * (beta * beta) * dP * (delta - 1.0f) * delta + rho * edelta * (f * Ri) *(f * Ri)) / (2.0f * rho * Ri * sqrtf(4.0f * (beta * dP / rho) * delta * edelta + (f * Ri) *(f * Ri))));
      Zi = fabs(f) + (beta*beta * dP * (delta * delta) * edelta / (2.0f * rho * Ri) - beta*beta * dP * delta * edelta / (2.0f * rho * Ri) + Ri * f * f / 4.0f) / sqrtf(beta * dP * delta * edelta / rho + (Ri * f / 2)*(Ri * f / 2)) + (sqrtf(beta * dP * delta * edelta / rho + (Ri * f / 2)*(Ri * f / 2))) / Ri;
    }
    VZ(i,0) = Vi * fs;
    VZ(i,1) = Zi * fs;
  }
  return VZ;
}

//' @title Holland Pressure Profile
//' @description Pressure profile at grid points
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param cP TC central pressure in hPa
//' @param beta exponential term for Holland vortex
//' @param R vector of distances from grid points to TC centre in km
//' @return vector of pressures.
//' //@example HollandPressureProfile(20,20,980,1.2,50)
// [[Rcpp::export]]
NumericVector HollandPressureProfile(float rMax, float dP, float cP, float beta, NumericVector R)
{
  //Holland pressure profile
  int n = R.size();
  NumericVector P(n);


  float Ri;

  for(int i = 0;i < n; i++){
    //
    Ri = R[i];
    P[i] = cP + dP*exp(-1.0f*pow(rMax / Ri, beta));
  }
  return P;

}

//' @title New Holland Wind Profile Time Series
//' @description Wind profile time series at a grid point. Holland et al. 2010.  In this version, the exponent is allowed to vary linearly outside the radius of maximum wind. I.e. rather than take the square root, the exponent varies around 0.5.Currently this version does not have a corresponding vorticity profile set up in wind Vorticity, so it cannot be applied in some wind field modelling.
//' @param f single coriolis parameter at the centre of TC in hz
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param rho density of air in Kg/m3
//' @param R vector of distances from grid points to TC centre in km
//' @param vMax maximum wind velocity calculation in m/s
//' @param beta exponential term for Holland vortex
//' @return array with two columns for velocity and then vorticity.
//' //@example NewHollandWindProfile(-1e-4,20,20,1.15,-14,50,1.3)
// [[Rcpp::export]]
NumericMatrix NewHollandWindProfile(float f, float rMax, float dP, float rho, NumericVector R, float vMax, float beta)
{
  //Holland et al. 2010.  In this version, the exponent is allowed to
  //vary linearly outside the radius of maximum wind.i.e.rather than
  //	take the sqare root, the exponent varies around 0.5.Currently
  //	this version does not have a corresponding vorticity profile set up
  //	in windVorticity, so it cannot be applied in some wind field modelling.
  int n = R.size();
  NumericMatrix VZ(n,2);
  float Ri;
  float delta, edelta;

  float Bs, deltag, edeltag, rgterm, xn, xx;
  float sf;
  sf = (f / fabs(f));
  Bs = beta;
  float rGale = 250.0; // Radius for gale force wind. This should be user defined
  for(int i = 0;i < n; i++){
    //
    Ri = R[i];
    //(-0.000044f * powf(dPi / 100.0f, 2.0f) + 0.01 * (dPi / 100.0f) - 0.014f * fabs(TClati) + 1.0);
    deltag = powf(rMax / rGale, Bs);
    edeltag = exp(-1.0f * deltag);
    rgterm = Bs * 100.0f * dP * deltag * edeltag / rho;
    xn = log(17.0f) / log(rgterm);
    xx = 0.5;

    if (Ri > rMax)
    {
      xx = (0.5 + (Ri - rMax) * (xn - 0.5) / (rGale - rMax));
    }

    delta = powf(rMax / Ri, Bs);
    edelta = exp(-delta);

    VZ(i,0) = sf * vMax*pow(delta * edelta, xx);
    VZ(i,1) = 0.0f;// Warning dummy value

  }
  return(VZ);
}


//' @title Double Holland Wind Profile
//' @description McConochie *et al*'s double Holland vortex model based on Cardone *et al*, 1994.This application is the Coral Sea adaptation of the double vortex model and it can also be used for concentric eye - wall configurations.
//' @param f single coriolis parameter at the centre of TC in hz
//' @param vMax maximum wind velocity calculation in m/s
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param cP TC central pressure in hPa
//' @param rho density of air in Kg/m3
//' @param beta exponential term for Holland vortex
//' @param R vector of distances from grid points to TC centre in km
//' @return array with two columns for velocity and then vorticity.
//' //@example DoubleHollandWindProfile(-1e-4,20,20,10,980,1.15,1.2,50)
// [[Rcpp::export]]
NumericMatrix DoubleHollandWindProfile(float f, float vMax, float rMax, float dP, float cP, float rho, float beta, NumericVector R)
{
  //McConochie *et al*'s double Holland vortex model (based on Cardone *et
  //al*, 1994).This application is the Coral Sea adaptation of the
  //double vortex model(it can also be used for concentric eye - wall
  //confgurations).
  //

  int n = R.size();
  NumericMatrix VZ(n,2);
  float Vi, Ri,rMax1;
  float E, d2Vm, aa, bb, cc;
  float cubic = 0.0f; //cubic profile (cubic == 0.1f) to avoid the barotropic instability mentioned in Kepert 2001

  float beta1, beta2;

  float rMax2 = 150.0f;
  //float rMax1 = rMax;
  float gradientV1, gradientV2;

  float chi, psi;

  float dp2,dp1,nu,mu,enu,emu;
  dP = dP*100.0f; //pa not hPa
  cP = cP*100.0f;
  rMax1 = rMax;

    if (dP < 1500.0f)
    {
      dp2 = ((dP / 1500.0f) * (800.0f + (dP - 800.0f) / 2000.0f));
    }
    else
    {
      dp2 = 800.0f + (dP - 800.0f) / 2000.0f;
    }

    dp1 = dP - dp2;


    //Second derivative of the Profile
    beta1 = beta;
    //beta2 = 7.2f - cP / 16000.0f;

    //beta2 = 7.2f - cP / 16000.0f;
    beta2 = beta - 0.1f;
    E = exp(1.0f);
    if(cubic == 1.0f){
    nu = pow((rMax2 / rMax1), beta2);
    //missing powf on second line
    d2Vm = (-1.0f / (8.0f * pow(4.0f * beta1 * dp1 / (rho * E) + (4.0f * beta2 * dp2 / rho) * nu * exp(-nu) + powf(rMax1 * f, 2.0f), 1.5f))*
      powf(-(4.0f * (beta1 *beta1) * dp1 / (rho * rMax1 * E)) + (4.0f * (beta1 * beta1) * dp1 / (rho * rMax1 * E)) - (4 * (beta2 *beta2) * dp2 / rho) *
      (nu / rMax1) * exp(-nu) + (4.0f * (beta2 *beta2) * dp2 / rho) *((nu *nu) / rMax1) * exp(-nu) + 2.0f * rMax1 * f *f, 2.0f) +
      1.0f / (4.0f * sqrt((4 * beta1 * dp1 / (rho * E)) +
      (4.0f * beta2 * dp2 / rho) * nu * 2.0f +
      exp(-nu) + pow(rMax1 * f,2.0f)))*
      ((4.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))+
      (4.0f * (beta1 *beta1) * dp1 / (rho * (rMax1 *rMax1) * E))-
      (12.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))-
      (4.0f * (beta1 *beta1) * dp1 / (rho * (rMax1 *rMax1) * E))+
      (4.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))+
      (4.0f * (beta2 *beta2*beta2) * dp2 / rho) *
      (nu / (rMax1 *rMax1)) * exp(-nu)+
      (4.0f * (beta2 *beta2) * dp2 / rho) *
      (nu / (rMax1 *rMax1)) * exp(-nu)-
      (12.0f * (beta2 *beta2*beta2) * dp2 / rho) *
      (nu *nu) / (rMax1 *rMax1) * exp(-nu)-
      (4.0f * (beta2 *beta2) * dp2 / rho) *
      (nu *nu) / (rMax1 *rMax1) * exp(-nu)+
      (4.0f * (beta2 *beta2*beta2) * dp2 / rho) *
      (nu *nu*nu) / (rMax1 *rMax1) * exp(-nu)+
      2.0f * f *f));

      aa = (d2Vm / 2.0f - (-vMax / rMax1) / rMax1) / rMax1;
      bb = (d2Vm - 6.0f * aa * rMax1) / 2.0f;
      cc = -3.0f * aa * rMax1 * rMax1 - 2.0f * bb * rMax1;
    }
  for(int i = 0;i < n; i++){
    //
    Ri = R[i];

    mu = powf(rMax1 / Ri, beta1);
    nu = powf(rMax2 / Ri, beta2);
    emu = exp(-mu);
    enu = exp(-nu);

    chi = beta1 * dp1 / rho;
    psi = beta2 * dp2 / rho;

    gradientV1 = (chi) * mu * emu;
    gradientV2 = (psi) * nu * enu;

    Vi = (f / fabs(f) * sqrt(gradientV1 + gradientV2 + (Ri *f / 2.0f) *(Ri *f / 2.0f)) - Ri * fabs(f) / 2.0f);

    //cubic profile to avoid the barotropic instability
    if(cubic == 0.1f){
    if (dP >= 1500.0f && Ri <= rMax1){
    {
      Vi = (f / fabs(f) * Ri * (Ri * (Ri * aa + bb) + cc));
    }
    }
    }

    VZ(i,0) = Vi;
    VZ(i,1) = 0.0f;

    //(f / fabs(f) * sqrtf(chi * delta * edelta + psi * nu * enu + (f * Ri / 2.0f) *(f * Ri / 2.0f)) / Ri -
    //fabs(f) + (0.5f) *
    //(chi * ddelta * edelta * (1 - delta) +
    //psi * dgamma * egamma * (1 - gamma) +
    //R * self.f ** 2) /
    //np.sqrt(chi * delta * edelta + psi * gamma *
    //egamma + (self.f * R / 2) ** 2))


  }
  return VZ;
}

//' @title Double Holland Pressure Profile
//' @description Pressure profile at grid points
//' @param rMax radius of maximum winds in km
//' @param dP pressure differential, environmental less TC central pressure in hPa
//' @param cP TC central pressure in hPa
//' @param beta exponential term for Holland vortex
//' @param R vector of distances from grid points to TC centre in km
//' @return vector of pressures.
//' //@example DoubleHollandPressureProfile(20,20,980,1.2,50)
// [[Rcpp::export]]
NumericVector DoubleHollandPressureProfile(float rMax, float dP, float cP,  float beta, NumericVector R)
{
  //Holland pressure Profile
  int n = R.size();
  NumericVector P(n);

  float Ri;
  float dp1,dp2;
  float beta1, beta2;
  float nu, mu, enu, emu;

  float rMax2 = 150.0f;
  float rMax1 = rMax;
  dP = dP*100;
  cP = cP*100;
  for(int i = 0;i < n; i++){

    if (dP < 1500.0f)
    {
      dp2 = (dP / 1500.0f)*(800.0f + (dP - 800.0f) / 2000.0f);
    }
    else
    {
      dp2 = 800.0f + (dP - 800.0f) / 2000.0f;
    }

    dp1 = dP - dp2;

    beta1 = beta;
    //beta1 = 7.3f - cP / 16000.0f;
    //beta2 = 7.2f - cP / 16000.0f;
    beta2 = beta - 0.1f;

    //
    Ri = R[i];
    mu = powf(rMax1 / Ri,beta1);
    nu = powf(rMax2 / Ri,beta2);
    emu = exp(-mu);
    enu = exp(-nu);
    P[i] = (cP + dp1*emu + dp2*enu)/100.0f;
  }
  return P;
}

//' @title Hubbert Wind Field
//' @description Grid point vortex Wind field, wind vectors. Hubbert, G.D., G.J.Holland, L.M.Leslie and M.J.Manton, 1991: A Real - Time System for Forecasting Tropical Cyclone Storm Surges. *Weather and Forecasting*, **6 * *, 86 - 97
//' @param f single coriolis parameter at the centre of TC in hz
//' @param rMax radius of maximum winds in km
//' @param vFm input forward velocity of TC
//' @param thetaFm input forward direction of TC
//' @param Rlam two columns for distances and direction from grid points to TC centre in km
//' @param V velocity profile
//' @param surface equals one if winds are reduced from the gradient level to the surface, otherwise gradient winds.
//' 
//' @return array with two columns for zonal and meridional wind speed vector-components.
//' //@example HubbertWindField(-1e-4,20,2,10,rbind(c(50,35),c(45,40)),c(20,20))
// [[Rcpp::export]]
NumericMatrix HubbertWindField(float f, float rMax, float vFm, float thetaFm, NumericMatrix Rlam, NumericVector V,float surface)
{
  //
  //Hubbert, G.D., G.J.Holland, L.M.Leslie and M.J.Manton, 1991:
  //A Real - Time System for Forecasting Tropical Cyclone Storm Surges.
  //	*Weather and Forecasting*, **6 * *, 86 - 97

  int n = V.size();
  NumericMatrix UwVw(n,2);
  float Km = 0.70f; // gradient to surface wind reduction factor
  if(surface < 1.0f) Km = 1.0f; //no reduction

  float inflow;
  float Ri,Vi;
  float lami;
  float thetaMax = 70.0f; //70.0f;
  float pi = 3.141592f;
  float thetaFmRAD;
  float sf;
  float thetaMaxAbsolute, asym, Vsf, phi;
  sf = f/fabs(f);
  if(sf > 0) thetaMax = thetaMax+180.0f; //correct for northern hemisphere
  float piOn180 = pi / 180.0f;

  thetaFmRAD = thetaFm * piOn180;
  thetaMaxAbsolute = thetaFmRAD + thetaMax*-1.0 * sf * piOn180;

  for(int i = 0;i < n; i++){

    //V = self.velocity(R)
    Ri = Rlam(i,0);
    lami = Rlam(i,1) * piOn180;
    Vi = V[i];
    //CB didn't have -1.0f*sign(f)
    inflow = -1.0f*sf*25.0f;
    if (Ri < rMax)
    {
      inflow = 0;
    }

    inflow = inflow * piOn180;
    //CB didn't have * pi / 180.0f

    asym = vFm * cosf(thetaMaxAbsolute - lami + pi);
    Vsf =  Km * (Vi + asym);
    phi = inflow - lami;
    //TCRM has an additional factor 1.069 to convert to 1-munite sustained wind speed
    UwVw(i,0) = Vsf * sinf(phi);
    UwVw(i,1) = Vsf * cosf(phi);

  }
  return UwVw;
}

//' @title McConochie Wind Field
//' @description Grid point vortex Wind field, wind vectors. McConochie, J.D., T.A.Hardy and L.B.Mason, 2004: Modelling tropical cyclone over - water wind and pressure fields. Ocean Engineering, 31, 1757 - 1782.
//' @param rMax radius of maximum winds in km
//' @param vMax maximum wind velocity calculation in m/s
//' @param vFm input forward velocity of TC
//' @param thetaFm input forward direction of TC
//' @param Rlam two columns for distances and direction from grid points to TC centre in km
//' @param V velocity profile
//' @param f coriolis parameter at the centre of TC in hz
//' @param surface equals one if winds are reduced from the gradient level to the surface, otherwise gradient winds.
//' @return array with two columns for zonal and meridional wind speed vector-components.
//' //@example McConochieWindField(-1e-4,20,2,10,rbind(c(50,35),c(45,40)),c(20,20))
// [[Rcpp::export]]
NumericMatrix McConochieWindField(float rMax, float vMax, float vFm, float thetaFm, NumericMatrix Rlam, NumericVector V,float f,float surface)
{
  //
  //McConochie, J.D., T.A.Hardy and L.B.Mason, 2004:
  //Modelling tropical cyclone over - water wind and pressure fields.
  //	Ocean Engineering, 31, 1757 - 1782.
  int n = V.size();
  NumericMatrix UwVw(n,2);

  //float Km = 0.70;
  float inflow;
  float Ri, Vi;
  float lami;
  float thetaMax = 70.0f;
  float pi = 3.141592f;
  float thetaMaxAbsolute, asym, Vsf, phi,swrf;
  float thetaFmRAD;
  float piOn180 = pi / 180.0f;
  thetaFmRAD = thetaFm * piOn180;
  rMax = rMax;
  float sf;
  sf = f/fabs(f);
  thetaMaxAbsolute = thetaFmRAD + thetaMax*-1.0f*sf;
  for(int i = 0;i < n; i++){
    //
    //V = self.velocity(R)
    Ri = Rlam(i,0);
    lami = Rlam(i,1) * piOn180;
    Vi = V[i];
    inflow = 25.0f;
    if (Ri < 1.2f*rMax)
    {
      // !!!TCRM = 35 at Ri == rMax, it should equal 25
      //inflow = 10.0f + 75.0f * (Ri / rMax - 1.0f);
      inflow = 75.0f * (Ri / rMax ) - 65.0f; //actual eq from McConochie
    }
    if (Ri < rMax)
    {
      //
      inflow = 10.0f * Ri / rMax;
    }
    inflow = inflow * piOn180;  //missing -sign(f) i.e. f/fabs(f)


    phi = inflow - lami;

    asym = (0.5f * (1.0f + cosf(thetaMaxAbsolute - lami)) * vFm * (Vi / vMax));
    Vsf = Vi + asym;

    //gradient to surface wind reduction factor
    swrf = 0.81f;
    // had an extra ;
    if (fabs(Vsf) >= 6.0f) swrf = 0.81f - (2.93f * (fabs(Vsf) - 6.0f) / 1000.0f);
    if (fabs(Vsf) >= 19.5f) swrf = 0.77f - (4.31f * (fabs(Vsf) - 19.5f) / 1000.0f);
    if (fabs(Vsf) >= 45.0f) swrf = 0.66f;
    if(surface < 1.0f) swrf = 1.0f;  //no reduction
    // TCRM has an additional factor to convert to 1-munite sustained wind speed:
    UwVw(i,0) = swrf * Vsf * sinf(phi);
    UwVw(i,1) = swrf * Vsf * cosf(phi);

  }
  return UwVw;
}

//' @title Kepert Wind Field
//' @description Grid point vortex Wind field, wind vectors. Kepert, J., 2001: The Dynamics of Boundary Layer Jets within the Tropical Cyclone Core.Part I : Linear Theory.J.Atmos.Sci., 58, 2469 - 2484
//' @param rMax radius of maximum winds in km
//' @param vMax maximum wind velocity calculation in m/s
//' @param vFm input forward velocity of TC
//' @param thetaFm input forward direction of TC
//' @param f single coriolis parameter at the centre of TC in hz
//' @param Rlam two columns for distances and Cartesian direction clocwise from the x axis from grid points to TC centre in km
//' @param VZ array two columns velocity then vorticity
//' @param surface equals one if winds are reduced from the gradient level to the surface, otherwise gradient winds.
//' @return array with two columns for zonal and meridional wind speed vector-components.
//' //@example KepertWindField(20,20,2,10,-1e-4,rbind(c(50,35),c(45,40)),rbind(c(20,2),c(22,3)))
// [[Rcpp::export]]
NumericMatrix KepertWindField(float rMax, float vMax, float vFm, float thetaFm, float f, NumericMatrix Rlam, NumericMatrix VZ,float surface)
{
  // Kepert, J., 2001: The Dynamics of Boundary Layer Jets within the
  //Tropical Cyclone Core.Part I : Linear Theory.J.Atmos.Sci., 58,
  //	2469 - 2484
  //Orginal code { Written Jeff Kepert, Bureau of Meteorology, 1998-2000.
  //Copyright the Bureau of Meteorology.
  //Please do not distribute orignal work without my knowledge.
  
  //The model is, so far as I know, robust, except if the storm is
  //close to inertially neutral (e.g. b too big). Note that it was written
  //to understand the dynamics, not to make accurate predictions -  
  //the constants (C, K, etc) have not been tuned to observations.
  //Note also that because of a linearisation in the derivation, the
  //model does not produce the correct limit in the limit r -> infinity.
  //This code uses a Holland (1980) parametric profile, but should also
  //work with other reasonable parametric profiles.
  //Created on Fri Oct  2 09:53:29 2015. Port from matlab code.
  //@author: Jeff }
  
  NumericVector V = VZ( _ , 0 );
  int n = V.size();
  NumericMatrix UwVw(n,2);

  float Ri, Vi, Zi,fs;
  float lami;
  float K = 50.0f; //diffusivity
  float Cd = 0.002f; // Constant drag coeff

  float Vt,al,be,gam, albe;
  float chi, eta, psi;

  float A0r, A0i,u0s,v0s,Amr,Ami,ums,vms,Apr,Api,ups,vps,Umod;

  float us, vs, usf, vsf,phi;
  float pi = 3.141592f;
  float piOn180 = pi / 180.0f;
  rMax  = rMax*1000.0f;
  float lami2,Ks;

  fs = f/fabs(f);
  Umod = vFm;
  thetaFm = thetaFm * piOn180;
  if((vFm > 0) and ((vMax/vFm) < 5.0f)){
    Umod = vFm*fabs(1.25f * (1.0f - (vFm / vMax)));
  }


  for(int i = 0;i < n; i++){
    lami = Rlam(i,1);
    lami = lami * piOn180;

    Ri = Rlam(i,0)*1000.0f;

    Vi = VZ(i,0);
    Zi = VZ(i,1);

    //if(Vi*fs < 0) Vi = Vi*fs; //Vi and f must have the same sign
    //if(Zi*fs < 0) Zi = Zi*fs;

    Vt = Umod;
    if (Ri >= (2.0f * rMax))  //this is not the core https://github.com/GeoscienceAustralia/tcrm/blob/1916233f7dfdecf6a1b5f5b0d89f3eb1d164bd3e/wind/windmodels.py#L1063
    {
      Vt = Umod * expf(-powf((Ri / (2.0f * rMax)) - 1.0f, 2.0f));
    }


    al = ((2.0f * Vi / Ri) + f) / (2.0f * K);
    be = (f + Zi) / (2.0f * K);
    gam = Vi / (2.0f * K * Ri);


    albe = sqrtf(al / be);

    chi = fabs((Cd / K) * Vi / sqrtf(sqrt(al * be)));
    eta = fabs((Cd / K) * Vi / sqrtf(sqrt(al * be) + fabs(gam)));
    psi = fabs((Cd / K) * Vi / sqrtf(fabs(sqrt(al * be) - fabs(gam))));

    Ks = ((chi * chi) + 2.0f * chi + 2.0f)/(2.0f * chi * chi + 3.0f * chi + 2.0f); //gradient to surface wind reduction factor Eq 30 in Kepert 2001
    if(surface < 1.0f) Ks = 1.0f;

    // converted from complex number formula to this
    A0r = -(chi  * (1.0f + 0.0f * (1.0f + chi)) * Vi) / (2.0f * chi * chi + 3.0f * chi + 2.0f);
    A0i = -(chi  * (1.0f + 1.0f * (1.0f + chi)) * Vi) / (2.0f * chi * chi + 3.0f * chi + 2.0f);

    //Symmetric surface wind component
    u0s = A0r *albe * fs;
    v0s = A0i;


    // converted from complex number formula to this
    Amr = -(psi * (1.0f + 2.0f * albe + (1.0f+0.0f)*(1.0f + albe) * eta) * Vt) / (albe * ((2.0f + 2.0f * 0.0f) * (1.0f + eta * psi) + 3.0f * psi + 3.0f * 0.0f * eta));
    Ami = -(psi * (1.0f + 2.0f * albe + (1.0f+1.0f)*(1.0f + albe) * eta) * Vt) / (albe * ((2.0f + 2.0f * 1.0f) * (1.0f + eta * psi) + 3.0f * psi + 3.0f * 1.0f * eta));

    if (fabs(gam) > sqrtf(al*be))
    {
      Amr =-(psi * (1.0f + 2.0f * albe + (1.0f+0.0f)*(1.0f + albe) * eta) * Vt) / (albe * ((2.0f - 2.0f * 0.0f + 3.0f * (eta + psi) + (2.0f + 2.0f * 0.0f) * eta * psi)));
      Ami =-(psi * (1.0f + 2.0f * albe + (1.0f+1.0f)*(1.0f + albe) * eta) * Vt) / (albe * ((2.0f - 2.0f * 1.0f + 3.0f * (eta + psi) + (2.0f + 2.0f * 1.0f) * eta * psi)));

    }

    //frst asymmetric surface component
    //ums = (Amr * expf(-0.0f * lami * fs)) * albe;
    //vms = (Ami * expf(-1.0f * lami * fs)) * fs;
    lami2 = lami-thetaFm;//- (pi/2.0f-thetaFm+pi);
    //lami2 = lami-(pi/2.0f-thetaFm);//- (pi/2.0f-thetaFm+pi);

    //ums = (Amr * (cosf(lami2 * fs) - 0.0f * sinf(lami2 * fs))) * albe;  //https://en.wikipedia.org/wiki/Euler%27s_formula
    //vms = (Ami * (cosf(lami2 * fs) - 1.0f * sinf(lami2 * fs))) * fs;

    ums = (Amr * (cosf(lami2 * fs) )) * albe;  //https://en.wikipedia.org/wiki/Euler%27s_formula
    vms = (Ami * (- 1.0f * sinf(lami2 * fs))) * fs;

    Apr = -(eta * (1.0f - 2.0f * albe + (1.0f + 0.0f) * (1.0f - albe) * psi) * Vt) / (albe * ((2.0f + 2.0f * 0.0f) * (1.0f + eta * psi) + 3.0f * eta + 3.0f * 0.0f * psi));
    Api = -(eta * (1.0f - 2.0f * albe + (1.0f + 1.0f) * (1.0f - albe) * psi) * Vt) / (albe * ((2.0f + 2.0f * 1.0f) * (1.0f + eta * psi) + 3.0f * eta + 3.0f * 1.0f * psi));

    if (fabs(gam) > sqrtf(al*be))
    {
      Apr = -(eta * (1.0f - 2.0f * albe + (1.0f - 0.0f) * (1.0f - albe) * psi) * Vt) / (albe * (2.0f + 2.0f * 0.0f + 3.0f * (eta + psi) + (2.0f - 2.0f * 0.0f) * eta * psi));
      Api = -(eta * (1.0f - 2.0f * albe + (1.0f - 1.0f) * (1.0f - albe) * psi) * Vt) / (albe * (2.0f + 2.0f * 1.0f + 3.0f * (eta + psi) + (2.0f - 2.0f * 1.0f) * eta * psi));
    }

    //Second asymmetric surface component
    //ups = (Apr * expf(0.0f * lami * fs)) * albe;
    //vps = (Api * expf(1.0f * lami * fs)) * fs;
    //ups = (Apr * (cosf(lami2 * fs) + 0.0f * sinf(lami2 * fs))) * albe;  //https://en.wikipedia.org/wiki/Euler%27s_formula
    //vps = (Api * (cosf(lami2 * fs) + 1.0f * sinf(lami2 * fs))) * fs;//albe;

    ups = (Apr * (cosf(lami2 * fs))) * albe;  //https://en.wikipedia.org/wiki/Euler%27s_formula
    vps = (Api * (1.0f * sinf(lami2 * fs))) * fs;//albe;

    //Total surface wind in (moving coordinate system)
    us = u0s + ups + ums;
    vs = v0s + vps + vms + Vi;

    usf = us + Vt * cosf(lami2);
    vsf = vs - Vt * sinf(lami2);
    phi = atan2f(usf, vsf);
    //if(phi > pi){
    //  phi = pi - phi;
    //}
    float square_term = sqrtf(usf*usf + vsf *vsf);
    UwVw(i,0) = square_term * sinf(phi - lami)*Ks;
    UwVw(i,1) = square_term * cosf(phi - lami)*Ks;

  }
  return UwVw;
}
