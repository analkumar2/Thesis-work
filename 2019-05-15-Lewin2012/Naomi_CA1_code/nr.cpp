/*
 *  nr.cpp
 *  
 *
 *  Created by Naomi Lewin on 1/16/09--funtions from Jun's ca3_hippocampus.cpp
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include "nr.h"


void NR::sprsin(Mat_I_DP &a, const DP thresh, Vec_O_DP &sa, Vec_O_INT &ija)  // converts matrix a into a row-indexed sparse storage mode in vectors sa and ija

{

	int i,j,k;



	int n=a.nrows();

	int nmax=sa.size();

	for (j=0;j<n;j++) sa[j]=a[j][j];

	ija[0]=n+1;

	k=n;

	for (i=0;i<n;i++) {

		for (j=0;j<n;j++) {

			if (fabs(a[i][j]) >= thresh && i != j) {

				if (++k > nmax) nrerror("sprsin: sa and ija too small");

				sa[k]=a[i][j];

				ija[k]=j;

			}

		}

		ija[i+1]=k+1;

	}

}

void NR::linbcg(Vec_I_DP &b, Vec_IO_DP &x, const int itol, const DP tol,

	const int itmax, int &iter, DP &err, Vec_INT * ija_p, Vec_DP * sa_p)

{

	DP ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;

	const DP EPS=1.0e-14;

	int j;



	int n=b.size();

	Vec_DP p(n),pp(n),r(n),rr(n),z(n),zz(n);

	iter=0;

	atimes(x,r,0, ija_p, sa_p);

	for (j=0;j<n;j++) {

		r[j]=b[j]-r[j];

		rr[j]=r[j];

	}

	//atimes(r,rr,0);

	if (itol == 1) {

		bnrm=snrm(b,itol);

		asolve(r,z,0, ija_p, sa_p);

	}

	else if (itol == 2) {

		asolve(b,z,0, ija_p, sa_p);

		bnrm=snrm(z,itol);

		asolve(r,z,0, ija_p, sa_p);

	}

	else if (itol == 3 || itol == 4) {

		asolve(b,z,0, ija_p, sa_p);

		bnrm=snrm(z,itol);

		asolve(r,z,0, ija_p, sa_p);

		znrm=snrm(z,itol);

	} else nrerror("illegal itol in linbcg");

	cout << fixed << setprecision(6);

	while (iter < itmax) {

		++iter;

		asolve(rr,zz,1, ija_p, sa_p);

		for (bknum=0.0,j=0;j<n;j++) bknum += z[j]*rr[j];

		if (iter == 1) {

			for (j=0;j<n;j++) {

				p[j]=z[j];

				pp[j]=zz[j];

			}

		} else {

			bk=bknum/bkden;

			for (j=0;j<n;j++) {

				p[j]=bk*p[j]+z[j];

				pp[j]=bk*pp[j]+zz[j];

			}

		}

		bkden=bknum;

		atimes(p,z,0, ija_p, sa_p);

		for (akden=0.0,j=0;j<n;j++) akden += z[j]*pp[j];

		ak=bknum/akden;

		atimes(pp,zz,1, ija_p, sa_p);

		for (j=0;j<n;j++) {

			x[j] += ak*p[j];

			r[j] -= ak*z[j];

			rr[j] -= ak*zz[j];

		}

		asolve(r,z,0, ija_p, sa_p);

		if (itol == 1)

			err=snrm(r,itol)/bnrm;

		else if (itol == 2)

			err=snrm(z,itol)/bnrm;

		else if (itol == 3 || itol == 4) {

			zm1nrm=znrm;

			znrm=snrm(z,itol);

			if (fabs(zm1nrm-znrm) > EPS*znrm) {

				dxnrm=fabs(ak)*snrm(p,itol);

				err=znrm/fabs(zm1nrm-znrm)*dxnrm;

			} else {

				err=znrm/bnrm;

				continue;

			}

			xnrm=snrm(x,itol);

			if (err <= 0.5*xnrm) err /= xnrm;

			else {

				err=znrm/bnrm;

				continue;

			}

		}

//		cout << "iter=" << setw(4) << iter+1 << setw(12) << err << endl;

		if (err <= tol) break;

	}

}



DP NR::snrm(Vec_I_DP &sx, const int itol)

{

	int i,isamax;

	DP ans;



	int n=sx.size();

	if (itol <= 3) {

		ans = 0.0;

		for (i=0;i<n;i++) ans += sx[i]*sx[i];

		return sqrt(ans);

	} else {

		isamax=0;

		for (i=0;i<n;i++) {

			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;

		}

		return fabs(sx[isamax]);

	}

}



//extern Vec_INT *ija_p;

//extern Vec_DP *sa_p;



void NR::atimes(Vec_I_DP &x, Vec_O_DP &r, const int itrnsp, Vec_INT* ija_p, Vec_DP* sa_p)

{

	if (itrnsp) sprstx(*sa_p,*ija_p,x,r);

	else sprsax(*sa_p,*ija_p,x,r);

}





//extern Vec_INT *ija_p;

//extern Vec_DP *sa_p;



void NR::asolve(Vec_I_DP &b, Vec_O_DP &x, const int itrnsp, Vec_INT* ija_p, Vec_DP* sa_p)

{

	int i;



	int n=b.size();

	for(i=0;i<n;i++) x[i]=((*sa_p)[i] != 0.0 ? b[i]/(*sa_p)[i] : b[i]);

}





void NR::sprstx(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b)

{

	int i,j,k;



	int n=x.size();

	if (ija[0] != (n+1))

		nrerror("mismatched vector and matrix in sprstx");

	for (i=0;i<n;i++) b[i]=sa[i]*x[i];

	for (i=0;i<n;i++) {

		for (k=ija[i];k<ija[i+1];k++) {

			j=ija[k];

			b[j] += sa[k]*x[i];

		}

	}

}

void NR::sprsax(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b)
{
 int i, k;
 int n = x.size();
 if(ija[0] != n+1)
  nrerror("sprsax: mismatched vector and matrix");
 for (i = 0; i<n; i++) {
  b[i] = sa[i]*x[i];
  for(k=ija[i]; k<ija[i+1]; k++) {
     b[i] += sa[k]*x[ija[k]];
     }
 }

}


