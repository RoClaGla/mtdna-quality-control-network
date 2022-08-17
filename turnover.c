#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define RND drand48()
#define PI 3.1415927

int MAXN = 10000;
int MAXM = 100;

// structure to store simulations stats
typedef struct{
	double mwc, mwn, mmc, mmn, mh, mu, md, mpnet, mpcyt, mmprop;
	double vwc, vwn, vmc, vmn, vh, vu, vd, vpnet, vpcyt, vmprop;
} SumStats;

// structure to store stats from distinct simulations
typedef struct{
	double wc, mc, wn, mn;
	double u, d, het, pnet, pcyt, mprop;
} Stats;

// GSL routine for gaussian random number with sd sigma
double gsl_ran_gaussian(const double sigma)
{
	double x, y, r2;

	do
	{
		/* choose x,y in uniform square (-1,-1) to (+1,+1) */

		x = -1 + 2 * RND;
		y = -1 + 2 * RND;

		/* see if it is in the unit circle */
		r2 = x * x + y * y;
	} while (r2 > 1.0 || r2 == 0);

	/* Box-Muller transform */
	return sigma * y * sqrt(-2.0 * log(r2) / r2);
}

// simple function returning max(a,b)
double mymax(double a, double b)
{
  if(a > b) return a;
  else return b;
}

// simple function returning min(a,b)
double mymin(double a, double b)
{
  if(a < b) return a;
  else return b;
}

// this set of functions is taken from https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/ and computes whether two line segments intersect
// this one computes whether q lies on segment pr
int onSegment(double px, double py, double qx, double qy, double rx, double ry)
{
  if(qx <= mymax(px, rx) && qx >= mymin(px, rx) && qy <= mymax(py, ry) && qy >= mymin(py, ry)) return 1;
  else return 0;
}

// this one computes the orientation of a set of three points (or whether they are colinear)
int orientation(double px, double py, double qx, double qy, double rx, double ry)
{
  double val = (qy-py)*(rx-qx) - (qx-px)*(ry-qy);
  if(val == 0) return 0;
  return (val > 0 ? 1 : 2);
}

// this one returns 1/0 intersection given two line segments described by start and end points 
int doIntersect(double p1x, double p1y, double q1x, double q1y, double p2x, double p2y, double q2x, double q2y)
{
  int o1 = orientation(p1x, p1y, q1x, q1y, p2x, p2y);
  int o2 = orientation(p1x, p1y, q1x, q1y, q2x, q2y);
  int o3 = orientation(p2x, p2y, q2x, q2y, p1x, p1y);
  int o4 = orientation(p2x, p2y, q2x, q2y, q1x, q1y);
  
  // General case
  if (o1 != o2 && o3 != o4)
    return 1;
  
  // Special Cases
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1
  if (o1 == 0 && onSegment(p1x, p1y, p2x, p2y, q1x, p1y)) return 1;
  
  // p1, q1 and q2 are colinear and q2 lies on segment p1q1
  if (o2 == 0 && onSegment(p1x, p1y, q2x, q2y, q1x, q1y)) return 1;
  
  // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  if (o3 == 0 && onSegment(p2x, p2y, p1x, p1y, q2x, q2y)) return 1;
  
  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  if (o4 == 0 && onSegment(p2x, p2y, q1x, q1y, q2x, q2y)) return 1;
  
  return 0;
}

int Output(double *xs, double *ys, double *xe, double *ye, double *mx, double *my, int *mt, int n, int nsegs, int nseed, double mass, double p, double q, double halo, double rho){
  int i;
  FILE *fp;
  char str[200];
  
  sprintf(str,"network-%i-%i-%.f-%.2f-%.2f-%.2f-%.2f.csv",n,nseed,mass,p,q,halo,rho);
  fp = fopen(str,"w");
  fprintf(fp,"xs,ys,xe,ye\n");
  for(i=0;i<nsegs;i++)
    fprintf(fp,"%f,%f,%f,%f\n",xs[i],ys[i],xe[i],ye[i]);
  fclose(fp);
  
  sprintf(str,"mtdna-%i-%i-%.f-%.2f-%.2f-%.2f-%.2f.csv",n,nseed,mass,p,q,halo,rho);
  fp = fopen(str,"w");
  fprintf(fp,"x,y,type\n");
  for(i=0;i<n;i++)
    fprintf(fp,"%f,%f,%i\n",mx[i],my[i],mt[i]);
  fclose(fp);
  return(0);
}

int BuildNetwork(double *xs,double *ys,double *xe,double *ye, double mass, int nseed, double seglength, double branchprob, double theta, int *nsegs, double *actual_mass){
  double newx, newy, angle, phaseshift;
  int i,j,k,nactive, doesintersect, outofbounds;
  double total, thism;
  int *active;
  double *thetas;
  
  active = (int *)malloc(sizeof(int)*MAXN);
  thetas = (double *)malloc(sizeof(double)*MAXN);
  
  j = 0;
  nactive = 0; // if starting off, or all branches become inactivated, we uniformly reseed the perimeter
  for(;total<mass;){
    // make seeds regularly spaced along the perimeter of the unit disc, starting at a random angle
    if(nactive == 0){
      phaseshift = 2*PI*RND;
      for(i=0;i<nseed;i++){
        xs[j] = xe[j] = cos(2*PI/nseed*i+phaseshift);
        ys[j] = ye[j] = sin(2*PI/nseed*i+phaseshift);
        thetas[j] = PI + 2*PI/nseed*i+phaseshift; // angles of growth are rotated by PI, pointing inwards
        active[j] = 1; nactive++; j++;
      }
    }
    // for each branch
    for(i=0;i<j;i++){
      // if it is active
      if(active[i]==1){
        //decide if the branch should split
        if(RND>branchprob){
          // no split, so grow current branch by seglength
          newx = xe[i] + seglength*cos(thetas[i]);
          newy = ye[i] + seglength*sin(thetas[i]);
          // if out of bounds, or intersects another segment, terminate the branch:
          doesintersect = 0;
          for(k=0;k<j;k++){
						//printf("Test for k - %i\n",doIntersect(xs[i],ys[i],newx,newy,xs[k],ys[k],xe[k],ye[k]));
            if(doIntersect(xs[i],ys[i],newx,newy,xs[k],ys[k],xe[k],ye[k])){
              doesintersect = 1;
            } 
          }
					outofbounds = 0;
          if(newx*newx+newy*newy>1){
						outofbounds = 1;
					}
					if(outofbounds == 1 || doesintersect == 1){
						nactive--;
            active[i] = 0;
          }else{
            // change the endpoint of the ith branch from its current one to this one
            xe[i] = newx;
            ye[i] = newy;
						
            total += seglength;
            if(total>mass){
              *nsegs = j;
              *actual_mass = total;
              free(thetas);
              free(active);
              return(0);
            }
          }
        }else{
          // kill the current branch, whose angle is thetas[i], and spawn two branches with new direction
          active[i] = 0;
          xs[j] = xe[j] = xe[i];
          ys[j] = ye[j] = ye[i];
          thetas[j] = theta + gsl_ran_gaussian(0.1);
          active[j] = 1;
          j++;

          xs[j] = xe[j] = xe[i];
          ys[j] = ye[j] = ye[i];
          thetas[j] = theta + gsl_ran_gaussian(0.1);
          active[j] = 1;
          
          // increment nactive by 1 (2- 1 net branches are spawned)
          nactive++;
          j++;
          // avoid processing new segments this time step
          i+=2;;
        }
      }
    }
  }
  
  *nsegs = j;
  *actual_mass = total;
  free(thetas);
  free(active);
  return(0);
}

int PlaceDNA(double *xs, double *ys, double *xe, double *ye, double *mx, double *my, int *mt, int *mnetworked, int n, double h, double p, double q, int nsegs, double halo){
  int i,j,k, placed, failed, counter, notdoneyet;
  double *cumsum;
  double mtx, mty, r, ball;
  
  cumsum = (double *)malloc(sizeof(double)*nsegs);
  
  cumsum[0] = sqrt((xs[0]-xe[0])*(xs[0]-xe[0])+(ys[0]-ye[0])*(ys[0]-ye[0]));
  for(i=1;i<nsegs;i++){
    cumsum[i] = cumsum[i-1] + sqrt((xs[i]-xe[i])*(xs[i]-xe[i])+(ys[i]-ye[i])*(ys[i]-ye[i]));
  }  
  
  // allocate mtDNA types:
  for(i=0;i<n;i++)mt[i] = 0; // WT
  for(i=0;i<h*n;i++){
    do{k = RND*n;}while(mt[k] == 1);
    mt[k] = 1;// MU
  }
  
  // allocate mtDNA molecules
  notdoneyet = 0;
  for(i=0;i<n;i++){
    placed = 0;
    if((RND<p && mt[i] == 0)||(RND<q && mt[i] == 1)){
      // place in the network
      mnetworked[i] = 1;
      counter = 0;
      while(placed == 0 && notdoneyet == 0){
        ball = RND*cumsum[nsegs-1]; // roulette wheel with cumulative masses to choose a strand;
        for(j=0;cumsum[j]<ball;j++);
        
        r = RND; // random distance along the strand
        
        // store the proposed position
        mtx = xs[j] + (xe[j]-xs[j])*r;
        mty = ys[j] + (ye[j]-ys[j])*r;

        // and check if is consistent with our rules
        failed = 0;
        // for each (networked) mtDNA molecules already placed
        for(k=0;k<i;k++){
          if(mnetworked[k]==1){
            if((mx[k]-mtx)*(mx[k]-mtx)+(my[k]-mty)*(my[k]-mty)<halo*halo){
              failed = 1;
              break;
            }
          }
        }
        // if not failed, place it here; otherwise try again
        if(failed == 0){
          mx[i] = mtx;
          my[i] = mty;
          placed = 1;
        }
        counter++;
        // if we have tried too many times with this mtDNA
        if(counter>10000){
          printf("Stuggling to place mito %i!\n",i);
          notdoneyet = 1;
          break;
        }
      }
    }else{
      mnetworked[i] = 0;
      // place randomly in the cytoplasm
      while(placed == 0){
        mtx = -1+2*RND;
        mty = -1+2*RND;
        if(mtx*mtx+mty*mty<1){
          mx[i] = mtx;
          my[i] = mty;
          placed = 1;
        }
      }
    }
    //printf("i,mt[i],mnetworked[i]=%i,%i,%i\n",i,mt[i],mnetworked[i]);
  }
  free(cumsum);
  if(notdoneyet == 1) return(1);
  return(0);
}

// function spatially correlates K DNA mutants with a random point P on the perimeter
void correlateDNA(double *mx, double *my, int *mt, int n, int K){
	double x,y,r;
	double thisdist,mindist,ball;
	int i,j,k,l, mini, iinset, *mmi;
	
	// keep track of the recorded indices of the closest ones
	mmi = (int *)malloc(sizeof(int)*K);
	
	// Pick a random point P on the perimeter
	r = 2*PI*RND;
	x = cos(r); y = sin(r);
	
	// make sure the K mitos closest to P are mutants
	l = 0; // keep track of the number that "mutate"
	while(k<K){
		// loop through the population
		mindist = 10;
		for(i=0;i<n;i++){
			// check whether this i is already recorded, in which case, do nothing
			iinset = 0;
			for(l=0;l<k;l++){
				if(mmi[l] == i)iinset = 1;
			}
			if(iinset==0){
				// if not in set, record its distance from P
				thisdist = sqrt((x-mx[i])*(x-mx[i])+(y-my[i])*(y-my[i]));
				if(thisdist<mindist){
					// record this as minimal, with its index
					mindist = thisdist;
					mini = i;
				}
			}
		}
		// record its index
		mmi[k] = mini;
		// if not mutant, turn it mutant, and record that we did so
		if(mt[mini] == 0){
			mt[mini] = 1;
			j++;
		}
		k++;
	}
	// now turn j of the remaining mutants into wildtypes
	l = 0;
	while(l<j){
		// select a random index between 0 and n-1
		ball = RND*n;
		for(i=0;i<ball;i++);
		// check whether this i is already recorded, in which case, do nothing
		iinset = 0;
		for(l=0;l<k;l++){
			if(mmi[l] == i)iinset = 1;
		}
		if(iinset==0 && mt[i]==1){
			mt[i] = 0;
			l++;
		}
	}
	free(mmi);
}

// counts the number of mutants with a wildtype within rho
int getMutantProp(double rho, double *mx, double *my, int *mt, int n, double *propwithin){
	int i, j;
	double mus, pxDNA;
	
	pxDNA = 0;
	mus = 0;
	for(i=0;i<n;i++){
		// if mutant, loop through all the others
		if(mt[i] == 1){
			mus++;
			for(j=1;j<n;j++){
				// is this a wildtype and is the distance smaller than rho?
				if(mt[j]==0 && ((mx[i]-mx[j])*(mx[i]-mx[j])+(mx[i]-mx[j])*(mx[i]-mx[j])<rho*rho)){
					pxDNA++;
					break;
				}
			}
		}
	}
	
	*propwithin = (double) pxDNA/((double) mus);
	return(0);
}

// function calculates the mean number of DNAs within rho of each other in the top daughter
int getSeparateProximalDNA(double rho, double *mx, double *my, int *mt, int *mnetworked, int n, double *mproxnet, double *mproxcyt){
  int i,j, pxDNA, k1, k2;
	double *netDNA, *cytDNA;
	double meannet, meancyt;
	
	netDNA = (double *)malloc(sizeof(double)*n);
	cytDNA = (double *)malloc(sizeof(double)*n);
	
	// compute DNAs within a radius of influence of rho of each DNA in the network of the top daughter
	k1 = 0;
  for(i=0;i<n;i++){
		pxDNA = 0;
    for(j=0;j<n;j++){
			if(i!=j){
				if(my[i]>0&&my[j]>0){
					if(mnetworked[i] == 1 && mnetworked[j] == 1){
						if((mx[j]-mx[i])*(mx[j]-mx[i])+(my[j]-my[i])*(my[j]-my[i])<rho*rho){
							pxDNA++;
							k1++;
						}
					}
				}
			}
    }
		netDNA[i] = pxDNA;
  }
	
	// compute DNAs within a radius of influence of rho of each DNA in the cytoplasm of the top daughter
	k2 = 0;
  for(i=0;i<n;i++){
		pxDNA = 0;
    for(j=0;j<n;j++){
			if(i!=j){
				if(my[i]>0&&my[j]>0){
					if(mnetworked[i] == 0 && mnetworked[j] == 0){
						if((mx[j]-mx[i])*(mx[j]-mx[i])+(my[j]-my[i])*(my[j]-my[i])<rho*rho){
							pxDNA++;
							k2++;
						}
					}
				}
			}
    }
		cytDNA[i] = pxDNA;
  }
	
	// mean in the network
	meannet = 0;
	for(i=0;i<k1;i++){
		meannet += netDNA[i];
	}
	if(k1>0)meannet /= k1;
	else meannet = 0;
	
	// mean in the cytoplasm
	meancyt = 0;
	for(i=0;i<k2;i++){
		meancyt += cytDNA[i];
	}
	if(k2>0)meancyt/=k2;
	else meancyt =0;
	
	*mproxnet = meannet;
	*mproxcyt = meancyt;
	
	free(netDNA);
	free(cytDNA);
	return(0);
}

// function to go through k "cell cycles": "fragment" network, do turnover, "refuse" (identical network for now)
int Cycle(double *mx, double *my, int *mt, int *mnetworked, int n, int K, double rho, double mut_rate, double to_rate){
	int i,j,k;
	int turnover;
	
	for(k=0;k<K;k++){
		// for each mutant, subject to turnover if no wildtype is within rho of it
		for(i=0;i<n;i++){
			if(mt[i] == 1){
				turnover = 1;
				for(j=0;j<n;j++){
					if(i!=j){
						if((mx[i]-mx[j])*(mx[i]-mx[j])+(my[i]-my[j])*(my[i]-my[j])<rho*rho){
							if(mt[j] == 0){
								turnover = 0;
							}
						}
					}
				}
				if(RND<to_rate && turnover == 1){
						mt[i] = 0;
				}
			}
		}
		// for each wildtype, chance of becoming mutant if not in the network (supposing it needs more work when not complemented by other mtDNA)
		for(i=0;i<n;i++){
			if(mt[i] == 0 && RND<mut_rate && mnetworked[i] == 0){
				mt[i] = 1;
			}
		}
	}
	return(0);
}

// fetches stats for a particular simulation
int getStats(double *mx, double *my, int *mt, int *mnetworked, int n, int *wc, int *mc, int *wn, int *mn, double *het){
	int i,wcyt,mcyt,wnet,mnet;
	double cellhet;
	
	wcyt = mcyt = wnet = mnet = 0;
	for(i=0;i<n;i++){
		if(my[i]>0){
			if(mnetworked[i] == 0 && mt[i] == 0){
				wcyt++;
			}
			if(mnetworked[i] == 0 && mt[i] == 1){
				mcyt++;
			}
			if(mnetworked[i] == 1 && mt[i] == 0){
				wnet++;
			}
			if(mnetworked[i] == 1 && mt[i] == 1){
				mnet++;
			}
		}
	}
	
	cellhet = ((double)mcyt+(double)mnet)/((double)mcyt+(double)mnet+(double)wcyt+(double)wnet);
	*wc = wcyt; *mc = mcyt; *wn = wnet; *mn = mnet;
	*het = cellhet;
	return(0);
}

int getNetworkProp(double *xs, double *ys, double *xe, double *ye, double nsegs, double *u){
	int i;
	double propupper;
	double mass;
	
	mass = 0;
	for(i=0;i<nsegs;i++){
		if(ys[i]>0&&ye[i]>0)propupper = 1;
		else if(ys[i]>0&&ye[i]<0)propupper = ys[i]/(ys[i]-ye[i]);
		else if(ys[i]<0&&ye[i]>0)propupper = ye[i]/(ye[i]-ys[i]);
		else propupper = 0;
		
		if(propupper > 1 || propupper < 0){
			printf("wtf!");
		}
		mass += sqrt((xs[i]-xe[i])*(xs[i]-xe[i])+(ys[i]-ye[i])*(ys[i]-ye[i]))*propupper;
	}
	*u = mass;
	return(0);
}

int getMinDNASeparation(double *mx, double *my, int n, double *d){
	int i, j;
	double mindist, thisdist;
	
	mindist = 10;
	for(i=0;i<n;i++){
		for(j=0;i<n;i++){
			if(i!=j){
				thisdist = sqrt((mx[i]-mx[j])*(mx[i]-mx[j])+(my[i]-my[j])*(my[i]-my[j]));
				if(thisdist<mindist)mindist = thisdist;
			}
		}
	}
	
	*d = mindist;
	return(0);
}

int computeStats(Stats *s, SumStats *ss, int nsims){
	int i;
	double mnet,vnet,mcyt,vcyt,het,hetvar;
	double mwc,mmc,mwn,mmn;
	double vwc,vmc,vwn,vmn;
	double mmass,vmass;
	double mmd,vmd;
	double mmprop,vmprop;
	
	// calculate means and variances of the number of nearby neighbors for this DNA distribution
	mnet = mcyt = vnet = vcyt = 0;
	for(i=0;i<nsims;i++){
		mnet+=s[i].pnet;
		mcyt+=s[i].pcyt;
	}
	mnet/=nsims;
	ss->mpnet = mnet;
	mcyt/=nsims;
	ss->mpcyt = mcyt;
	
	for(i=0;i<nsims;i++){
		vnet+=(s[i].pnet-mnet)*(s[i].pnet-mnet);
		vcyt+=(s[i].pcyt-mcyt)*(s[i].pcyt-mcyt);
	}
	
	ss->vpnet = vnet/(nsims-1);
	ss->vpcyt = vcyt/(nsims-1);

	// calculate number of wts, number of mus, heteroplasmy level in the upper half
	mwc=mmc=mwn=mmn=0;
	for(i=0;i<nsims;i++){
		mwc += s[i].wc;
		mmc += s[i].mc;
		mwn += s[i].wn;
		mmn += s[i].mn;
	}
	ss->mwc = mwc/nsims;
	ss->mmc = mmc/nsims;
	ss->mwn = mwn/nsims;
	ss->mmn = mmn/nsims;
	
	vwc=vmc=vwn=vmn=0;
	for(i=0;i<nsims;i++){
		vwc += (s[i].wc-mwc)*(s[i].wc-mwc);
		vmc += (s[i].mc-mmc)*(s[i].mc-mmc);
		vwn += (s[i].wn-mwn)*(s[i].wn-mwn);
		vmn += (s[i].mn-mmn)*(s[i].mn-mmn);
	}
	
	ss->vwc = vwc/(nsims-1);
	ss->vmc = vmc/(nsims-1);
	ss->vwn = vwn/(nsims-1);
	ss->vmn = vmn/(nsims-1);
	
	het = 0;
	for(i=0;i<nsims;i++){
		het+=s[i].het;
	}
	het/=nsims;
	ss->mh = het;
	
	hetvar = 0;
	for(i=0;i<nsims;i++){
		hetvar+=(s[i].het-het)*(s[i].het-het);
	}
	ss->vh = hetvar/(nsims-1);
	
	mmass = vmass = 0;
	for(i=0;i<nsims;i++){
		mmass += s[i].u;
	}
	mmass /= nsims;
	ss->mu = mmass;
	
	for(i=0;i<nsims;i++){
		vmass += (s[i].u-mmass)*(s[i].u-mmass);
	}
	ss->vu = vmass/(nsims-1);
	
	mmd = vmd = 0;
	for(i=0;i<nsims;i++){
		mmd += s[i].d;
	}
	mmd/=nsims;
	ss->md = mmd;
	for(i=0;i<nsims;i++){
		vmd+=(s[i].d-mmd)*(s[i].d-mmd);
	}
	ss->vd = vmd/(nsims-1);
	
	mmprop = vmprop = 0;
	for(i=0;i<nsims;i++){
		mmprop+=s[i].mprop;
	}
	mmprop/=nsims;
	ss->mmprop = mmprop;
	for(i=0;i<nsims;i++){
		vmprop+=(s[i].mprop-mmprop)*(s[i].mprop-mmprop);
	}
	ss->vmprop=vmprop/(nsims-1);
	
	return(0);
}

int main(int argc, char *argv[]){
  double *xs,*xe,*ys,*ye, *mx, *my;
  int *mt, *mnetworked;
  double rho, h, p, q, seglength, branchprob, halo, theta;
  double target_mass, actual_mass;
	double mproxnet, mproxcyt;
	int wc,mc,wn,mn;
	double het, mhet, vhet, proxnet,proxcyt,vproxnet,vproxcyt;
	double mwc, mmc, mwn, mmn,vwc,vmc,vwn,vmn,u;
	double mprop, vprop, mmind, vmind;
	//double mut_rate, to_rate; // mutation and turnover rates; 
	//int t,tmax; // maximal number of "cycles" of turnover/mtDNA mutation before cell division
	Stats *S;
	SumStats Ss;
  int error, nsims, nsim, n, nseed, nsegs, output, notdoneyet, K;
  char str[200];
  FILE *fp;
  
  error = 1;
  printf("argc = %i\n",argc);
  if(argc == 15 || argc == 8){
    // computational parameters and parameters in common
    nsims = atoi(argv[2]);
    // network parameters
    target_mass = atof(argv[3]);
    seglength = atof(argv[4]);
    branchprob = atof(argv[5]);
		theta = atof(argv[6]);
    // genetic parameters
    n = atoi(argv[7]);
    if(argc == 15 && strcmp(argv[1],"--snapshots\0")==0){
      // generate snapshot
      error = 0;
      output = 1;
			h = atof(argv[8]);
      nseed = atoi(argv[9]);
      halo = atof(argv[10]);
      p = atof(argv[11]);
      q = atof(argv[12]);
      rho = atof(argv[13]);
			K = atoi(argv[14]);
      printf("Parameters nsims, nseed, target_mass, seglength, branchprob, h, n, halo, p, q: %i, %i, %.f, %.2f, %.2f, %.2f, %i, %.2f, %.2f, %.2f\n", nsims,nseed,target_mass,seglength,branchprob,h,n,halo,p,q);
    }
    if(argc == 7 && strcmp(argv[1],"--simulate\0")==0){
      // no visual output
      error = 0;
      output = 0;
      printf("Parameters nsims, target_mass, seglength, branchprob, h, n: %i, %.f, %.2f, %.2f, %.2f, %i\n", nsims, target_mass,seglength,branchprob,h,n);
    }
  }
  //printf("error = %i\n",error);
  if(error == 1){
    printf("Error. Usage ./network.ce --snapshots [nsims] [nseed] [seglength] [branchprob] [p] [q]\n");
    printf("./network.ce --simulate [nsims] [nseed] [seglength] [branchprob]\n");
    exit(0);
  }
	
	//tmax = 4;
	theta = PI/3;
  
  xs = (double *)malloc(sizeof(double)*MAXN);
  ys = (double *)malloc(sizeof(double)*MAXN);
  xe = (double *)malloc(sizeof(double)*MAXN);
  ye = (double *)malloc(sizeof(double)*MAXN);

  mx = (double *)malloc(sizeof(double)*MAXM);
  my = (double *)malloc(sizeof(double)*MAXM);
  mt = (int *)malloc(sizeof(int)*MAXM);
  mnetworked = (int *)malloc(sizeof(int)*MAXM);
	
	S = (Stats *)malloc(sizeof(Stats)*nsims);
	
  // decide whether these are cell snapshots or simulations (below)
  if(output == 1){
    notdoneyet = 1;
    while(notdoneyet == 1){
      BuildNetwork(xs,ys,xe,ye,target_mass,nseed,seglength,branchprob,theta,&nsegs,&actual_mass);
      notdoneyet = PlaceDNA(xs,ys,xe,ye,mx,my,mt,mnetworked,n,h,p,q,nsegs,halo);
    }
		correlateDNA(mx,my,mt,n,K);
    Output(xs,ys,xe,ye,mx,my,mt,n,nsegs,nseed,target_mass,p,q,halo,rho);
    return(0);
  }else{		
		sprintf(str,"output2-%i.csv",n);
    fp = fopen(str,"w");
    fprintf(fp,"h,n,nseed,p,q,halo,rho,mpnet,mpcyt,vpnet,vpcyt,mmprop,vmprop,mwc,vwc,mmc,vmc,mwn,vwn,mmn,vmn,mh,vh,mu,vu,md,vd\n");
    for(h=0.1;h<=0.5;h+=0.4){
			for(nseed=4;nseed<=64;nseed*=4){
				for(p=0.0;p<=1.0;p+=0.1){
					for(q=0.0;q<=1.0;q+=0.1){
						for(halo=0;halo<=0.1;halo+=0.1){
							for(rho=0.0;rho<=.3;rho+=0.01){
								//for(K=0;K<15;K+=5){
								//for(mut_rate=0.0;mut_rate<=0.05;mut_rate+=0.025){
								//for(to_rate=0.0;to_rate<=0.05;to_rate+=0.025){
								//for(t=0;t<tmax;t++){
								nsim = 0;
								while(nsim<nsims){
									//printf("nsim,nseed,p,q,halo,rho = %i,%i,%.2f,%.2f,%.2f,%.2f\n",nsim,nseed,p,q,halo,rho);
									notdoneyet = 1;
									while(notdoneyet == 1){
										//printf("New attempt\n");
										BuildNetwork(xs,ys,xe,ye,target_mass,nseed,seglength,branchprob,theta,&nsegs,&actual_mass);
										notdoneyet = PlaceDNA(xs,ys,xe,ye,mx,my,mt,mnetworked,n,h,p,q,nsegs,halo);
									}

									// turnover according to parameterisation and number of turnover occasions:
									//Cycle(mx,my,mt,mnetworked,n,t,rho,mut_rate,to_rate);
									// correlate DNA according to cluster size K
									//correlateDNA(mx,my,mt,n,K);
									// get DNA stats
									getStats(mx,my,mt,mnetworked,n,&wc,&mc,&wn,&mn,&het);
									getMutantProp(rho,mx,my,mt,n,&mprop);
									// fix so that all functions below are called, pass Stats directly to getStats
									getSeparateProximalDNA(rho, mx, my, mt, mnetworked, n, &mproxnet, &mproxcyt);								
									getNetworkProp(xs,ys,xe,ye,nsegs,&u);
									getMinDNASeparation(mx,my,n,&mmind);
									S[nsim].wc = wc;
									S[nsim].mc = mc;
									S[nsim].wn = wn;
									S[nsim].mn = mn;
									S[nsim].het = het;
									S[nsim].u = u;
									S[nsim].d = mmind;
									S[nsim].pnet = mproxnet;
									S[nsim].pcyt = mproxcyt;
									S[nsim].mprop = mprop;
									nsim++;
								}
								// compute stats for the given parameterisation:
								computeStats(S,&Ss,nsims);
								// for later: chose if the network remains equally heterogeneous throughout, or if we randomly draw heterogeneity of network
								// bump to output file
								fprintf(fp,"%.2f,%i,%i,%.2f,%.2f,%.2f,%.3f,%f,%f,%f,%f,%f,%f,%f,%.2e,%f,%.2e,%f,%.2e,%f,%.2e,%f,%f,%f,%f,%f,%f\n",h,n,nseed,p,q,halo,rho,Ss.mpnet,Ss.mpcyt,Ss.vpnet,Ss.vpcyt,Ss.mmprop,Ss.vmprop,Ss.mwc,Ss.vwc,Ss.mmc,Ss.vmc,Ss.mwn,Ss.vwn,Ss.mmn,Ss.vmn,Ss.mh,Ss.vh,Ss.mu,Ss.vu,Ss.md,Ss.vd);
								fflush(fp);
								printf("Should print!\n");
								//}
								//}
								//}
								//}
							}
						}
					}
				}
			}
		}
    fclose(fp);
    return(0);
  }
}
