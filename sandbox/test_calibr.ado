version 12

mata:
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       MANUAL INPUT      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

/*                           generate artificial covariates                                 */


// number of observations
n = 2000

// inflated category
infcat = 0

x1 = rnormal(n,1,0,1):+2	// mimics gdp growth rate
x2 = rnormal(n,1,0,1)		// mimics change in inflation 
x3 = runiform(n,1)
x3 = (x3:>0.7) - (x3:<0.3)	// mimics CB policy bias
x4 = rnormal(n,1,0,1)

/*                   enter all variables WITHOUT a constant                    */
    /* no overlap       
xv		= (1,0,0,0);
znv		= (0,1,0,0);
z0v		= (0,0,1,0);
zpv		= (0,0,0,1);

beta	= (.6, .4, .8, .3);
mu      = (0.7, 1.7);

betan   = (.2, .8, .3, .7);
beta0   = (.3, .1, .7, .5);
betap   = (.4, .3, .9, .2);
mun     = (-1.02, -0.25);
mu0     = (-0.91, 0.91);
mup     = (0.20, 0.82);		*/
    // partial overlap  
xv		= (1,1,0,0);
znv		= (1,0,1,0);
z0v		= (0,1,1,0);
zpv		= (0,0,1,1);

beta	= (.6, .4, .8, .3);
mu      = (0.65, 1.75);

betan   = (.2, .8, .3, .7);
beta0   = (.3, .1, .7, .5);
betap   = (.4, .3, .9, .2);
mun     = (-0.55, 0.10);
mu0     = (-0.95, 0.95);
mup     = (0.25, 1.00);	
    /* complete overlap 
xv		= (1,1,1,1);
znv		= (1,1,1,1);
z0v		= (1,1,1,1);
zpv		= (1,1,1,1);

beta	= (.6, .4, .8, .3);
mu      = (0.6, 1.8);

betan   = (.2, .8, .3, .7);
beta0   = (.3, .1, .7, .5);
betap   = (.4, .3, .9, .2);
mun     = (-1.35, -0.5);
mu0     = (-0.35, 1.55);
mup     = (1.72, 2.45);		
*/
// building model

beta 	= beta :* xv			// select coefficients in the model
betan 	= betan :* znv
beta0 	= beta0 :* z0v
betap 	= betap :* zpv

x 	= (x1, x2, x3, x4); 
/* allvars in GAUSS file consists of pointers; in contrast, x is a real matrix */

// check for multicollinearity among covariates
corrx	= correlation(x);
vnames	= strofreal(range(1, cols(x),1) * J(1,cols(x),1));	//column of strings 1, 2, ...
vnames	= "x" :+ vnames :+ ", x" :+ vnames';				// matrix
corrvec = abs(vech(corrx-I(cols(x)))):>0.1

if(sum(corrvec)>0){
	display ("CORRELATION IS greater than 0.1 BETWEEN FOLLOWING VARIABLES:");
	display(select(vech(vnames), corrvec));
}


//

eps  	= rnormal(n,1,0,1);
epsp 	= rnormal(n,1,0,1);
epsn 	= rnormal(n,1,0,1);
eps0 	= rnormal(n,1,0,1);

rs 		= x * beta' + eps;					// first layer of y
r 		= J(n, 3, 0);

r[.,1]  = rs :< mu[1];
r[.,3]  = rs :>= mu[2];
r[.,2]  = 1 :- (r[.,1] + r[.,3]);

rvec	= rowsum(r*(-1,0,1)');				// indicates -1, 0 or 1.

ysn     = x * betan' + epsn;				// second layer of y
ysp     = x * betap' + epsp;
ys0     = x * beta0' + eps0;

yn                	= J(n, cols(mun)+1,0);
yn[., 1]           	= ysn :< mun[1];
yn[., cols(mun)+1]	= ysn :>= mun[cols(mun)];

if (cols(mun) > 1){
    for (i=2; i<=cols(mun); i++){
        yn[.,i]   = (ysn :>= mun[i-1]) :* (ysn :< mun[i]);
	}
}


yp                	= J(n, cols(mup)+1,0);
yp[.,1]           	= ysp :< mup[1];
yp[., cols(mup)+1]	= ysp :>= mup[cols(mup)];

if (cols(mup) > 1){
    for (i=2; i<=cols(mup); i++){
        yp[.,i]   = (ysp :>= mup[i-1]) :* (ysp :< mup[i]);
	}
}


y0                = J(n, 3, 0);
y0[.,1]           = ys0 :< mu0[1];
y0[., 3] 		  = ys0 :>= mu0[2];
y0[., 2] 		  = 1 :- (y0[.,1] + y0[.,3]);

yn      = yn*range(-cols(mun),0,1);
yp      = yp*range(0,cols(mup),1);
y0      = y0*range(-1,1,1);

yr      = (yn,y0,yp);
y       = rowsum(r :* yr);	// observed y 

ycat	= uniqrows(y);		// column -2 -1 0 1 2 
infcat	= select(range(1,rows(ycat),1), ycat :==infcat);	// implementation of indexOf

capJ	= rows(ycat);
capJp   = capJ - (infcat - 1);
capJn   = infcat;

yo      = y;            	// to keep original values of y 

// recoding y to positive values
for (i=1; i<=capJ;i++){
    for(i2=1; i2 <= rows(y); i2++){
        if(yo[i2] == ycat[i]){
            y[i2] = i; 
        }
    }
}


yocat   = uniqrows(yo);
st_numscalar("infcat", yocat[infcat]);				// send scalar to stata

// here data must have been saved on disk, but I instead load it to the Stata dataset

end							// now switch from mata to stata!

getmata y=yo x1=x1 x2=x2 x3=x3 x4=x4 r=rvec, replace force	// copy all variables to Stata interface
disp "	New variables y x1 x2 x3 x4 r have been created"
disp "         Sample descriptive statistics"
sum y x1 x2 x3 x4
disp "Frequency distribution of discrete categories of y"
tab y
disp "Inflated category of y is "	infcat
disp "                 Decomposition of zeros"
tab r if y==0
disp "                 Decomposition of -1"
tab r if y==-1
disp "                 Decomposition of 1"
tab r if y==1
