version 13

drop _all
mata: mata clear

/* Load all functions needed to generate data and estimate it. The functions will be eventually compiled */

cd "C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox"

/* I commented tis file because its contents are already in the library */
run C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\stata_wrappers.ado

run generate_CNOP.ado

// use mata only to generate dataset
mata

rseed(42+15)

n = 200; k = 5; ncat = 5; infcat = 3;
beta 	= (2 \ 3 \ 0 \ 0 \ 0.0);	 	a = (-0.4 \ 1.2);
gammap 	= (0 \ 1 \ 2 \ 0 \ 0.1);	
gamman 	= (0 \ 0 \ 3 \ 5 \ 0.2);	
mup = (-3 \ 2) 
mun = (-2 \ 1) 
genCNOP(x=., y=., q=., n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun)
tmp = st_addvar("double", ("x":+strofreal(1..k), "y"))
st_addobs(n)
st_view(data = ., ., .)
data[,] = (x,y)

end


/*
// now try to estimate

cnop y x1 x2 x3 x4 x5, infcat(3)

//cnop y x1 x2, infcat(3) zp(x2 x3 x5) zn(x3 x4 x5)

//cnop y x1 x2, infcat(3) zp(x2 x3 x5) zn(x3 x4 x5) correlated

// estimate ME at means (for the latest estimate model)
cnopmargins

// estimate ME, when some coordinates are given by user (the rest is at means)
cnopmargins, at(x1 = 0.5, x3 = -0.1)

// will create variables p_1, ..., p_5 with predicted probabilities for each observation
predict p
*/

miop y x1 x2 x3 x4 x5, infcat(3) correlated
/**/

