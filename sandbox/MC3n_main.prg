/* сей код может генерировать данные одной из трёх моделей (OP, NOP or CNOP),
затем оценивает все три модели, выдаёт все результаты на экран и в
файл MC3***.out и сохраняет dataset in the last MC iteration в файле MC3_data */

new; 
closeall; cls;

library cml, pgraph;
#include cml.ext;
#include C:\Users\David\YandexDisk\hsework\gauss-mata\mc\MC3n_aux.prc;

cmlset;
workdir     = cdir(0);
load path 	= ^workdir;
save path 	= ^workdir;

t0     = timeutc;
cv     = cdfni(0.975);

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       MANUAL INPUT      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

No          = 1000;    /* choose sample size */
rep         = 1;   /* number of generated samples */
meNo        = 1;    /* enter No of observations to compute marginal effects for their values of covariates or 0 if you don't want to compute ME  */

/*                                  choose the dgp                                          */

let model	= nop;     /* 2: Nested Ordered Probit */
let model	= op;      /* 1: Ordinary Ordered Probit */
let model	= tiop;    /* 4: Triple-Zero-Inflated Ordered Probit (CNOP) */

let case    = no;
let case    = complete;
let case    = partial;

let overlap = no partial complete;
overlap     = indcv(case,overlap); 
let dgp     = op nop tiop; 
dgp         = indcv(model,dgp);

/*screen off;*/
screen on;

cJ     = 5;      /* enter the number of outcome categories */
inf    = 3;      /* indicate which category of y (1,2,3,...) is to be inflated */
xbar   = 1;      /* enter values of covariates (observations in rows) to compute MEs; set to 0 for zeros, 1 for population means  */
vsop   = 1;      /* enter 0 if you don't request searching for higher likelihood for larger nested models than for smaller ones, and 1 otherwise */
dummy  = 3;      /* enter 0 if you don't have the dummy (discrete x), or n if Xn is a dummy variable */
toler  = 1e-5;   /* set tolerance criteria for cml procedure */
bound  = 1000;   /* set bounds for parameters in cml procedure */

let sv_op   = 0; /* enter starting values, or 0 for default ones */
let sv_nop  = 0;
let sv_tiop = 0; 

/*                       load or generate artificial covariates                     */

/*if  No le 2000;
    load x_2000;
    x1 = x_2000[1:No,1]; x2 = x_2000[1:No,2]; x3 = x_2000[1:No,3]; x4 = x_2000[1:No,4];
    else; */
    x1   = rndn(No,1)+2;
    x2   = rndn(No,1);
    x3   = rndu(No,1);
    x3   = (x3 .gt 0.70) - (x3 .lt 0.30);
    x4   = rndn(No,1);
/*endif;*/

/* parameters are calibrated to produce in case of 5 outcome categories approximately:
y: -2: 7%, -1: 14%, 0: 58%, 1: 14%, 2: 7% observations
Decomposition of zeros:   r=-1: 1/3, r=0: 1/3, r=1: 1/3 
*/

if dgp eq 3; /*  TIOP */

    output file = C:\Users\David\YandexDisk\hsework\gauss-mata\mc\MC3_TIOP.out reset;
    lpwidth 256; outwidth 256; output on;

    if overlap eq 1;  @ no overlap @

    let xv   = 1 0 0 ;   
    let znv  = 0 1 0 ;
    let zpv  = 0 0 1 ;

    let beta    = .6 .4 .8 ;
    let mu      = 0.95 1.45;
    let betan   = .2 .8 .3 ;
    let betap   = .4 .3 .9 ;
    let mun     = -1.22 0.03;
    let mup     = -.03 1.18;

    elseif overlap eq 2;  @ partial overlap @

    let xv   = 1 1 0  ;   
    let znv  = 1 0 1  ;
    let zpv  = 0 1 1  ;

    let beta    = .6 .4 .8;
    let mu      = .9 1.5;
    let betan   = .2 .8 .3;
    let betap   = .4 .3 .9;
    let mun     = -.67 .36;
    let mup     = .02 1.28;

    elseif overlap eq 3;  @ complete overlap @

    let xv   = 1 1 1 ;   
    let znv  = 1 1 1 ;
    let zpv  = 1 1 1 ;

    let beta    = .6 .4 .8 ;
    let mu      = 0.85 1.55;
    let betan   = .2 .8 .3 ;
    let betap   = .4 .3 .9 ;
    let mun     = -1.2 0.07;
    let mup     = 1.28 2.5;

    endif;

    beta        = selif(beta,  xv);
    betan       = selif(betan, znv);
    betap       = selif(betap, zpv);

    par_true    = beta|mu|betan|mun|betap|mup;
    save xv, znv, zpv, beta, mu, betan, betap, mun, mup;


elseif dgp eq 2; /* NOP */

    output file = C:\Users\David\YandexDisk\hsework\gauss-mata\mc\MC3_NOP.out reset;
    lpwidth 256; outwidth 256; output on;

    if overlap eq 1;  @ no overlap @

    let beta    = .6 .4 .8 ;
    let mu      = 0.26 2.14;
    let betan   = .2 .8 .3 ;
    let betap   = .4 .3 .9 ;
    let mun     = -0.54;
    let mup     =  0.54;

    let xv   = 1 0 0 ;   
    let znv  = 0 1 0 ;
    let zpv  = 0 0 1 ;

    elseif overlap eq 2;  @ partial overlap @

    let beta    = .6 .4 .8 ;
    let mu      = 0.21 2.19;
    let betan   = .2 .8 .3 ;
    let betap   = .4 .3 .9 ;
    let mun     = -0.17 ;
    let mup     = 0.68;

    let xv      = 1 1 0 ;   
    let znv     = 1 0 1 ;
    let zpv     = 0 1 1 ;

    elseif overlap eq 3;  @ complete overlap @

    let xv   = 1 1 1 ;   
    let znv  = 1 1 1 ;
    let zpv  = 1 1 1 ;

    let beta    = .6 .4 .8 ;
    let mu      = 0.09 2.32;
    let betan   = .2 .8 .3 ;
    let betap   = .4 .3 .9 ;
    let mun     = -.72;
    let mup     =  2.12;

    endif;

    beta        = selif(beta,  xv);
    betan       = selif(betan, znv);
    betap       = selif(betap, zpv);

    par_true    = beta|mu|betan|mun|betap|mup;
    save xv, znv, zpv, beta, mu, betan, betap, mun, mup;

elseif dgp eq 1; /* OP */

    output file = C:\Users\David\YandexDisk\hsework\gauss-mata\mc\MC3_OP.out reset;
    lpwidth 256; outwidth 256; output on;

    let beta    = .6 .4 .8;
    let mu      = -1.83 -1.01  1.01  1.83;

    let xv     = 0 1 1  ; 
    let znv    = 0 0 0  ;
    let zpv    = 0 0 0  ;
  
    beta       = selif(beta, xv);
    par_true   = beta|mu;
    save xv, znv, zpv, beta, mu;

endif;

/*******************************   END of MANUAL INPUT   ********************************/

let allvars  = x1 x2 x3;
index        = (xv + zpv + znv) .ge 1;
dummy        = sumc(index[1:dummy]);
allvars      = selif(allvars,index);
xv           = selif(xv,index);
znv          = selif(znv,index);
zpv          = selif(zpv,index);

xvars        = selif(allvars,xv); 
zpvars       = selif(allvars,zpv); 
znvars       = selif(allvars,znv);
xzvars       = selif(allvars,(xv + zpv + znv ) .ge 1);

alldata      = x1~x2~x3;
seldata      = selif(alldata',index)';
xzdata       = selif(seldata', (xv + zpv + znv ) .ge 1)';
xdata        = selif(seldata',xv)';
zpdata       = selif(seldata',zpv)';
zndata       = selif(seldata',znv)';

kxz     = sumc((xv + zpv + znv ) .ge 1);
kx      = sumc(xv);
kzp     = sumc(zpv);
kzn     = sumc(znv);

/* values of covariates to compute marginal effects */

if (cols(xbar) eq 1) and (xbar eq 1); 
        let mebar = 2 0 0 ; mebar = selif(mebar,index)'; /*meanc(xzdata)'; */      /* a row vector 1 by kxz */
    elseif (cols(xbar) eq 1) and (xbar eq 0); 
        mebar  = zeros(1,kxz);    /* a row vector 1 by kxz */
    else;
        mebar  = selif(xbar,index);
endif;

{MEtrue,Me1true} = ME_true(par_true,dgp,mebar,seldata,No,xv,znv,zpv,kxz,kx,kzn,kzp,cJ,inf,dummy,meNo);

save dgp, meNo, overlap, inf, dummy, toler, bound, MEtrue, Me1true, No, rep, cJ, par_true, xzdata;

failure      = zeros(rep,3);
ret          = zeros(rep,3);
est_time     = zeros(rep,3);
LL           = zeros(rep,3);
Hit          = zeros(rep,3);
sfreq        = zeros(rep,cJ);

b_op         = zeros(rep,kxz+cJ-1);
bse_op       = zeros(rep,kxz+cJ-1);
/*V_op         = arrayinit(rep|kxz+cJ-1|kxz+cJ-1,0);*/
ME1_op       = arrayinit(rep|kxz|cJ,0);
ME1_opse     = arrayinit(rep|kxz|cJ,0);
CPme1_op     = arrayinit(rep|kxz|cJ,0);
PRMSE_op     = zeros(rep,cJ);

if dgp eq 1;

b_nop        = zeros(rep,kxz+2+kxz+inf-2+kxz+cJ-inf-1);
bse_nop      = zeros(rep,kxz+2+kxz+inf-2+kxz+cJ-inf-1);
/*V_nop        = arrayinit(rep|kxz+2+kxz+inf-2+kxz+cJ-inf-1|kxz+2+kxz+inf-2+kxz+cJ-inf-1,0);*/
b_tiop       = zeros(rep,kxz+2+kxz+inf-1+kxz+cJ-inf);
bse_tiop     = zeros(rep,kxz+2+kxz+inf-1+kxz+cJ-inf);
/*V_tiop       = arrayinit(rep|kxz+2+kxz+inf-1+kxz+cJ-inf|kxz+2+kxz+inf-1+kxz+cJ-inf,0);*/

    else;

b_nop        = zeros(rep,kx+2+kzn+inf-2+kzp+cJ-inf-1);
bse_nop      = zeros(rep,kx+2+kzn+inf-2+kzp+cJ-inf-1);
/*V_nop        = arrayinit(rep|kx+2+kzn+inf-2+kzp+cJ-inf-1|kx+2+kzn+inf-2+kzp+cJ-inf-1,0);*/
b_tiop       = zeros(rep,kx+2+kzn+inf-1+kzp+cJ-inf);
bse_tiop     = zeros(rep,kx+2+kzn+inf-1+kzp+cJ-inf);
/*V_tiop       = arrayinit(rep|kx+2+kzn+inf-1+kzp+cJ-inf|kx+2+kzn+inf-1+kzp+cJ-inf,0);*/

endif;

ME1_nop      = arrayinit(rep|kxz|cJ,0);
ME1_nopse    = arrayinit(rep|kxz|cJ,0);
CPme1_nop    = arrayinit(rep|kxz|cJ,0);
PRMSE_nop    = zeros(rep,cJ);

ME1_tiop     = arrayinit(rep|kxz|cJ,0);
ME1_tiopse   = arrayinit(rep|kxz|cJ,0);
CPme1_tiop   = arrayinit(rep|kxz|cJ,0);
PRMSE_tiop   = zeros(rep,cJ);

if meNo gt 0;

me_op        = arrayinit(rep|meNo|kxz|cJ,0);
me_opse      = arrayinit(rep|meNo|kxz|cJ,0);
/*CPme_op      = arrayinit(rep|meNo|kxz|cJ,0);*/
me_nop       = arrayinit(rep|meNo|kxz|cJ,0);
me_nopse     = arrayinit(rep|meNo|kxz|cJ,0);
/*CPme_nop     = arrayinit(rep|meNo|kxz|cJ,0);*/
me_tiop      = arrayinit(rep|meNo|kxz|cJ,0);
me_tiopse    = arrayinit(rep|meNo|kxz|cJ,0);
/*CPme_tiop    = arrayinit(rep|meNo|kxz|cJ,0);*/

    else;

me_op        = arrayinit(rep|1|kxz|cJ,0);
me_opse      = arrayinit(rep|1|kxz|cJ,0);
/*CPme_op      = arrayinit(rep|1|kxz|cJ,0);*/
me_nop       = arrayinit(rep|1|kxz|cJ,0);
me_nopse     = arrayinit(rep|1|kxz|cJ,0);
/*CPme_nop     = arrayinit(rep|1|kxz|cJ,0);*/
me_tiop      = arrayinit(rep|1|kxz|cJ,0);
me_tiopse    = arrayinit(rep|1|kxz|cJ,0);
/*CPme_tiop    = arrayinit(rep|1|kxz|cJ,0);*/


endif;

Vuong_21     = zeros(rep,1);
Vuong_31     = zeros(rep,1);


/*@@@@@@@@@@@@@@@@@@@  begin MONTE-CARLO loop @@@@@@@@@@@@@@@@@@@@@@@@*/

for jrep (1,rep,1);

print "Replication # ";;print jrep;

if dgp eq 3; /* TIOP */

eps  = rndn(No,1);
epsp = rndn(No,1);
epsn = rndn(No,1);

rs      = xdata*beta + eps;
r       = zeros(No,3);

r[.,1]  = rs .lt mu[1];
r[.,3]  = rs .ge mu[2];
r[.,2]  = 1 - r[.,1] - r[.,3];

ysn     = zndata*betan + epsn;
ysp     = zpdata*betap + epsp;

yn                = zeros(No,rows(mun)+1);
yn[.,1]           = ysn .lt mun[1];
yn[.,rows(mun)+1] = ysn .ge mun[rows(mun)];

if rows(mun) gt 1;
    for i (2,rows(mun),1);
        yn[.,i]   = (ysn .ge mun[i-1]) .and (ysn .lt mun[i]);
    endfor;
endif;

yp                = zeros(No,rows(mup)+1);
yp[.,1]           = ysp .lt mup[1];
yp[.,rows(mup)+1] = ysp .ge mup[rows(mup)];

if rows(mup) gt 1;
    for i (2,rows(mup),1);
        yp[.,i]   = (ysp .ge mup[i-1]) .and (ysp .lt mup[i]);
    endfor;
endif;

yn      = yn*seqa(-rows(mun),1,rows(mun)+1);
yp      = yp*seqa(0,1,rows(mup)+1);

yr      = yn~zeros(No,1)~yp;
y       = sumr(r .* yr); /* observed y */

elseif dgp eq 2; /* NOP */

eps  = rndn(No,1);
epsp = rndn(No,1);
epsn = rndn(No,1);

rs      = xdata*beta + eps;
r       = zeros(No,3);

r[.,1]  = rs .lt mu[1];
r[.,3]  = rs .ge mu[2];
r[.,2]  = 1 - r[.,1] - r[.,3];

ysn     = zndata*betan + epsn;
ysp     = zpdata*betap + epsp;

yn                = zeros(No,rows(mun)+1);
yn[.,1]           = ysn .lt mun[1];
yn[.,rows(mun)+1] = ysn .ge mun[rows(mun)];

if rows(mun) gt 1;
    for i (2,rows(mun),1);
        yn[.,i]   = (ysn .ge mun[i-1]) .and (ysn .lt mun[i]);
    endfor;
endif;

yp                = zeros(No,rows(mup)+1);
yp[.,1]           = ysp .lt mup[1];
yp[.,rows(mup)+1] = ysp .ge mup[rows(mup)];

if rows(mup) gt 1;
    for i (2,rows(mup),1);
        yp[.,i]   = (ysp .ge mup[i-1]) .and (ysp .lt mup[i]);
    endfor;
endif;

yn      = yn*seqa(-rows(mun)-1,1,rows(mun)+1);
yp      = yp*seqa(1,1,rows(mup)+1);

yr      = yn~zeros(No,1)~yp;
y       = sumr(r .* yr); /* observed y */

elseif dgp eq 1; /* OP */

eps  = rndn(No,1);

ys      = xdata*beta + eps;

y       = zeros(No,cJ);
y[.,1]  = ys .lt mu[1];
y[.,cJ] = ys .ge mu[cJ-1];

if cJ gt 2;
    for i (2,cJ-1,1);
        y[.,i]   = (ys .ge mu[i-1]) .and (ys .lt mu[i]);
    endfor;
endif;

y       = y*seqa(1-inf,1,cJ); /* observed y */

endif;

/* recoding y as 1,2, ... , cJ-1, cJ */

yo  = y;              /* to keep original values of y */
y   = y+inf;
yxz = y~xzdata;

    {op,   opse,    ME1op,   ME1opse,   CPme1op,   MEop,   MEopse,    maxfop,   retop,   op_et,   PRMSEop,   Hitop, 
    nop,   nopse,   ME1nop,  ME1nopse,  CPme1nop,  MEnop,  MEnopse,   maxfnop,  retnop,  nop_et,  PRMSEnop,  Hitnop,
    tiop,  tiopse,  ME1tiop, ME1tiopse, CPme1tiop, MEtiop, MEtiopse,  maxftiop, rettiop, tiop_et, PRMSEtiop, Hittiop,
    Vuong21, Vuong31, freq, fail, fail_n, fail_t}
    = MC(No,dgp,sv_op,sv_nop,sv_tiop,xv,zpv,znv,
         y,seldata,inf,vsop,cJ,xzvars,xvars,zpvars,znvars,cv,par_true,mebar,ME1true,MEtrue,dummy,meNo,toler,bound);

sfreq[jrep,.]        = freq';
/*yv[.,jrep]           = y;*/

b_op[jrep,.]         = op';
bse_op[jrep,.]       = opse';
/*V_op[jrep,.,.]       = covop;*/
ME1_op[jrep,.,.]     = ME1op;
ME1_opse[jrep,.,.]   = ME1opse;
CPME1_op[jrep,.,.]   = CPME1op;
ME_op[jrep,.,.,.]    = MEop;
ME_opse[jrep,.,.,.]  = MEopse;
/*CPME_op[jrep,.,.,.]  = CPMEop;*/
/*p_op[jrep,.,.]       = popb;*/
PRMSE_op[jrep,.]     = PRMSEop';

b_nop[jrep,.]           = nop';
bse_nop[jrep,.]         = nopse';
/*V_nop[jrep,.,.]         = covnop;*/
ME1_nop[jrep,.,.]       = ME1nop;
ME1_nopse[jrep,.,.]     = ME1nopse;
CPME1_nop[jrep,.,.]     = CPME1nop;
ME_nop[jrep,.,.,.]      = MEnop;
ME_nopse[jrep,.,.,.]    = MEnopse;
/*CPME_nop[jrep,.,.,.]    = CPMEnop;*/
/*p_nop[jrep,.,.]         = pnop;*/
PRMSE_nop[jrep,.]       = PRMSEnop';

b_tiop[jrep,.]          = tiop';
bse_tiop[jrep,.]        = tiopse';
/*V_tiop[jrep,.,.]        = covtiop ;*/
ME1_tiop[jrep,.,.]      = ME1tiop;
ME1_tiopse[jrep,.,.]    = ME1tiopse;
CPME1_tiop[jrep,.,.]    = CPME1tiop;
ME_tiop[jrep,.,.,.]     = MEtiop;
ME_tiopse[jrep,.,.,.]   = MEtiopse;
/*CPME_tiop[jrep,.,.,.]   = CPMEtiop;*/
/*p_tiop[jrep,.,.]        = ptiop ;*/
PRMSE_tiop[jrep,.]      = PRMSEtiop';

ret[jrep,.]        = retop~retnop~rettiop;
failure[jrep,.]    = fail~fail_n~fail_t;
est_time[jrep,.]   = op_et~nop_et~tiop_et;
LL[jrep,.]         = maxfop~maxfnop~maxftiop;
/*LL0[jrep,.]        = maxfop0~maxfnop0~maxftiop0;*/
Hit[jrep,.]        = Hitop~Hitnop~Hittiop;

Vuong_21[jrep]     = Vuong21;
Vuong_31[jrep]     = Vuong31;

endfor; 

save sfreq, ret, failure, est_time, LL, Hit, Vuong_21, Vuong_31,
b_op,  bse_op,  ME1_op,  ME1_opse,  CPme1_op,  me_op,  me_opse,   PRMSE_op, 
b_nop, bse_nop, ME1_nop, ME1_nopse, CPme1_nop, me_nop, me_nopse,  PRMSE_nop,
b_tiop,bse_tiop,ME1_tiop,ME1_tiopse,CPme1_tiop,me_tiop,me_tiopse, PRMSE_tiop
;


/*********** The end of Monte Carlo loop *************/

/* to save the data set used in the latest MC iteration */

    let yxzvars = y x1 x2 x3;
    if not saved(yxz,"MC3_data",yxzvars);
        errorlog "Disk full! Cannot save the data set in MC3n_main.prg. Program terminated"; stop;
        end;
    endif;

a_sfreq       = meanc(sfreq);
indeks        = ret .eq 0;
ind_op        = indeks[.,1];
ind_nop       = indeks[.,2];
ind_tiop      = indeks[.,3];
ind           = (ind_op + ind_nop + ind_tiop) .eq 3;
NR            = sumc(ind);

me1_op_       = arrayinit(NR|kxz|cJ,0);
me1_nop_      = arrayinit(NR|kxz|cJ,0);
me1_tiop_     = arrayinit(NR|kxz|cJ,0);

j = 0;
for i (1,rep,1);

    if ind[i] eq 1;

        j = j + 1;

        me1_op_[j,.,.]       = me1_op[i,.,.];
        me1_nop_[j,.,.]      = me1_nop[i,.,.];
        me1_tiop_[j,.,.]     = me1_tiop[i,.,.];

    endif;
endfor;

clear me1_op, me1_nop, me1_tiop;

me1_opse_     = arrayinit(NR|kxz|cJ,0);
me1_nopse_    = arrayinit(NR|kxz|cJ,0);
me1_tiopse_   = arrayinit(NR|kxz|cJ,0);

j = 0;
for i (1,rep,1);

    if ind[i] eq 1;

        j = j + 1;

        me1_opse_[j,.,.]     = real(me1_opse[i,.,.]);
        me1_nopse_[j,.,.]    = real(me1_nopse[i,.,.]);
        me1_tiopse_[j,.,.]   = real(me1_tiopse[i,.,.]);

    endif;
endfor;

clear me1_opse, me1_nopse, me1_tiopse;

CPme1_op_     = arrayinit(NR|kxz|cJ,0);
CPme1_nop_    = arrayinit(NR|kxz|cJ,0);
CPme1_tiop_   = arrayinit(NR|kxz|cJ,0);

j = 0;
for i (1,rep,1);

    if ind[i] eq 1;

        j = j + 1;

        CPme1_op_[j,.,.]     = CPme1_op[i,.,.];
        CPme1_nop_[j,.,.]    = CPme1_nop[i,.,.];
        CPme1_tiop_[j,.,.]   = CPme1_tiop[i,.,.];

    endif;
endfor;

clear CPme1_op, CPme1_nop, CPme1_tiop;

se_op         = zeros(kxz,cJ);
se_nop        = zeros(kxz,cJ);
se_tiop       = zeros(kxz,cJ);

ss_op         = zeros(kxz,cJ);
ss_nop        = zeros(kxz,cJ);
ss_tiop       = zeros(kxz,cJ);

a_me1_op      = arraytomat(amean(me1_op_,3));
a_me1_nop     = arraytomat(amean(me1_nop_,3));
a_me1_tiop    = arraytomat(amean(me1_tiop_,3));

for i (1,NR,1);

    se_op    = (arraytomat(me1_op_[i,.,.]) - ME1true).^2 + se_op;
    ss_op    = (arraytomat(me1_op_[i,.,.]) - a_me1_op).^2 + ss_op;

    se_nop   = (arraytomat(me1_nop_[i,.,.]) - ME1true).^2  + se_nop;
    ss_nop   = (arraytomat(me1_nop_[i,.,.]) - a_me1_nop).^2 + ss_nop;

    se_tiop  = (arraytomat(me1_tiop_[i,.,.]) - ME1true).^2   + se_tiop;
    ss_tiop  = (arraytomat(me1_tiop_[i,.,.]) - a_me1_tiop).^2 + ss_tiop;

endfor;

a_rmse_me1_op    = sqrt(se_op./NR);
a_rmse_me1_nop   = sqrt(se_nop./NR);
a_rmse_me1_tiop  = sqrt(se_tiop./NR);

std_me1_op       = sqrt(ss_op./NR);
std_me1_nop      = sqrt(ss_nop./NR);
std_me1_tiop     = sqrt(ss_tiop./NR);

a_me1_opse      = arraytomat(amean(me1_opse_,3));
a_CPme1_op      = arraytomat(amean(CPme1_op_,3));

a_me1_nopse     = arraytomat(amean(me1_nopse_,3));
a_CPme1_nop     = arraytomat(amean(CPme1_nop_,3));

a_me1_tiopse    = arraytomat(amean(me1_tiopse_,3));
a_CPme1_tiop    = arraytomat(amean(CPme1_tiop_,3));

op_ratio1       = a_me1_opse./std_me1_op;
nop_ratio1      = a_me1_nopse./std_me1_nop;
tiop_ratio1     = a_me1_tiopse./std_me1_tiop;

a_PRMSE_op     = meanc(selif(PRMSE_op,ind));
a_PRMSE_nop    = meanc(selif(PRMSE_nop,ind));
a_PRMSE_tiop   = meanc(selif(PRMSE_tiop,ind));

a_PRMSE        = a_PRMSE_op~a_PRMSE_nop~a_PRMSE_tiop;

if meNo gt 0;

    ME_op_       = arrayinit(NR|meNo|kxz|cJ,0);
    ME_nop_      = arrayinit(NR|meNo|kxz|cJ,0);
    ME_tiop_     = arrayinit(NR|meNo|kxz|cJ,0);

    j = 0;
    for i (1,rep,1);

        if ind[i] eq 1;

            j = j + 1;

            ME_op_[j,.,.,.]       = ME_op[i,.,.,.];
            ME_nop_[j,.,.,.]      = ME_nop[i,.,.,.];
            ME_tiop_[j,.,.,.]     = ME_tiop[i,.,.,.];

        endif;
    endfor;
 
    clear ME_op, ME_nop, ME_tiop;

    ME_opse_     = arrayinit(NR|meNo|kxz|cJ,0);
    ME_nopse_    = arrayinit(NR|meNo|kxz|cJ,0);
    ME_tiopse_   = arrayinit(NR|meNo|kxz|cJ,0);

    j = 0;
    for i (1,rep,1);

        if ind[i] eq 1;

            j = j + 1;

            ME_opse_[j,.,.,.]     = real(ME_opse[i,.,.,.]);
            ME_nopse_[j,.,.,.]    = real(ME_nopse[i,.,.,.]);
            ME_tiopse_[j,.,.,.]   = real(ME_tiopse[i,.,.,.]);

        endif;
    endfor;

    clear ME_opse, ME_nopse, ME_tiopse;

    a_me_op       = amean(me_op_,4);
    a_me_opse     = amean(me_opse_,4);

    aa_me_op      = arraytomat(amean(amean(me_op_,4),3));
    aa_me_opse    = arraytomat(amean(amean(me_opse_,4),3));

    a_me_nop      = amean(me_nop_,4);
    a_me_nopse    = amean(me_nopse_,4);

    aa_me_nop     = arraytomat(amean(amean(me_nop_,4),3));
    aa_me_nopse   = arraytomat(amean(amean(me_nopse_,4),3));

    a_me_tiop     = amean(me_tiop_,4);
    a_me_tiopse   = amean(me_tiopse_,4);

    aa_me_tiop    = arraytomat(amean(amean(me_tiop_,4),3));
    aa_me_tiopse  = arraytomat(amean(amean(me_tiopse_,4),3));

    rmse_me_op    = arrayinit(meNo|kxz|cJ,0);
    rmse_me_nop   = arrayinit(meNo|kxz|cJ,0);
    rmse_me_tiop  = arrayinit(meNo|kxz|cJ,0);

    std_me_op     = arrayinit(meNo|kxz|cJ,0);
    std_me_nop    = arrayinit(meNo|kxz|cJ,0);
    std_me_tiop   = arrayinit(meNo|kxz|cJ,0);

    op_ratio     = arrayinit(meNo|kxz|cJ,0);
    nop_ratio    = arrayinit(meNo|kxz|cJ,0);
    tiop_ratio   = arrayinit(meNo|kxz|cJ,0);

    CPme_op       = zeros(kxz,cJ);
    CPme_nop      = zeros(kxz,cJ);
    CPme_tiop     = zeros(kxz,cJ);

    for j (1,meNo,1);

    se_op         = zeros(kxz,cJ);
    se_nop        = zeros(kxz,cJ);
    se_tiop       = zeros(kxz,cJ);

    ss_op         = zeros(kxz,cJ);
    ss_nop        = zeros(kxz,cJ);
    ss_tiop       = zeros(kxz,cJ);

        for i (1,NR,1);

            se_op    = (arraytomat(me_op_[i,j,.,.]) - arraytomat(MEtrue[j,.,.])).^2 + se_op;
            ss_op    = (arraytomat(me_op_[i,j,.,.]) - arraytomat(a_me_op[1,j,.,.])).^2 + ss_op;

            se_nop   = (arraytomat(me_nop_[i,j,.,.]) - arraytomat(MEtrue[j,.,.])).^2  + se_nop;
            ss_nop   = (arraytomat(me_nop_[i,j,.,.]) - arraytomat(a_me_nop[1,j,.,.])).^2 + ss_nop;

            se_tiop  = (arraytomat(me_tiop_[i,j,.,.]) - arraytomat(MEtrue[j,.,.])).^2   + se_tiop;
            ss_tiop  = (arraytomat(me_tiop_[i,j,.,.]) - arraytomat(a_me_tiop[1,j,.,.])).^2 + ss_tiop;

            CPme_op    =  ( ( ( ( arraytomat(me_op_[i,j,.,.]) -    cv*arraytomat(me_opse_[i,j,.,.]) )    .le arraytomat(MEtrue[j,.,.]) ) +
                              ( ( arraytomat(me_op_[i,j,.,.]) +    cv*arraytomat(me_opse_[i,j,.,.]) )    .ge arraytomat(MEtrue[j,.,.]) ) ) .eq 2 ) + CPme_op;
            CPme_nop   =  ( ( ( ( arraytomat(me_nop_[i,j,.,.]) -   cv*arraytomat(me_nopse_[i,j,.,.]) )   .le arraytomat(MEtrue[j,.,.]) ) +
                              ( ( arraytomat(me_nop_[i,j,.,.]) +   cv*arraytomat(me_nopse_[i,j,.,.]) )   .ge arraytomat(MEtrue[j,.,.]) ) ) .eq 2 ) + CPme_nop;
            CPme_tiop  =  ( ( ( ( arraytomat(me_tiop_[i,j,.,.]) -  cv*arraytomat(me_tiopse_[i,j,.,.]) )  .le arraytomat(MEtrue[j,.,.]) ) +
                              ( ( arraytomat(me_tiop_[i,j,.,.]) +  cv*arraytomat(me_tiopse_[i,j,.,.]) )  .ge arraytomat(MEtrue[j,.,.]) ) ) .eq 2 ) + CPme_tiop;

        endfor;

        rmse_me_op[j,.,.]    = sqrt(se_op./NR);
        rmse_me_nop[j,.,.]   = sqrt(se_nop./NR);
        rmse_me_tiop[j,.,.]  = sqrt(se_tiop./NR);

        std_me_op[j,.,.]     = sqrt(ss_op./NR);
        std_me_nop[j,.,.]    = sqrt(ss_nop./NR);
        std_me_tiop[j,.,.]   = sqrt(ss_tiop./NR);

        op_ratio[j,.,.]     = arraytomat(a_me_opse[1,j,.,.])./arraytomat(std_me_op[j,.,.]);
        nop_ratio[j,.,.]    = arraytomat(a_me_nopse[1,j,.,.])./arraytomat(std_me_nop[j,.,.]);
        tiop_ratio[j,.,.]   = arraytomat(a_me_tiopse[1,j,.,.])./arraytomat(std_me_tiop[j,.,.]);

    endfor;

    a_CPme_op       = CPme_op./(NR*meNo);
    a_CPme_nop      = CPme_nop./(NR*meNo);
    a_CPme_tiop     = CPme_tiop./(NR*meNo);

    a_op_ratio      = arraytomat(amean(op_ratio,3));
    a_nop_ratio     = arraytomat(amean(nop_ratio,3));
    a_tiop_ratio    = arraytomat(amean(tiop_ratio,3));

    a_std_me_op     = arraytomat(amean(std_me_op,3));
    a_std_me_nop    = arraytomat(amean(std_me_nop,3));
    a_std_me_tiop   = arraytomat(amean(std_me_tiop,3));

    a_rmse_me_op    = arraytomat(amean(rmse_me_op,3));
    a_rmse_me_nop   = arraytomat(amean(rmse_me_nop,3));
    a_rmse_me_tiop  = arraytomat(amean(rmse_me_tiop,3));

endif;

a_time = meanc(selif(est_time,ind));
LL     = selif(LL,ind);
a_LL   = meanc(LL);

K      = cols(b_op)~cols(b_nop)~cols(b_tiop);

LR_21    = -2*(LL[.,1] - LL[.,2]);
LR_31    = -2*(LL[.,1] - LL[.,3]);
LR_32    = -2*(LL[.,2] - LL[.,3]);

a_LR_21  = meanc(LR_21);
a_LR_31  = meanc(LR_31);
a_LR_32  = meanc(LR_32);

a_Vuong_21 = meanc(selif(Vuong_21,ind));
a_Vuong_31 = meanc(selif(Vuong_31,ind));

df21      = cols(b_nop)   - cols(b_op);
df31      = cols(b_tiop)  - cols(b_op);
df32      = cols(b_tiop)  - cols(b_nop);

if df21 gt 0;

TestLR21    = meanc(LR_21 .gt cdfchii(0.95,df21))|meanc(LR_21 .le 0-2*toler*No)|meanc((LR_21 .lt cdfchii(0.95,df21)) .and (LR_21 .gt 0-2*toler*No));

endif;

if df31 gt 0;

TestLR31    = meanc(LR_31 .gt cdfchii(0.95,df31))|meanc(LR_31 .le 0-2*toler*No)|meanc((LR_31 .lt cdfchii(0.95,df31)) .and (LR_31 .gt 0-2*toler*No));

endif;

TestLR32    = meanc(LR_32 .gt cdfchii(0.95,df32))|meanc(LR_32 .le 0-2*toler*No)|meanc((LR_32 .lt cdfchii(0.95,df32)) .and (LR_32 .gt 0-2*toler*No));

TestVuong21 = meanc(selif(Vuong_21,ind) .gt cv)|meanc(abs(selif(Vuong_21,ind)) .le cv)|meanc(selif(Vuong_21,ind) .lt -cv); 
TestVuong31 = meanc(selif(Vuong_31,ind) .gt cv)|meanc(abs(selif(Vuong_31,ind)) .le cv)|meanc(selif(Vuong_31,ind) .lt -cv); 

AIC    = -2*LL + 2*K;  
BIC    = -2*LL + ln(No)*K;  
cAIC   = -2*LL + (1+ln(No))*K;  
AICcor =   AIC + 2*K.*(K+1)./(No-K-1);
HQIC   = -2*LL + 2*K*ln(ln(No));
/*McFa   = 1 - (LL - K)./LL0;
McF    = 1 - LL./LL0;
CU     = (1 - exp(LL0*2/No)./exp(LL*2/No))./(1 - exp(LL0*2/No));*/

a_AIC     = meanc(AIC);
a_BIC     = meanc(BIC);
a_cAIC    = meanc(cAIC);
a_AICcor  = meanc(AICcor);
a_HQIC    = meanc(HQIC);
/*a_McF     = meanc(McF);
a_McFa    = meanc(McFa);
a_CU      = meanc(CU);*/
a_Hit     = meanc(Hit);

AICind       = minindc(AIC');
s_AIC        = meanc(AICind .eq 1)|meanc(AICind .eq 2)|meanc(AICind .eq 3);

BICind       = minindc(BIC');
s_BIC        = meanc(BICind .eq 1)|meanc(BICind .eq 2)|meanc(BICind .eq 3);

cAICind      = minindc(cAIC');
s_cAIC       = meanc(cAICind .eq 1)|meanc(cAICind .eq 2)|meanc(cAICind .eq 3);

AICcorind    = minindc(AICcor');
s_AICcor     = meanc(AICcorind .eq 1)|meanc(AICcorind .eq 2)|meanc(AICcorind .eq 3);

HQICind      = minindc(HQIC');
s_HQIC       = meanc(HQICind .eq 1)|meanc(HQICind .eq 2)|meanc(HQICind .eq 3);

/*McFind       = maxindc(McF');
s_McF        = meanc(McFind .eq 1)|meanc(McFind .eq 2)|meanc(McFind .eq 3);

McFaind      = maxindc(McFa');
s_McFa       = meanc(McFaind .eq 1)|meanc(McFaind .eq 2)|meanc(McFaind .eq 3);

CUind        = maxindc(CU');
s_CU         = meanc(CUind .eq 1)|meanc(CUind .eq 2)|meanc(CUind .eq 3);*/

Hitind       = maxindc(Hit');
s_Hit        = meanc(Hitind .eq 1)|meanc(Hitind .eq 2)|meanc(Hitind .eq 3);

if dgp eq 1;

    b_op_        = selif(b_op,ind_op);
    bse_op_      = selif(bse_op,ind_op);
    a_b_op       = meanc(b_op_);
    a_bse_op     = meanc(bse_op_);
    std_b_op     = sqrt(meanc((b_op_-ones(rows(b_op_),1)*a_b_op').^2));
    a_rmse_b_op  = sqrt(meanc((b_op_-ones(rows(b_op_),1)*par_true').^2));

    cp_b_op      = zeros(rows(b_op_),rows(par_true));

    for i (1,rows(b_op_),1);

        cp_b_op[i,.] = (((b_op_[i,.] - cv*bse_op_[i,.]) .le par_true' ) +
                        ((b_op_[i,.] + cv*bse_op_[i,.]) .ge par_true' )) .eq 2;

    endfor;

    a_cp_b_op    = meanc(cp_b_op);
    save   a_b_op, a_rmse_b_op, a_cp_b_op, a_bse_op, std_b_op;

elseif dgp eq 2;

    b_nop_        = selif(b_nop,ind_nop);
    bse_nop_      = selif(bse_nop,ind_nop);
    a_b_nop       = meanc(b_nop_);
    a_bse_nop     = meanc(bse_nop_);
    std_b_nop     = sqrt(meanc((b_nop_-ones(rows(b_nop_),1)*a_b_nop').^2));
    a_rmse_b_nop  = sqrt(meanc((b_nop_-ones(rows(b_nop_),1)*par_true').^2));

    cp_b_nop      = zeros(rows(b_nop_),rows(par_true));

    for i (1,rows(b_nop_),1);

        cp_b_nop[i,.] = (((b_nop_[i,.] - cv*bse_nop_[i,.]) .le par_true' ) +
                         ((b_nop_[i,.] + cv*bse_nop_[i,.]) .ge par_true' )) .eq 2;

    endfor;

    a_cp_b_nop    = meanc(cp_b_nop);
    save   a_b_nop, a_rmse_b_nop, a_cp_b_nop, a_bse_nop, std_b_nop;

elseif dgp eq 3;

    b_tiop_        = selif(b_tiop,ind_tiop);
    bse_tiop_      = selif(bse_tiop,ind_tiop);
    a_b_tiop       = meanc(b_tiop_);
    a_bse_tiop     = meanc(bse_tiop_);
    std_b_tiop     = sqrt(meanc((b_tiop_-ones(rows(b_tiop_),1)*a_b_tiop').^2));
    a_rmse_b_tiop  = sqrt(meanc((b_tiop_-ones(rows(b_tiop_),1)*par_true').^2));

    cp_b_tiop      = zeros(rows(b_tiop_),rows(par_true));

    for i (1,rows(b_tiop_),1);

        cp_b_tiop[i,.] = (((b_tiop_[i,.] - cv*bse_tiop_[i,.]) .le par_true' ) +
                         ((b_tiop_[i,.] + cv*bse_tiop_[i,.]) .ge par_true' )) .eq 2;

    endfor;

    a_cp_b_tiop    = meanc(cp_b_tiop);
    save   a_b_tiop, a_rmse_b_tiop, a_cp_b_tiop, a_bse_tiop, std_b_tiop;

endif;

/****************************/
  
screen on;

print "================================================================================";
print "Monte-Carlo simulations";
print;
print "       Sample size      No of iterations       No of used iterations";
print No~rep~NR;
print;
print;
print "Overlap (1-no, 2-partial, 3-complete):";; overlap;
print;
print "Average time per MC replication in minutes";
t1 = ((timeutc - t0)/60)/rep;     /* ML estimation time in minutes */
print t1;

if dgp eq 3;

    print;
    print "dgp is TIOP model";
    print;
    print "xv  =";; xv';
    print "znv =";; znv';
    print "zpv =";; zpv';
    print;
    print "       a_b_tiop          b_true          a_rmse_b_tiop       a_cp_b_tiop     a_bse_tiop/std_b_tiop";
    print a_b_tiop~par_true~a_rmse_b_tiop~a_cp_b_tiop~a_bse_tiop./std_b_tiop;
    print;
    print "       aa_rmse_b_tiop       aa_cp_b_tiop     aa_bse_tiop/std_b_tiop";
    print meanc(a_rmse_b_tiop)~meanc(a_cp_b_tiop)~meanc(a_bse_tiop./std_b_tiop);
    print;

    elseif dgp eq 2;

    print;
    print "dgp is NOP model";
    print;
    print "xv  =";; xv';
    print "znv =";; znv';
    print "zpv =";; zpv';
    print;
    print "       a_b_nop          b_true          a_rmse_b_nop          a_cp_b_nop     a_bse_nop/std_b_nop";
    print a_b_nop~par_true~a_rmse_b_nop~a_cp_b_nop~a_bse_nop./std_b_nop;
    print;
    print "       aa_rmse_b_nop       aa_cp_b_nop     aa_bse_nop/std_b_nop";
    print meanc(a_rmse_b_nop)~meanc(a_cp_b_nop)~meanc(a_bse_nop./std_b_nop);
    print;


    elseif dgp eq 1;

    print;
    print "dgp is OP model";
    print "xv  =";; xv';
    print;
    print "       a_b_op          b_true          a_rmse_b_op        a_cp_b_op     a_bse_op/std_b_op";
    print a_b_op~par_true~a_rmse_b_op~a_cp_b_op~a_bse_op./std_b_op;
    print;
    print "       aa_rmse_b_op       aa_cp_b_op     aa_bse_op/std_b_op";
    print meanc(a_rmse_b_op)~meanc(a_cp_b_op)~meanc(a_bse_op./std_b_op);
    print;

endif;

print;
print "Number of observations per parameter";
print "         OP              NOP            TIOP";  
print No/cols(b_op)~No/cols(b_nop)~No/cols(b_tiop);
print cols(b_op)~cols(b_nop)~cols(b_tiop);
print;
print "Average RMSE of Marginal Effects";
print "          OP              NOP              TIOP";  
print "ME: ";; print meanc(meanc(a_rmse_me_op))~meanc(meanc(a_rmse_me_nop))~meanc(meanc(a_rmse_me_tiop));
print;
Nomin_op   = 0.95*ones(kxz,cJ); Nomin_nop  = 0.95*ones(kxz,cJ); Nomin_tiop = 0.95*ones(kxz,cJ);
if dgp ne 1;
    if overlap eq 1;
        Nomin_nop   = 0.95*ones(1,cJ)|(0.95*ones(1,inf-1)~ones(1,cJ-inf+1))|(ones(1,inf)~0.95*ones(1,cJ-inf));
        Nomin_tiop  = 0.95*ones(1,cJ)|(0.95*ones(1,inf)~ones(1,cJ-inf))|(ones(1,inf-1)~0.95*ones(1,cJ-inf+1));
        elseif overlap eq 2;
        Nomin_nop  = 0.95*ones(2,cJ)|(0.95*ones(1,inf-1)~1~ones(1,cJ-inf));
    endif;
endif;
print "Average Coverage Probabilities of MEs:";
print "         OP              NOP              TIOP";  
print "ME: ";; print meanc(delif(vecr(a_CPme_op),vecr(Nomin_op-0.95)))~meanc(delif(vecr(a_CPme_nop),vecr(Nomin_nop-0.95)))~meanc(delif(vecr(a_CPme_tiop),vecr(Nomin_tiop-0.95)));
print;
print "Average ratio of estimated st. errors of MEs and st. deviation of estimated MEs";
print "         OP              NOP              TIOP";  
print "ME  ";; print meanc(packr(vecr(a_op_ratio)))~meanc(packr(vecr(a_nop_ratio)))~meanc(packr(vecr(a_tiop_ratio)));
print;
print "        Model selection results by information criteria and measures of fit";
print "                OP              NOP         TIOP";  
print "AIC   ";; print s_AIC';
print "BIC   ";; print s_BIC';
print "cAIC  ";; print s_cAIC';
print "AICcor";; print s_AICcor';
print "HQIC  ";; print s_HQIC';
/*print "McF   ";; print s_McF';
print "McFa  ";; print s_McFa';
print "CU    ";; print s_CU';*/
print "Hi    ";; print s_Hit';
print;
print "                    Average RMSE of probabilities";
print "         OP              NOP              TIOP";  
print meanc(a_PRMSE)';
print;
print "       a_Hit_op        a_Hit_nop       a_Hit_tiop";
print a_Hit';
print;
print "            Problems with overal estimation:";
print "         OP              NOP              TIOP";  
print meanc(ret .ne 0)';
print;
print "Vuong non-nested test:";
print "         NOP           Indecisive            OP";
print TestVuong21';
print "Average Vuong test statistic:";
print a_Vuong_21;
print;
print "Vuong non-nested test:";
print "         TIOP           Indecisive           OP";
print TestVuong31';
print "Average Vuong test statistic:";
print a_Vuong_31;
print;
print "LR test:";
print "       Ha: TIOP        Indecisive        Ho: NOP";
print TestLR32';
print "Average LR test statistic:";
print a_LR_32;
print;
print "Average RMSE of Marginal Effects";
print "             OP              NOP            TIOP";  
print "ME1:";; print meanc(meanc(a_rmse_me1_op))~meanc(meanc(a_rmse_me1_nop))~meanc(meanc(a_rmse_me1_tiop));
print;
print "Average Coverage Probabilities of MEs:";
print "             OP              NOP            TIOP";  
print "ME1:";; print meanc(delif(vecr(a_CPme1_op),vecr(Nomin_op-0.95)))~meanc(delif(vecr(a_CPme1_nop),vecr(Nomin_nop-0.95)))~meanc(delif(vecr(a_CPme1_tiop),vecr(Nomin_tiop-0.95)));
print;
print "Average ratio of estimated st. errors of MEs and st. deviation of estimated MEs";
print "             OP              NOP            TIOP";  
print "ME1 ";; print meanc(packr(vecr(op_ratio1)))~meanc(packr(vecr(nop_ratio1)))~meanc(packr(vecr(tiop_ratio1)));
print;
print "      a_time_op        a_time_nop       a_time_tiop";
print a_time';
print;
print "                       Problems with convergence:";
print "         OP              NOP            TIOP";  
print meanc(ret .lt 0 .or ret .gt 0.5)';
print;
print "                Problems with covariance of parameters:";
print "         OP              NOP            TIOP";  
print meanc(ret .eq 0.5)';
print;
print "                   Problems with covariance of MEs:";
print "         OP              NOP            TIOP";  
print meanc(ret .eq 0.25)';
print;
print "Coverage Probabilities of MEs: Average deviation from nominal level";
print "             OP              NOP            TIOP";  
print "ME: ";; print meanc(meanc(abs(a_CPme_op-Nomin_op)))~meanc(meanc(abs(a_CPme_nop-Nomin_nop)))~meanc(meanc(abs(a_CPme_tiop-Nomin_tiop)));
print "ME1:";; print meanc(meanc(abs(a_CPme1_op-Nomin_op)))~meanc(meanc(abs(a_CPme1_nop-Nomin_nop)))~meanc(meanc(abs(a_CPme1_tiop-Nomin_tiop)));
print;
print "     a_PRMSE_op      a_PRMSE_nop      a_PRMSE_tiop";
print a_PRMSE;
print;

if df21 gt 0;
print "LR test:";
print "       Ha: NOP         Indecisive       Ho: OP";
print TestLR21';
print "Average LR test statistic:";
print a_LR_21;
print;
endif;
print "LR test:";
print "       Ha: TIOP      Indecisive         Ho:  NOP";
print TestLR32';
print "Average LR test statistic:";
print a_LR_32;
print;
print "LR test:";
print "       Ha: TIOP        Indecisive       Ho: OP";
print TestLR31';
print "Average LR test statistic:";
print a_LR_31;
print;
print "a_rmse_me_op";
print a_rmse_me_op;
print;
print "a_rmse_me_nop";
print a_rmse_me_nop;
print;
print "a_rmse_me_tiop";
print a_rmse_me_tiop;
print;
print "a_rmse_me1_op";
print a_rmse_me1_op;
print;
print "a_rmse_me1_nop";
print a_rmse_me1_nop;
print;
print "a_rmse_me1_tiop";
print a_rmse_me1_tiop;
print;
print "a_CPme_op";
print a_CPme_op;
print;
print "a_CPme_nop";
print a_CPme_nop;
print;
print "a_CPme_tiop";
print a_CPme_tiop;
print;
print "a_CPme1_op";
print a_CPme1_op;
print;
print "a_CPme1_nop";
print a_CPme1_nop;
print;
print "a_CPme1_tiop";
print a_CPme1_tiop;
print;
print "a_me1_op";
print a_me1_op;
print;
print "a_me1_nop";
print a_me1_nop;
print;
print "a_me1_tiop";
print a_me1_tiop;
print;
print "a_me1_opse";
print a_me1_opse;
print;
print "a_me1_nopse";
print a_me1_nopse;
print;
print "a_me1_tiopse";
print a_me1_tiopse;
print;
print "ME1true";
print ME1true;
print;
print "      op_ratio1";
print op_ratio1;
print;
print "      nop_ratio1";
print nop_ratio1;
print;
print "      tiop_ratio1";
print tiop_ratio1;
print;
print "      op_ratio";
print a_op_ratio;
print;
print "      nop_ratio";
print a_nop_ratio;
print;
print "      tiop_ratio";
print a_tiop_ratio;
print;
print "       a_LL_op         a_LL_nop         a_LL_tiop";
print a_LL';
print;
print "       a_AIC_op        a_AIC_nop        a_AIC_tiop";
print a_AIC';
print;
print "       a_BIC_op        a_BIC_nop        a_BIC_tiop";
print a_BIC';
print;
print "       a_CAIC_op       a_CAIC_nop       a_CAIC_tiop";
print a_CAIC';
print;
/*print "      a_McF_op         a_McF_nop        a_McF_tiop";
print a_McF';
print;
print "      a_McFa_op        a_McFa_nop       a_McFa_tiop";
print a_McFa';
print;
print "       a_CU_op          a_CU_nop         a_CU_tiop";
print a_CU';
print;*/
number = seqa(1,1,rep);
print "ret";
print number~ret;
print;
print "        Iteration        fail           fail_n          fail_t";
print number~failure;
print;

    /* saveall H:\Andrey\MC3_op_500_1000; */ /* Gauss shuts down after this command */