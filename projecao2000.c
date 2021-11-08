/* -------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* -------------------------------------------------------- */
/* // Date     : 04/03/95 10:20am                           */
/* // Code By  : Alvaro G. S. Barcellos                     */
/* -------------------------------------------------------- */
/* --------------------------------------------------------
roteiro de utilizacao das funcoes:

       1) definir o elipsoide: (SEMPRE....)
           int elipsoid (int k, double cm, double bl)
               k = numero do elipsoide
                   CLA_1866 = 0, HAY_1924 = 1, SAD_1969 = 2, WGS_1984 = 3
               cm = meridiano central (graus e decimos  W < 0)
               bl = latitude base     (graus e decimos  S < 0)

       2) converter geografica para utm
           void geo2utm ( double xdeg, double ydeg, double *xkm, double *ykm)
               xdeg = longitude (graus e decimos  W < 0)
               ydeg = latitude  (graus e decimos  S < 0)
               xkm  = este  UTM (metros)
               ykm  = norte UTM (metros)

       3) converter utm para geografica
           void utm2geo (double xkm, double ykm, double * xdeg, double * ydeg)
               xkm  = este  UTM (metros)
               ykm  = norte UTM (metros)
               xdeg = longitude (graus e decimos  W < 0)
               ydeg = latitude  (graus e decimos  S < 0)

       4) calculo da convergencia do meridiano em relacao ao meridiano central
           void lm2cm (double lon, double lat, double mc, double * ecm)
               lon = longitude do ponto (graus e decimos  W < 0)
               lat = latitude do ponto  (graus e decimos  S < 0)
               mc  = meridiano central  (graus e decimos  W < 0)
               ecm = convergencia meridiana

       5) converte distancias Utm para distancias locais (met) e vice-versa
           void utm2mtu (double *dx, double *dy, double ecm)
               dx = coordenada este   (metros)
               dy = coordenada norte  (metros)
               ecm = convergencia meridiana

       6) converter geografica para metrica policonica
           void geo2pol ( double xdeg, double ydeg, double *xkm, double *ykm)
               xdeg = longitude (graus e decimos  W < 0)
               ydeg = latitude  (graus e decimos  S < 0)
               xkm  = este  UTM (metros)
               ykm  = norte UTM (metros)

       7) converter metrica policonica para geografica
           void pol2geo (double xkm, double ykm, double * xdeg, double * ydeg)
               xkm  = este  UTM (metros)
               ykm  = norte UTM (metros)
               xdeg = longitude (graus e decimos  W < 0)
               ydeg = latitude  (graus e decimos  S < 0)

referencias :

   1) elipsoid(...), geo2utm(...), utm2geo(...)

      J.P.Snyder, 1982, Map projections used by the USGS;
      USGS Bulletin no. 1532, p 129-131, p 256-258.

      valores corretos ate metro e validos apenas ate 4 graus do meridiano central

   2) lm2cm(...), utm2mtu(...)

       manual de geodesia ...

   3) phi4(...), adjlon(...), pol2geo(...), geo2pol(...)

      modificadas a partir de:

      GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.0
      U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/18/92


*/
/**********************************************************/
/* -------------------------------------------------------- */
#define M_PI  3.14159265358979323846

/* parametros de elipsoide */

static double a, f, e, esq, epsq;
#define k0 0.9996

/* constantes de conversao */

const double deg2rad = (M_PI/180.0);
const double rad2deg = (180.0/M_PI);

const double tol = 0.0000001;

/* parametros comuns */

static double lam0, phi0;
static double xp0, yp0, xu0, yu0;
static double c0, c2, c4, c6, m0;
static double w0, w1, w2, w3, w4, w5, w6;
static double t, m, n, c;

/* elipsoides disponiveis */

enum  { CLA_1866 = 0, HAY_1924, SAD_1969, WGS_1984, GRS_1980 };

/* ---------------------------------------------------------------------- */
double mfcn (double phi)
{
    phi = a * (c0*phi-c2*sin(2.0*phi)+c4*sin(4.0*phi)-c6*sin(6.0*phi));
    return (phi);
}

/* ---------------------------------------------------------------------- */
void etmset (double cm, double bl)
{
    lam0 = cm * deg2rad;
    phi0 = bl * deg2rad;
    w2 = e * e;
    w4 = w2 * w2;
    w6 = w2 * w4;
    c0 =  1.0 - w2/4.0 - 3.0 * w4/64.0 - 5.0 * w6/256.0;
    c2 =  3.0 * w2/8.0 + 3.0 * w4/32.0 + 45.0 * w6/1024.0;
    c4 = 15.0 * w4/256.0 + 45.0 * w6/1024.0;
    c6 = 35.0 * w6/3072.0;
    m0 = mfcn (phi0);

/* metrica origem UTM */
    xu0 =   500000.0;
    yu0 = 10000000.0;

/* metrica origem Policonica */

    xp0 = 10000000.0;
    yp0 = 10000000.0;
}

/* ---------------------------------------------------------------------- */
/*  define os parametros em relacao ao ellipsoide */
int elipsoid (int k, double cm, double bl)
{
int err;
    err = 0;
    switch ( k) {
        case CLA_1866 :   /* 0 */
        /* parametros do elipsoide de clarke 1866 */
                a = 6378.206;
                f = 1.0/294.979;
                break;
        case HAY_1924 :   /* 1 */
        /* parametros do elipsoide HAYFORD 1924  */
                a = 6378.388;
                f = 1.0/297.0;
                break;
        case SAD_1969:   /* 2 */
        /* parametros do elipsoide SAD-69  */
                a = 6378.160;
                f = 1.0/298.25;
                break;
        case WGS_1984 :  /* 3 */
        /* parametros do elipsoide WGS-84 */
                a = 6378.137;
                f = 1.0/298.257223563;
                break;
        case GRS_1980 :  /* 4 */
        /* parametros do elipsoide GRS-80 */
                a = 6378.137;
                f = 1.0/298.257222101;
                break;
        default :
                err = 1;
                break;
    }
    if (err) return (err);

    esq =  (f * ( 2.0 - f));
    e = sqrt(esq);
    epsq = esq / (1.0 - esq);
    etmset(cm,bl);
    return (0);
}

/* ---------------------------------------------------------------------- */
void geo2utm ( double xdeg, double ydeg, double *xkm, double *ykm)
{
    xdeg = xdeg * deg2rad;
    ydeg = ydeg * deg2rad;
    w0 = tan(ydeg);
    t = w0 * w0;
    w0 = sin(ydeg);
    n =  a / sqrt (1.0 - esq * w0 * w0);
    w0 = cos(ydeg);
    c = epsq * w0 * w0 ;
    w0 = w0 * (xdeg - lam0);  /* biga powers */
    w2 = w0 * w0;
    w3 = w0 * w2;
    w4 = w2 * w2;
    w5 = w0 * w4;
    w6 = w2 * w4;
    m = mfcn(ydeg) - m0;

    *xkm = k0 * n * (w0 + (1.0 - t + c) * w3 / 6.0 +
           (5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * epsq) * w5 / 120.0);
    *ykm = k0 * ( m + n * tan(ydeg) * ( w2 / 2.0 +
           (5.0 - t + 9.0 * c + 4 * c * c) * w4 / 24.0) +
           (61.0 - 58.0 * t + t * t  + 600.0 * c - 330.0 * epsq) * w6 / 720.0);

    /* metrica normal */
    {
        *xkm =   xu0 + *xkm * 1000.0;
        *ykm =   yu0 + *ykm * 1000.0;
    }

}

/* ---------------------------------------------------------------------- */
void utm2geo (double xkm, double ykm, double * xdeg, double * ydeg)
{
    /* metrica normal */
    {
        xkm = (xkm - xu0) / 1000.0;
        ykm = (ykm - yu0) / 1000.0;
    }

    /* e1 */
    w0 = sqrt (1.0 - esq);
    w1 = (1.0 - w0) / (1.0 + w0);
    w2 = w1 * w1;
    w3 = w1 * w2;
    w4 = w2 * w2;

    /* mu */
    m = m0 + ykm / k0;
    m = m / (a * (1.0 - esq / 4.0 - 3.0 * esq * esq / 64.0 -
             5.0 * esq * esq * esq / 256.0));

    /* phi */
    m = m + (3.0 * w1 / 2.0 - 27.0 * w3 / 32.0) * sin (2.0 * m) +
                    (21.0 * w2 / 16.0 - 55.0 * w4 / 32.0) * sin (4.0 * m) +
                    (151.0 * w3 / 96.0) * sin (6.0 * m);

    w0 = tan(m);
    t = w0 * w0;

    w0 = cos(m);
    c = esq * w0 * w0;

    w0 = sin(m);
    w1 = 1.0 - esq * w0 * w0;
    n = a / sqrt (w1);

    /* r */
    w0 = a * (1.0 - esq) / (w1 * sqrt (w1));

    /* d */
    w1 = xkm / (n * k0);
    w2 = w1 * w1;
    w3 = w1 * w2;
    w4 = w2 * w2;
    w5 = w2 * w3;

    ykm = m - (n * tan (m) / w0) * (w2/2.0
            - (5.0 + 3.0*t + 10.0*c - 4.0*c*c - 9.0*epsq)
            * w4 / 24.0
            + (61.0 + 90.0*t + 298.0*c + 45.0*t*t -  252.0*epsq - 3.0*c*c)
            * w6/720.0);

    xkm = lam0 + (w1 - (1.0 + 2.0*t + c) * w3/6.0
               + (5.0 - 2.0*c + 28.0*t - 3.0*c*c
               + 8.0*epsq + 24.0*t*t) * w5/120.0)
               / cos(m);

    *xdeg = xkm * rad2deg;
    *ydeg = ykm * rad2deg;
}

/* ---------------------------------------------------------------------- */
double hms2seg (double v)
{
   double h, m, s, d;
   d = 1;
   if (v < 0) d = -1.0;
   v = v * d;
   h = (int) v;
   m = (int) ((v - h) * 60);
   s = (((v - h) * 60) - m) * 60;
   v = d * (h*10000+m*100+s);
   return (v);
}

/* ---------------------------------------------------------------------- */
/* calculo da convergencia do meridiano em relacao ao meridiano central */
int lm2cm (double lon, double lat, double mc, double * ecm)
{
double p0, p1, p3, p5, X1, X3, X5;
double tg, sn, cn, ec, s1, cn2;

       p0 =  mc - lon;
       p1 = fabs(p0*3600.0) / 10000.0;
       p3 = p1 * p1 * p1;
       p5 = p3 * p1 * p1;

       ec = fabs(lat);
       sn = sin(ec*deg2rad);
       cn = cos(ec*deg2rad);
       tg = sn/cn;
       cn2 = cn * cn;
       ec =  w2 * cn2;

       s1 = sin (1.0/3600.0 * deg2rad);
       s1 = s1 * s1;

       X1 = sn * 1.0E+04;
       X3 = s1*sn*cn2/3.0 * ( 1 + 3*ec + 2*ec*ec ) * 1.0E+12;
       X5 = s1*s1*sn*cn2*cn2/15.0 * (2.0 - tg*tg)  * 1.0E+20;

       ec = p1 * X1 + p3 * X3 + p5 * X5;

       if (p0 < 0)  ec = ec * -1.0;
       ec = ec / 3600.0;

       *ecm = ec;
       return (0);
}

/* ---------------------------------------------------------------------- */
// converte distancias Utm para distancias locais (met) e vice-versa
void utm2mtu (double *dx, double *dy, double ecm)
{
double px, py, tg, ds;

       px = *dx;
       py = *dy;

       if (px == 0 && py == 0) {
           return ;
       }
       else if (px == 0) {
           if (py < 0) tg = 3.0/2.0*M_PI;
           if (py > 0) tg = 1.0/2.0*M_PI;
       }
       else if (py == 0) {
           if (px < 0) tg = M_PI;
           if (px > 0) tg = 2.0*M_PI;
       }
       else {
           tg = atan2 (px,py);
       }

       ds = sqrt (px*px + py*py);

       tg = tg + ecm*deg2rad;

       *dx = ds * sin (tg);
       *dy = ds * cos (tg);
}

/* -------------------------------------------------------- */
/*
C **********************************************************************
C ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.0 **
C ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/18/92 **
C **********************************************************************
*/
/**********************************************************/
double adjlon (double lon)

/*function to adjust longitude angle to modulo 180 degrees.*/

{
      double ss, pi2;
      ss = 1;
      if (lon < 0.0) ss = -1.0;
      pi2 = 2.0 * M_PI * ss;
      while (fabs(lon) > M_PI) {
         lon = lon - pi2;
         }
      return (lon);
}

/* -------------------------------------------------------- */
int phi4 (double al, double bl, double *cl, double * phic)
/* function to compute latitude angle (phi-4) by iteration.*/
{
extern double c0,c2,c4,c6;
extern double w0,w1,w2,w3;
double tol, sinphi, tanphi;
double cc, ml, mlp, phi;
int it, nit;

      tol = 0.0000000000001;
      nit = 20;
      phi = al;

      for (it = 0; it < nit; it++) {
         sinphi = sin (phi);
         tanphi = tan (phi);
         cc = tanphi * sqrt (1.0 - esq * sinphi * sinphi);
         ml = c0 * phi - c2 * sin (2.0 * phi) + c4 * sin (4.0 * phi)
            - c6 * sin (6.0 * phi);
         mlp = c0 - 2.0 * c2 * cos (2.0 * phi) +
            + 4.0 * c4 * cos (4.0 * phi) - 6.0 * c6 * cos (6.0 * phi);

         w0 = 2.0 * ml + cc * (ml * ml + bl)
              - 2.0 * al * (cc * ml + 1.0) ;
         w1 = esq * sin (2.0 * phi) * (ml * ml + bl - 2.0 * al * ml)
              / (2.0 * cc);
         w2 = 2.0 * (al - ml) * (cc * mlp - 2.0 / sin (2.0 * phi))
              - 2.0 * mlp;
         w3 = w0 / (w1 + w2);

         phi = phi + w3;
         if (fabs(w3) < tol) break;
         }

      *phic = phi;
      *cl   = cc;

      return(0);
}

/*......................................................................
             policonic forward transformation  .
  ......................................................................*/
int geo2pol (double lon,double lat,double *xpol,double *ypol)
{
extern double c0, c2, c4, c6;
extern double phi0, lam0, a, e, esq;
double  xp, yp, con, ml, ms;
double  sinphi, cosphi;

         lat=(lat)*deg2rad;
         lon=(lon)*deg2rad;
         con = adjlon (lon - lam0);
         if ( fabs(lat) < tol) {
             xp = a * con;
             yp = a * phi0;
             }
         else {
             sinphi = sin (lat);
             cosphi = cos (lat);
             ml = c0 * lat - c2 * sin (2.0 * lat)
                           + c4 * sin (4.0 * lat)
                           - c6 * sin (6.0 * lat);
             ms = e * sinphi;
             ms = cosphi / sqrt (1.0 - ms * ms);
             con = con * sinphi;
             xp = a * ms * sin (con) / sinphi;
             yp = a * (ml - phi0 + ms * (1.0 - cos(con)) / sinphi);
             }

         *xpol = xp0 + xp * 1000.0;
         *ypol = yp0 + yp * 1000.0;

         return (0);
}
/*......................................................................
                       .  inverse transformation  .
  .....................................................................*/
int pol2geo (double xpol,double ypol,double *lon,double *lat)
{
extern double lam0, phi0, xp0, yp0, a;
double lt,ln,al,bl,cl;

         xpol = (xpol - xp0) / 1000.0;
         ypol = (ypol - yp0) / 1000.0;
         al = phi0 + ypol / a;
         if (fabs (al) < tol) {
            ln = xpol / a + lam0;
            lt = 0.0;
            }
         else {
            bl = al * al + (xpol / a) * (xpol / a);
            phi4 (al, bl, &cl, &lt);
            ln = adjlon (asin (xpol * cl / a) / sin (lt) + lam0);
            }

      *lat=(lt) * rad2deg;
      *lon=(ln) * rad2deg;

      return(0);
}

double deg2gms (double deg) {
	double  g, m, s;

	g = (int) deg;
	m = (int)((deg - g) * 60);
	s = (((deg - g) - m / 60.0) * 3600) ;

	return ( g * 10000.0 + m * 100.0 + s ) ;
	
}	
	
main (int argc, char * argv[])
{
 double cm, bl, lat, lon, xpo, ypo;
 int hh, md, gm, xy ;
 
 if (argc < 5) {
	 fprintf(stderr,
	" Use: projecao hh cm bl md gm < input > output\n"
	" hh -> elipsoide\n"
	" CLA_1866 = 0, HAY_1924 = 1, SAD_1969 = 2," 
	" WGS_1984 = 3, GRS_1980 (SIRGAS) = 4\n"
	" cm = central meridiano, bl = base latitude\n"
	" md -> universal mercator = 0, policonica = 1\n" 
	" gm -> graus decimais = 0, graus minuto segundo = 1\n" 
	" xy -> xyz = 0 yxz = 1\n" 
	);
		   
	 return (1);
	 }

 hh = atoi(argv[1]);
 cm = atof(argv[2]);
 bl = atof(argv[3]);
 md = atoi(argv[4]);
 gm = atoi(argv[5]);

 elipsoid (hh,cm,bl);

 {

/* conversao (Alvaro Barcellos 07/2013	) */
	 
	 double deg, gms, g, m, s;
	 double shot,bat;
	 char nil[1024];
	 char line[1024];
	 char buff[1024];
	 char ht, hg;

while (!feof(stdin) && fgets(buff,1024,stdin)) {
	 
	 if (xy == 0) sscanf(buff,"%lf %lf %lf",&xpo,&ypo,&bat);
	 if (xy == 1) sscanf(buff,"%lf %lf %lf",&ypo,&xpo,&bat);

	if (md == 0) utm2geo(xpo,ypo,&lon,&lat);

	if (md == 1) pol2geo(xpo+5000000,ypo,&lon,&lat);

	if (gm == 1) {

		lat = deg2gms(lat);
		lon = deg2gms(lon);

		ht = 'N';
		if (lat < 0.0) {
			lat *= -1;
			ht = 'S';
			}
		hg = 'E';
		if (lon < 0.0) {
			lon *= -1;
			hg = 'W';
			}
		sprintf (buff,"%09.2lf%c %010.2lf%c %06.2lf\n", line,nil,shot,lat,ht,lon,hg,xpo,ypo,bat);	 
	}
        else {
		sprintf (buff,"%09.6lf %09.6lf %06.2lf\n", lat, lon, bat);
	}

	 fputs(buff,stdout);
 	}

 }


			 
 /*
	 	gms = deg2gms(deg,&g,&m,&s);
    {
 lat = atof (argv[4]);
 lon = atof (argv[5]);
 geo2pol(lon,lat,&xpo,&ypo);
 printf ("F %lf %lf == %lf %lf\n",lat,lon,ypo,xpo);
 pol2geo(xpo,ypo,&lon,&lat);
 printf ("I %lf %lf == %lf %lf\n",lat,lon,ypo,xpo);
    }
 */
			 
 return(0);
}


