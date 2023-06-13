/*********************************************************************************
**                                                                              **
**	Program	: solvate.c	(c) 1996-2010 by Helmut Grubmuller, Gottingen           **
**                                                                              **
*********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#define PI 3.14159265358979323846
#define XPLOR_SCRIPT_NAME "mkpsf.inp"
#define XPLOR_TOPOLOGY_FILE_WATER "toph19.sol"
#define XPLOR_PARAMETER_FILE_WATER "param19.sol"
#define XPLOR_TOPOLOGY_FILE_IONS "toph19.nacl"
#define XPLOR_PARAMETER_FILE_IONS "param19.nacl"
#define MAX_STRING_SIZE 255
#define MAX_NR_OF_PDB_ATOMS 999999
#define MAX_NR_OF_PDB_GROUPS 9999
#define MAX_NR_OF_PDB_SEGMENTS 999
#define SURFACE_LEVEL 0.25
#define TIMAX_FACTOR 10
#define SIG_H2O 3.04166
#define BOND_LENGTH_H2O 0.9572  /* as in XPLOR TIP3W; otherwise problems with SHAKE */
#define BOND_ANGLE_H2O (105.0*PI/180.0)
#define INITIAL_STEP_WIDTH_FOR_CALC_BOUNDING_SPHERE 1.0
#define STEP_WIDTH_DECREASE_FOR_CALC_BOUNDING_SPHERE 1.25
#define MINIMAL_STEP_WIDTH_FOR_CALC_BOUNDING_SPHERE 0.000001
#define MAX_NR_OF_GRID_POINTS_FACTOR 1.2
#define ANGLE_STEP_FOR_CUT_AWAY 5.0
#define NR_OF_DICE_DATA 100000
#define FACTOR_FOR_INITIAL_LAMBDA_FOR_LEARNS_GAUSSIANS 2.0
#define ENTROPY_AVERAGE_LENGTH 1000.0
#define NOISE_FOR_LEARN_CODEBOOK_POSITIONS 0.1
#define RHOAVERSTART_FACTOR 2.0
#define RHOAVEREND_FACTOR (1.0/20.0)
#define EPS_CONST_FOR_LEARN_POSITIONS 0.02
#define EPS_CONST_FOR_LEARN_DIRECTIONS 0.01
#define EPS_CONST_FOR_LEARN_VARIANCES 0.01
#define EPS_START_FOR_LEARN_DIRECTIONS 0.01
#define EPS_END_FOR_LEARN_DIRECTIONS 0.0001
#define MU_FOR_LEARN_VARIANCES_INDIVIDUAL 0.1
#define RESOLUTION_FOR_GET_MIN_DISTANCE 1.0
#define SCALE_FACTOR_FOR_GET_MIN_DISTANCE 100.0
#define ACCURACY_OF_GET_MIN_DISTANCE 0.01
#define INITIAL_RHOFAC 5.0
#define INITIAL_RHOFAC_SHRINK_FACTOR 0.5
#define ACCURACY_OF_ADJUST_SURFACE_DISTANCE 0.05
#define ADD_LAYER_FOR_CREATE_VOLUME 10.0
#define INITIAL_STEPSIZE_FOR_SURFACE_STATISTICS 0.1
#define ACCURACY_FOR_SURFACE_STATISTICS 0.00001
#define DISTANCE_TIMES_WATER_SHELL_THICKNESS 2.0
#define BOUNDARY_SURFACE_THICKNESS_FACTOR 0.05
#define INITIAL_STEPSIZE_FOR_PLACE_WATER 0.5
#define FINAL_STEPSIZE_FOR_PLACE_WATER 0.001
#define DEFAULT_WATER_SHELL_THICKNESS 10.0
#define DEFAULT_BOUNDARY_RADIUS 100000.0
#define DEFAULT_NR_OF_GAUSSIANS 1
#define MAXNWATER_FACTOR 1.0
#define KT (1.380662e-23 * 300.0)
#define ANGSTROM_PER_METER 1.0e10
#define ELEMENTARY_CHARGE 1.6021892e-19
#define AVOGADRO 6.022045e23
#define ISOTONIC_CONCENTRATION 0.154     /* [mol/liter] */
#define LITER_PER_CUBIC_METER 1000.0
#define DEBYE_HUECKEL_LENGTH (ANGSTROM_PER_METER/sqrt(2.0*AVOGADRO*ISOTONIC_CONCENTRATION*LITER_PER_CUBIC_METER*ELEMENTARY_CHARGE*ELEMENTARY_CHARGE/8.85418782e-12/80.0/KT))    /* [A] */
#define ISOTONIC_WATER_PER_ION 361.0
#define FOUR_PI_EPS_EPSR (4.0*PI*8.85418782e-12*80.0)
#define NR_OF_DEBYE_HUECKEL_MONTE_CARLO_MOVES 2000000

#define PDB_P_ATOMNR            5
#define PDB_L_ATOMNR            6
#define PDB_P_ATOMSTR           12
#define PDB_L_ATOMSTR           4
#define PDB_P_RESIDSTR          17
#define PDB_L_RESIDSTR          4
#define PDB_P_RESIDNR           22
#define PDB_L_RESIDNR           5
#define PDB_P_XKO               30
#define PDB_L_XKO               8
#define PDB_P_YKO               38
#define PDB_L_YKO               8
#define PDB_P_ZKO               46
#define PDB_L_ZKO               8
#define PDB_P_BETA              54
#define PDB_L_BETA              6
#define PDB_P_HARMON            60
#define PDB_L_HARMON            6
#define PDB_P_SEGIDSTR          72
#define PDB_L_SEGIDSTR          4

#define PSF_P_ATOMNR            0
#define PSF_L_ATOMNR            8
#define PSF_P_SEGIDSTR          9
#define PSF_L_SEGIDSTR          4
#define PSF_P_RESIDNR           13
#define PSF_L_RESIDNR           6
#define PSF_P_RESIDSTR          19
#define PSF_L_RESIDSTR          4
#define PSF_P_ATOMSTR           24
#define PSF_L_ATOMSTR           4
#define PSF_P_ATOMTYPSTR        29
#define PSF_L_ATOMTYPSTR        4
#define PSF_P_CHARGE            34
#define PSF_L_CHARGE            15
#define PSF_P_MASS              50
#define PSF_L_MASS              9


/**************** functions ***************************************************/

static void error(int err_nr);
static void terminate(void);
static void here(int here_nr);
static void skip_line(FILE *fp);
static int read_string(FILE *f,char *str);
static int strup(char *str);
static int erase_spaces(char *str);

static int  read_next_pdbdata(FILE *fh,int *atom_nr,double *ko,double *beta,double *harmon,char *atomstr,char *residstr,int *residnr,char *segidstr);
static int  patch_pdb_coords(FILE *fh,FILE *fh1,double x,double y,double z);
static void add_pdb_tip(FILE *fp,int flag,double x,double y,double z,long int atom_nr,int res_id,int seg_id);
static void add_pdb_ion(FILE *fp,int iontype,double x,double y,double z,long int atom_nr,int res_id);
static int goto_psf_atomdata(FILE *h_psf);
static int read_next_psfatom(FILE *fh,int *atom_nr,char *atom_typestr,double *charge,double *mass,char *atomstr,char *residstr,int *residnr,char *segidstr);
static void read_pdb_file(char *infile);
static void read_atom_charges_from_psf_file(char *infile);
static void calc_bounding_sphere(void);
static void cut_away(double theta, double phi);
static void create_convex_volume(void);
static double calc_maxdist(double xr, double yr, double zr);
static void dice_data(void);
static void do_pca(void);
static void init_codebook_positions(void);
static void init_codebook_directions(void);
static void init_codebook_variances(void);
static void compute_activities(void);
static void learn_codebook_positions(double eps);
static void learn_codebook_directions(double eps);
static void learn_codebook_variances_soft_individual(double eps, double mu);
static void learn_codebook_variances_soft_global(double eps, double mu);
static void learn_codebook_variances_hard(double eps);
static void equalize_codebook_variances(void);
static void outputs(void);
static double get_rhoaver2(double rhoaverstart, double rhoaverend, long int timax, long int ti);
static double get_eps(double epsstart, double epsend, long int timax, long int ti);
static void learn_gaussians(void);
static double density_varrho(double x, double y, double z, double rhofac);
static double density(double x, double y, double z);
static double density_gradient(double x, double y, double z, double *gx, double *gy, double *gz);
static double density_grad_curve(double x, double y, double z, double *gx, double *gy, double *gz, double *g, double *cr, double *curvature);
static double get_min_distance(double rhofac, double level);
static void adjust_surface_distance(void);
static void create_volume(void);
static void save_gaussian_parameters(void);
static void save_boundary_parameters(void);
static void surface_statistics(void);
static void place_ions(void);
static void dice_hydrogens(double x0, double y0, double z0, double *x1, double *y1, double *z1, double *x2, double *y2, double *z2);
static double vdw_energy_force(double x, double y, double z, double *fx, double *fy, double *fz);
static void place_water_molecules(void);
static void check_max_atom_nr(long int nr);
static void check_max_segment_nr(int nr);
static void write_xplor_segment(FILE *fp, int segid, long int nr, int diff);
static void write_xplor_segment_na(FILE *fp, long int nr);
static void write_xplor_segment_cl(FILE *fp, long int nr);
static void write_water_molecules(char *infile, char *outfile);
static void group_water_molecules(void);

void  *malloc(size_t size);
void  *realloc(void* ptr, size_t size);
void   free(void *ptr);
	
static int  	*ivector(int nl, int nh);
static void	free_ivector(int *v, int nl, int nh);

static double	*dvector(int nl, int nh);
static void	free_dvector(double *v, int nl, int nh);

static double	**dmatrix(int nrl, int nrh, int ncl, int nch);
static void	free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

static double	***dcube(int nrl, int nrh, int ncl, int nch, int npl, int nph);
static void	free_dcube(double ***m, int nrl, int nrh, int ncl, int nch, int npl, int nph);

static void	tred2(double **a,int n,double *d,double *e);
static void	tqli(double *d,double *e,int n,double **z);
static void 	sort(int n, double *ra);
static void 	sort2(int n, double *ra, double *rb);
static void     sort5(int n, double *ra, double *rb, double *rc, double *rd, double *re);

static int      *alloc_ints(int len);
static double   *alloc_doubles(int len);
static int      **alloc_ptrs(int len);

static double   ran2(long *idum);
static double   rang(void);
static void     bad_input(void);

/**************** variables ***************************************************/

static  int     ndata;
static  int     ncodebook;
static  int     natoms;
static  int     nr_of_water_groups;
static  long int  nwater;
static  long int  maxnwater;
static  long int  nions;
static  long int  nbulkwater;
static  long int  ngrid,ngrid1;
static  long int  ti;
static  long int  timax;
static  long int  seed;
static  long int  max_atom_nr_in_pdb;
static  double  total_charge;
static  double  eps,epsstart,epsend,min_rhoaver2,rhoaver2,rhoav2,temp,distmin,rho0;
static  double  vol2,atot,prod,prod1,prod2,rhoaverstart,rhoaverend;
static  double  min_entropy,entropy,enttfac;
static 	double  a_vdw,b_vdw;
static  double  x_sphere; 
static  double  y_sphere;
static  double  z_sphere;
static  double  r_sphere;
static  double  boundary_radius;
static  double  water_shell_thickness;
static  double  *xdata,*ydata,*zdata;
static  double  *gauss_xm;
static  double  *lambda;
static  double  *gauss_dist;
static  double  **ev;
static	double	**cov;
static	double	**gauss_y;                /* codebook (r) vectors gauss_y[][r] */
static	double	***gauss_w;               /* codebook (r) directions w[][j][r] (j-th ev) */
static	double	**gauss_rho2;             /* codebook (r) variances rho2[][r] */
static	double	*gauss_vol;
static	double	*gauss_height;
static	double	*p;                 /* codebook (r) weight p[r] */
static	double	*g;                 /* unnormalized activities g[r] */
static	double	*a;                 /* normalized activities a[r] */
static	double	*xx,*xt,*xtt;
static	double  *xwater;
static	double  *ywater;
static	double  *zwater;
static	double  *xatom;
static	double  *yatom;
static	double  *zatom;
static	double  *qatom;
static	double  *sigatom;
static	double  *epsatom;
static	double  *a_vdw_atom;
static	double  *b_vdw_atom;
static	double  *xgrid;
static	double  *ygrid;
static	double  *zgrid;
static	double  *dgrid;
static	double  *densgrid;
static	int     *cutflags;
static  double  *distgrid;
static  double  *sdistgrid;
static 	double  *gradgrid;
static 	double  *curvegrid;
static 	double  *curvaturegrid;
static 	int     *replace_water_by_ion;
static	int     *groupflags;
static  int     write_surface_statistics_flag;
static  int     write_volume_statistics_flag;
static  int     use_gaussians_flag;
static  int     use_boundary_flag;
static  int     output_bulk_water_only_flag;
static  int     no_solute_present_flag;
static  int     output_no_solute_flag;
static  int     add_ions_flag;
static  int     output_mkpsf_script;
static  int     water_shell_thickness_flag;
static  int     ncodebook_flag;
static  int     enforced_total_charge_flag;
static  int     enforced_total_charge;


static void error(int errornr)
{
  printf("error nr. %d\n",errornr);
  exit(1);
}


static void terminate()
{
  fprintf(stderr,"\n\nSolvate terminated with fatal error.\n\n");
  exit(1);
}



static void here(int herenr)
{
  printf("here %d\n",herenr);
  return;
}


static void skip_line(FILE *fp)
{
  fgetc(fp);
  for ( ; fgetc(fp)!='\n' ; ) ;
  return;
}


static int     read_string(FILE *f,char *str)
{
        int     i,ch;
        i=0;
        str[0]=0;
        do {
                ch=fgetc(f);
                if ((ch==EOF)||(ch=='\n')) {
                        str[i]=0;
                        if (ch==EOF)    return -1;
                        else            return 0;
                } else {
                        str[i++]= (char)ch;
                }
                if (i>=MAX_STRING_SIZE-2) {
                        fprintf(stderr,"\n\nERROR: Line too long for string-buffer (read_string).\n");
                        fprintf(stderr,"Actual stringbuffer-size is %d.\n",MAX_STRING_SIZE);
                        str[i]=0;
                        return -1;
                }
        } while (ch!=EOF);
        str[i]=0;
        return -1;
}




static int     strup(char *str)
{
        int     i=0;
        int     toupper(int i);

        while (str[i]!=0) {
                str[i]=(char)toupper((int)str[i]);
                i++;
        }
        return 0;
}



static int     erase_spaces(char *str)
{
        int i,slen=(int)strlen(str);
        int j,k;

        for (i=0;i<slen;i++) {
                if (str[i]!=32) break;
        }
        if (str[i]==0) {
                str[0]=0;
                return 0;
        }
        if (i==0)       return 0;
        for (j=i,k=0;j<slen;j++,k++) {
                str[k]=str[j];
        }
        for(j=0;j<i;j++) {
                str[slen-1-j]=32;
        }
        return 0;
}



/* read_next_pdbdata is from the EGO_VIII-distribution, (c) M. Eichinger */

static int read_next_pdbdata(FILE *fh,int *atom_nr,double *ko,double *beta,double *harmon,char *atomstr,char *residstr,int *residnr,char *segidstr)
{
	int	res;
	double	resd;
	char    linestr[MAX_STRING_SIZE+1];
        char    tmpstr[MAX_STRING_SIZE+1];

                                        /* search for next line starting with 'ATOM' */
	while (read_string(fh,linestr)!=-1) {
		strup(linestr);
		strncpy(tmpstr,linestr,4);
		tmpstr[4]=0;
                if (strcmp(tmpstr,"ATOM")==0)   goto perform_atom_data; /* found one */
		tmpstr[3]=0;
                if (strcmp(tmpstr,"END")==0)    return -1;
	}
	return -1;

perform_atom_data:
				/********** get the atom-number ******/
	strncpy(tmpstr,linestr+PDB_P_ATOMNR,PDB_L_ATOMNR);
	tmpstr[PDB_L_ATOMNR]=0;
        if (sscanf(tmpstr,"%d",&res)!=1)        goto illegal_pdb_line;
	if (res<1)				goto illegal_pdb_line;
	*atom_nr=res;

				/********** get the atom-string *************/
	strncpy(tmpstr,linestr+PDB_P_ATOMSTR,PDB_L_ATOMSTR);
	tmpstr[PDB_L_ATOMSTR]=0;
	erase_spaces(tmpstr);
	strcpy(atomstr,tmpstr);

				/********** get the resid-string *************/
	strncpy(tmpstr,linestr+PDB_P_RESIDSTR,PDB_L_RESIDSTR);
	tmpstr[PDB_L_RESIDSTR]=0;
	erase_spaces(tmpstr);
	strcpy(residstr,tmpstr);

				/********** get the resid-number *************/
	strncpy(tmpstr,linestr+PDB_P_RESIDNR,PDB_L_RESIDNR);
	tmpstr[PDB_L_RESIDNR]=0;
        if (sscanf(tmpstr,"%d",&res)!=1)        goto illegal_pdb_line;
	if (res<1)				goto illegal_pdb_line;
	*residnr=res;

				/********* get the x-coordinate ************/
	strncpy(tmpstr,linestr+PDB_P_XKO,PDB_L_XKO);
	tmpstr[PDB_L_XKO]=0;
        if (sscanf(tmpstr,"%lf",&resd)!=1)      goto illegal_pdb_line;
	if (resd<-990.0)			goto unknown_coords;
	ko[0]=resd;
				/********* get the y-coordinate ************/
	strncpy(tmpstr,linestr+PDB_P_YKO,PDB_L_YKO);
	tmpstr[PDB_L_YKO]=0;
        if (sscanf(tmpstr,"%lf",&resd)!=1)      goto illegal_pdb_line;
	if (resd<-990.0)			goto unknown_coords;
	ko[1]=resd;
				/********* get the z-coordinate ************/
	strncpy(tmpstr,linestr+PDB_P_ZKO,PDB_L_ZKO);
	tmpstr[PDB_L_ZKO]=0;
        if (sscanf(tmpstr,"%lf",&resd)!=1)      goto illegal_pdb_line;
	if (resd<-990.0)			goto unknown_coords;
	ko[2]=resd;
				/********* get the beta-value ***************/
	strncpy(tmpstr,linestr+PDB_P_BETA,PDB_L_BETA);
	tmpstr[PDB_L_BETA]=0;
        if (sscanf(tmpstr,"%lf",&resd)!=1)      goto illegal_pdb_line;
	*beta=resd;
				/******* get the harmonic constant ***********/
	strncpy(tmpstr,linestr+PDB_P_HARMON,PDB_L_HARMON);
	tmpstr[PDB_L_HARMON]=0;
        if (sscanf(tmpstr,"%lf",&resd)!=1)      goto illegal_pdb_line;
	*harmon=resd;

				/********** get the segid-string *************/
	if ((int)strlen(linestr)<PDB_P_SEGIDSTR) {
		segidstr[0]=0;
	} else {
		strncpy(tmpstr,linestr+PDB_P_SEGIDSTR,PDB_L_SEGIDSTR);
		tmpstr[PDB_L_SEGIDSTR]=0;
		erase_spaces(tmpstr);
		strcpy(segidstr,tmpstr);
	}

	return 0;

illegal_pdb_line:
        fprintf(stderr,"\n\nFatal error: Illegal PDB-Line. Check row-positions\n");
        fprintf(stderr,"ATOM!atomn! !st! !re!!resnr!  !xko   !!yko   !!zko   !!beta!!harmon!    !se!\n");
	fprintf(stderr,linestr);
	return -2;

unknown_coords:
        fprintf(stderr,"\n\nFatal error: Unkown coordinates in PDB-File.\n");
	fprintf(stderr,linestr);
	return -2;
}




/* patch_pdb_coords is from the EGO_VIII-distribution, (c) M. Eichinger */


static int patch_pdb_coords(FILE *fh,FILE *fh1,double x,double y,double z)
{
	int	res;
	long	curpos;
	char    linestr[MAX_STRING_SIZE+1];
        char    tmpstr[MAX_STRING_SIZE+1];
                           /* search for next line starting with 'ATOM' */
	do {
		curpos=ftell(fh);
		res=read_string(fh,linestr);
		strup(linestr);
		strncpy(tmpstr,linestr,4);
		tmpstr[4]=0;
                if (strcmp(tmpstr,"ATOM")==0)   goto set_atom_data; /* found one */
		tmpstr[3]=0;
		fprintf(fh1,"%s\n",linestr);	
		fflush(fh1);
                if (strcmp(tmpstr,"END")==0)    return -1;
	} while(res!=-1);
	return -1;

set_atom_data:

        sprintf(tmpstr,"%8.3f        ",x);
	memcpy(linestr+PDB_P_XKO,tmpstr,PDB_L_XKO);

        sprintf(tmpstr,"%8.3f        ",y);
	memcpy(linestr+PDB_P_YKO,tmpstr,PDB_L_YKO);

        sprintf(tmpstr,"%8.3f        ",z);
	memcpy(linestr+PDB_P_ZKO,tmpstr,PDB_L_ZKO);

        fprintf(fh1,"%s\n",linestr);
	fflush(fh1);
	return 0;
}




static void add_pdb_tip(FILE *fp,int flag,double x,double y,double z,long int atom_nr,int res_id,int seg_id)
{
        char    tmpstr[MAX_STRING_SIZE+1];
        char    linestr[MAX_STRING_SIZE+1];

        if (flag==0) {
          strcpy(linestr,"ATOM 000000  OH2 TIP3    0       0.000   0.000   0.000  0.00  0.00      W000");
          }
        if (flag==1) {
          strcpy(linestr,"ATOM 000000  H1  TIP3    0       0.000   0.000   0.000  0.00  0.00      W000");
          }
        if (flag==2) {
          strcpy(linestr,"ATOM 000000  H2  TIP3    0       0.000   0.000   0.000  0.00  0.00      W000");
          }
        sprintf(tmpstr,"%6ld          ",atom_nr);
	memcpy(linestr+PDB_P_ATOMNR,tmpstr,PDB_L_ATOMNR);
        sprintf(tmpstr,"%8.3f        ",x);
	memcpy(linestr+PDB_P_XKO,tmpstr,PDB_L_XKO);
        sprintf(tmpstr,"%8.3f        ",y);
	memcpy(linestr+PDB_P_YKO,tmpstr,PDB_L_YKO);
        sprintf(tmpstr,"%8.3f        ",z);
	memcpy(linestr+PDB_P_ZKO,tmpstr,PDB_L_ZKO);
        sprintf(tmpstr,"%5d          ",res_id);
	memcpy(linestr+PDB_P_RESIDNR-1,tmpstr,PDB_L_RESIDNR);
        sprintf(tmpstr,"%3d          ",seg_id);
	memcpy(linestr+PDB_P_SEGIDSTR+1,tmpstr,PDB_L_SEGIDSTR-1);
        fprintf(fp,"%s\n",linestr);
	fflush(fp);
	return;
}



static void add_pdb_ion(FILE *fp,int iontype,double x,double y,double z,long int atom_nr,int res_id)
{
        char    tmpstr[MAX_STRING_SIZE+1];
        char    linestr[MAX_STRING_SIZE+1];

        if (iontype==1) {
          strcpy(linestr,"ATOM 000000  NA  INA     0       0.000   0.000   0.000  0.00  0.00      NA  ");
          }
        else if (iontype==-1) {
          strcpy(linestr,"ATOM 000000  CL  ICL     0       0.000   0.000   0.000  0.00  0.00      CL  ");
          }
	else {
	  fprintf(stderr,"\n\nERROR: Only sodium and chloride ions are supported here.\n");
	  terminate();
	  }
        sprintf(tmpstr,"%6ld          ",atom_nr);
	memcpy(linestr+PDB_P_ATOMNR,tmpstr,PDB_L_ATOMNR);
        sprintf(tmpstr,"%8.3f        ",x);
	memcpy(linestr+PDB_P_XKO,tmpstr,PDB_L_XKO);
        sprintf(tmpstr,"%8.3f        ",y);
	memcpy(linestr+PDB_P_YKO,tmpstr,PDB_L_YKO);
        sprintf(tmpstr,"%8.3f        ",z);
	memcpy(linestr+PDB_P_ZKO,tmpstr,PDB_L_ZKO);
        sprintf(tmpstr,"%5d          ",res_id);
	memcpy(linestr+PDB_P_RESIDNR-1,tmpstr,PDB_L_RESIDNR);
        fprintf(fp,"%s\n",linestr);
	fflush(fp);
	return;
}




static int goto_psf_atomdata(FILE *h_psf)
{
        static  char    lstr[MAX_STRING_SIZE];
        static  char    tmpstr[MAX_STRING_SIZE];
        int     nr_of_atoms;

        while (read_string(h_psf,lstr)!=-1) {
                strup(lstr);
                if (sscanf(lstr,"%d %s",&nr_of_atoms,tmpstr)==2) {
                        if (strcmp(tmpstr,"!NATOM")==0) {
                                if (nr_of_atoms!=natoms) {
	                          fprintf(stderr,"\n\nERROR: Atom counts in pdb- and psf-file differ.\n");
	                          terminate();
				  }
                                return 0;
                        }
                }
        }
	fprintf(stderr,"\n\nERROR: No atom-data-section (!NATOM) in PSF-file.\n");
	terminate();
        return -1;
}




static int read_next_psfatom(FILE *fh,int *atom_nr, char *atom_typestr, double *charge, double *mass, char *atomstr, char *residstr, int *residnr, char *segidstr)
{
	int	res;
	double	resd;
        static  char    lstr[MAX_STRING_SIZE];
        static  char    tmpstr[MAX_STRING_SIZE];

	if (read_string(fh,lstr)==-1)	return -1;
	strup(lstr);

	strncpy(tmpstr,lstr+PSF_P_ATOMNR,PSF_L_ATOMNR);
	tmpstr[PSF_L_ATOMNR]=0;
        if (sscanf(tmpstr,"%d",&res)!=1)        goto illegal_psf_atomline;
	if (res<0)				goto illegal_psf_atomline;
	*atom_nr=res;

				/********** get the segid-string *************/
	strncpy(tmpstr,lstr+PSF_P_SEGIDSTR,PSF_L_SEGIDSTR);
	tmpstr[PSF_L_SEGIDSTR]=0;
	erase_spaces(tmpstr);
	strcpy(segidstr,tmpstr);

				/********** get the resid-number *************/
	strncpy(tmpstr,lstr+PSF_P_RESIDNR,PSF_L_RESIDNR);
	tmpstr[PSF_L_RESIDNR]=0;
        if (sscanf(tmpstr,"%d",&res)!=1)        goto illegal_psf_atomline;
	if (res<1)				goto illegal_psf_atomline;
	*residnr=res;

				/********** get the resid-string *************/
	strncpy(tmpstr,lstr+PSF_P_RESIDSTR,PSF_L_RESIDSTR);
	tmpstr[PSF_L_RESIDSTR]=0;
	erase_spaces(tmpstr);
	strcpy(residstr,tmpstr);

				/********** get the atom-string *************/
	strncpy(tmpstr,lstr+PSF_P_ATOMSTR,PSF_L_ATOMSTR);
	tmpstr[PSF_L_ATOMSTR]=0;
	erase_spaces(tmpstr);
	strcpy(atomstr,tmpstr);

				/********** get the atomtype-string ********/
	strncpy(tmpstr,lstr+PSF_P_ATOMTYPSTR,PSF_L_ATOMTYPSTR);
	tmpstr[PSF_L_ATOMTYPSTR]=0;
	erase_spaces(tmpstr);
	strcpy(atom_typestr,tmpstr);

				/********** get the charge *************/
	strncpy(tmpstr,lstr+PSF_P_CHARGE,PSF_L_CHARGE);
	tmpstr[PSF_L_CHARGE]=0;
        if (sscanf(tmpstr,"%lf",&resd)!=1)      goto illegal_psf_atomline;
	*charge=resd;

				/********** get the mass *************/
	strncpy(tmpstr,lstr+PSF_P_MASS,PSF_L_MASS);
	tmpstr[PSF_L_MASS]=0;
        if (sscanf(tmpstr,"%lf",&resd)!=1)      goto illegal_psf_atomline;
	if (resd<0.5)				goto illegal_psf_atomline;
	*mass=resd;

	return 0;

illegal_psf_atomline:
        fprintf(stderr,"\n\nFatal error: Illegal PSF-atom-data-line. Check row-positions !\n");
        fprintf(stderr,"!atomnr! !se!!resd!!re! !as! !ts! !charge       ! !mass   !\n");
	fprintf(stderr,lstr);
	return -2;
}





static void read_pdb_file(char *infile)
{
	int i,res,resid_nr,atom_nr;
        double  beta,harmon;
        double  ko[3];
        char   atomstr[MAX_STRING_SIZE+1];
        char   residstr[MAX_STRING_SIZE+1];
        char   segidstr[MAX_STRING_SIZE+1];
        char   infile1[MAX_STRING_SIZE+1];
	FILE	*fp;


        strcpy(infile1,infile);
        strcat(infile1,".pdb");
	max_atom_nr_in_pdb=0;
	if (no_solute_present_flag==0) {
	  fprintf(stderr,"Reading solute from file %s ... ",infile1);
	  if ((fp=fopen(infile1,"r"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not read from file %s\n",infile1);
	    terminate();
	    }
          natoms=0;
          do {
            res=read_next_pdbdata(fp,&atom_nr,ko,&beta,&harmon,atomstr,residstr,&resid_nr,segidstr);
            if (res==-2) {
	      fprintf(stderr,"\n\nERROR: Error in pdb-file.\n");
              terminate();
              }
	    if (max_atom_nr_in_pdb<atom_nr) max_atom_nr_in_pdb=atom_nr;
            natoms++;
          } while (res!=-1);
          fclose(fp);
          natoms--;
	  }
	else {
	  natoms=0;
	  }
	if (natoms>=1) {
	  xatom=dvector(1,natoms);
	  yatom=dvector(1,natoms);
	  zatom=dvector(1,natoms);
	  qatom=dvector(1,natoms);
	  sigatom=dvector(1,natoms);
	  epsatom=dvector(1,natoms);
	  if ((fp=fopen(infile1,"r"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not read from file %s\n",infile1);
	    terminate();
	    }
	  for (i=1;i<=natoms;i++) {
            res=read_next_pdbdata(fp,&atom_nr,ko,&beta,&harmon,atomstr,residstr,&resid_nr,segidstr);
            if (res==-2) {
              terminate();
              }
	    xatom[i]=ko[0];
	    yatom[i]=ko[1];
	    zatom[i]=ko[2];
	    /* non-bonded vdw parameters are taken from XPLOR, param19.pro */
	    if (strncmp(atomstr,"H",1)==0) {
	      epsatom[i]=0.0498;
	      sigatom[i]=1.4254;
	      }
	    else if (strncmp(atomstr,"CH1",3)==0) {
	      epsatom[i]=0.0486;
	      sigatom[i]=4.2140;
	      }
	    else if (strncmp(atomstr,"CH2",3)==0) {
	      epsatom[i]=0.1142;
	      sigatom[i]=3.9823;
	      }
	    else if (strncmp(atomstr,"CH3",3)==0) {
	      epsatom[i]=0.1811;
	      sigatom[i]=3.8576;
	      }
	    else if (strncmp(atomstr,"C",1)==0) {
	      epsatom[i]=0.1200;
	      sigatom[i]=3.7418;
	      }
	    else if (strncmp(atomstr,"OC",2)==0) {
	      epsatom[i]=0.6469;
	      sigatom[i]=2.8509;
	      }
	    else if (strncmp(atomstr,"O",1)==0) {
	      epsatom[i]=0.1591;
	      sigatom[i]=2.8509;
	      }
	    else if (strncmp(atomstr,"N",1)==0) {
	      epsatom[i]=0.2384;
	      sigatom[i]=2.8509;
	      }
	    else if (strncmp(atomstr,"S",1)==0) {
	      epsatom[i]=0.0430;
	      sigatom[i]=3.3676;
	      }
	    else {
	      epsatom[i]=SIG_H2O;
	      sigatom[i]=0.15;
	      }
	    }
          fclose(fp);
	  if (natoms==1) fprintf(stderr,"%d atom read.\n",natoms);
	  if (natoms!=1) fprintf(stderr,"%d atoms read.\n",natoms);
	  }
	else {
	  if (no_solute_present_flag==0) {
	    fprintf(stderr,"No solute found; creating a pure water droplet at (0,0,0):\n\n");
	    }
	  else {
	    fprintf(stderr,"No solute given; creating a pure water droplet at (0,0,0):\n\n");
	    }
	  natoms=1;
	  xatom=dvector(1,natoms);
	  yatom=dvector(1,natoms);
          zatom=dvector(1,natoms);
	  sigatom=dvector(1,natoms);
	  epsatom=dvector(1,natoms);
	  xatom[1] = 0.0;
	  yatom[1] = 0.0;
	  zatom[1] = 0.0;
	  sigatom[1] = 0.0;
	  epsatom[1] = 0.15;
	  }
        return;
}


static void read_atom_charges_from_psf_file(char *infile)
{
	int    i,res,atom_index,residnr;
        double atom_charge,atom_mass;
        char   atomtypestr[MAX_STRING_SIZE+1];
        char   atomstr[MAX_STRING_SIZE+1];
        char   residstr[MAX_STRING_SIZE+1];
        char   segidstr[MAX_STRING_SIZE+1];
        char   infile1[MAX_STRING_SIZE+1];
	FILE *fp;


        strcpy(infile1,infile);
        strcat(infile1,".psf");
	if (no_solute_present_flag==0) {
	  total_charge = 0.0;
	  fprintf(stderr,"Reading atomic charges from file %s ... ",infile1);
	  if ((fp=fopen(infile1,"r"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not read from file %s\n",infile1);
	    terminate();
	    }
          goto_psf_atomdata(fp);
          for (i=1; i<=natoms; i++) {
            res=read_next_psfatom(fp,&atom_index,atomtypestr,&atom_charge,&atom_mass,atomstr,residstr,&residnr,segidstr);
            if (res==-2) {
	      fprintf(stderr,"\n\nERROR: Error in psf-file.\n");
              terminate();
              }
            qatom[i]=atom_charge;
	    total_charge += atom_charge;
            }
          fclose(fp);
	  fprintf(stderr,"total charge is %+5.2f e.\n",total_charge);
	  }
        return;
}


static double calc_maxdist(double xr, double yr, double zr)
{
	int ii;
	double dx,dy,dz,maxd,d2;

	maxd=0.0;
	for (ii=1; ii<=natoms; ii++) {
          dx=xr-xatom[ii];
          dy=yr-yatom[ii];
          dz=zr-zatom[ii];
	  d2=dx*dx+dy*dy+dz*dz;
	  if (d2>maxd) {
	    maxd=d2;
	    }
	  }
	if (maxd==0.0) maxd=1.0e-10;
        return sqrt(maxd);
}





static void calc_bounding_sphere()
{
	int i,flag;
	double stepwidth;
	double xmin,xmax;
	double ymin,ymax;
	double zmin,zmax;
	double xr,xrnew;
	double yr,yrnew;
	double zr,zrnew;
	double mr,mrnew;

	fprintf(stderr,"\nCalculating bounding sphere ... ");
        xmin=1.0e20;
        ymin=1.0e20;
        zmin=1.0e20;
        xmax=-1.0e20;
        ymax=-1.0e20;
        zmax=-1.0e20;
        for (i=1;i<=natoms;i++) {
          if (xatom[i]<xmin) xmin=xatom[i];
          if (xatom[i]>xmax) xmax=xatom[i];
          if (yatom[i]<ymin) ymin=yatom[i];
          if (yatom[i]>ymax) ymax=yatom[i];
          if (zatom[i]<zmin) zmin=zatom[i];
          if (zatom[i]>zmax) zmax=zatom[i];
          }
        xr=(xmin+xmax)/2.0;
        yr=(ymin+ymax)/2.0;
        zr=(zmin+zmax)/2.0;
	mr=calc_maxdist(xr,yr,zr);
	stepwidth=INITIAL_STEP_WIDTH_FOR_CALC_BOUNDING_SPHERE;
	do {
	  do {
	    flag=0;
	    for (i=0; i<50; i++) {
	      xrnew=xr-stepwidth*rang();
	      yrnew=yr-stepwidth*rang();
	      zrnew=zr-stepwidth*rang();
	      mrnew=calc_maxdist(xrnew,yrnew,zrnew);
	      if (mrnew<mr) {
	        mr=mrnew;
	        xr=xrnew;
	        yr=yrnew;
	        zr=zrnew;
	        flag=1;
		}
	      }
	    } while (flag==1);
	  stepwidth /= STEP_WIDTH_DECREASE_FOR_CALC_BOUNDING_SPHERE;
	  } while (stepwidth>=MINIMAL_STEP_WIDTH_FOR_CALC_BOUNDING_SPHERE);
	x_sphere=xr;
	y_sphere=yr;
	z_sphere=zr;
	r_sphere=mr;
	if (r_sphere >= boundary_radius) {
	  fprintf(stderr,"\n\nERROR: Selected maximum boundary radius (%f) is smaller\n",boundary_radius);
	  fprintf(stderr,"       than calculated bounding sphere radius (%f).\n",r_sphere);
	  fprintf(stderr,"       Use a larger boundary radius!\n");
	  terminate();
	  }
	fprintf(stderr,"done.\n");
        return;
}


static void cut_away(double theta, double phi)
{
	long int i,j,k;
	double xx,yy,zz,dis;
	double d11,d12,d13;
	double d21,d22,d23;
	double d31,d32,d33;
	double dis_closest;
	double **d1, **d2, **dd;

        d1=dmatrix(1,3,1,3);
        d2=dmatrix(1,3,1,3);
        dd=dmatrix(1,3,1,3);
	d1[1][1]=cos(theta/180.0*PI);
	d1[1][2]=0.0;
	d1[1][3]=sin(theta/180.0*PI);
	d1[2][1]=0.0;
	d1[2][2]=1.0;
	d1[2][3]=0.0;
	d1[3][1]=-sin(theta/180.0*PI);
	d1[3][2]=0.0;
	d1[3][3]=cos(theta/180.0*PI);
	d2[1][1]=cos(phi/180.0*PI);
	d2[1][2]=sin(phi/180.0*PI);
	d2[1][3]=0.0;
	d2[2][1]=-sin(phi/180.0*PI);
	d2[2][2]=cos(phi/180.0*PI);
	d2[2][3]=0.0;
	d2[3][1]=0.0;
	d2[3][2]=0.0;
	d2[3][3]=1.0;
	for (i=1;i<=3;i++) {
	  for (j=1;j<=3;j++) {
	    dd[i][j]=0.0;
	    for (k=1;k<=3;k++) {
	      dd[i][j] += (d1[i][k]*d2[k][j]);
	      }
            } 
          }
	d11=dd[1][1];
	d12=dd[1][2];
	d13=dd[1][3];
	d21=dd[2][1];
	d22=dd[2][2];
	d23=dd[2][3];
	d31=dd[3][1];
	d32=dd[3][2];
	d33=dd[3][3];
        free_dmatrix(d1,1,3,1,3);
        free_dmatrix(d2,1,3,1,3);
        free_dmatrix(dd,1,3,1,3);
	dis_closest=-1.0e20;
        for (i=1; i<=natoms; i++) {
	  xx = d11*xatom[i] + d12*yatom[i] + d13*zatom[i];
	  yy = d21*xatom[i] + d22*yatom[i] + d23*zatom[i];
	  zz = d31*xatom[i] + d32*yatom[i] + d33*zatom[i];
	  dis=zz+(xx*xx+yy*yy)/(2.0*boundary_radius); 
          if (dis >= dis_closest) dis_closest=dis;
	  }
	dis_closest += 0.5;
        for (i=1; i<=ngrid; i++) {
	  if (cutflags[i]==1) {
	    xx = d11*xgrid[i] + d12*ygrid[i] + d13*zgrid[i];
	    yy = d21*xgrid[i] + d22*ygrid[i] + d23*zgrid[i];
	    zz = d31*xgrid[i] + d32*ygrid[i] + d33*zgrid[i];
	    dis=zz+(xx*xx+yy*yy)/(2.0*boundary_radius); 
	    if (dis>dis_closest) cutflags[i]=0;
	    }
	  }
	return;
}



static void create_convex_volume(void)
{
	long int i,ii,maxn;
	double x,y,z;
	double gridmax;
	double theta,phi,r2,phistep;
	double *xgrid1;
	double *ygrid1;
	double *zgrid1;

	fprintf(stderr,"Creating convex volume (max. radius is %f) ... ",boundary_radius);
        for (i=1; i<=natoms; i++) {
          xatom[i] -= x_sphere;
          yatom[i] -= y_sphere;
          zatom[i] -= z_sphere;
	  }
	gridmax=floor(r_sphere)+2.0;
	maxn=floor(MAX_NR_OF_GRID_POINTS_FACTOR*4.0/3.0*PI*gridmax*gridmax*gridmax);
	xgrid = dvector(1,maxn);
	ygrid = dvector(1,maxn);
	zgrid = dvector(1,maxn);
	cutflags = ivector(1,maxn);
	i=1;
	r2=gridmax*gridmax;
	for (x=-gridmax; x<gridmax+0.1; x++) {
	  for (y=-gridmax; y<gridmax+0.1; y++) {
	    for (z=-gridmax; z<gridmax+0.1; z++) {
	      if ((x*x+y*y+z*z)<=r2) {
		xgrid[i]=x;
		ygrid[i]=y;
		zgrid[i]=z;
		cutflags[i]=1;
		i++;
		if (i>=maxn) {
		  fprintf(stderr,"\n\nERROR: Maximum number of grid points exceeded.\n");
		  fprintf(stderr,"       Increase MAX_NR_OF_GRID_POINTS_FACTOR!\n");
		  terminate();
		  }
		}
	      }
	    }
          }
	ngrid=i-1;

	cut_away(0.0,0.0);
	cut_away(90.0,0.0);
	cut_away(90.0,90.0);
	cut_away(90.0,180.0);
	cut_away(90.0,270.0);
	cut_away(180.0,0.0);
	ngrid1=0;
	for (i=1; i<=ngrid; i++) {
	  if (cutflags[i]==1) {
	    ngrid1++;
	    }
	  }
	xgrid1=dvector(1,ngrid1);
	ygrid1=dvector(1,ngrid1);
	zgrid1=dvector(1,ngrid1);
	ii=1;
	for (i=1; i<=ngrid; i++) {
	  if (cutflags[i]==1) {
	    xgrid1[ii] = xgrid[i];
	    ygrid1[ii] = ygrid[i];
	    zgrid1[ii] = zgrid[i];
	    ii++;
	    }
	  }
	free_dvector(xgrid,1,maxn);
	free_dvector(ygrid,1,maxn);
	free_dvector(zgrid,1,maxn);
	free_ivector(cutflags,1,maxn);
	xgrid=xgrid1;
	ygrid=ygrid1;
	zgrid=zgrid1;
	ngrid=ii-1;
	cutflags = ivector(1,ngrid);
	for (i=1; i<=ngrid; i++) cutflags[i]=1;
	for (theta=ANGLE_STEP_FOR_CUT_AWAY; theta<180.1-ANGLE_STEP_FOR_CUT_AWAY; theta+=ANGLE_STEP_FOR_CUT_AWAY) {
	  phistep=360.0/floor((360.0-10.0*ANGLE_STEP_FOR_CUT_AWAY)/ANGLE_STEP_FOR_CUT_AWAY*sin(theta*PI/180.0)+10.0);
	  for (phi=0.0; phi<360.1-phistep; phi+=phistep) {
	    cut_away(theta,phi);
	    }
	  }
	ngrid1=0;
	for (i=1; i<=ngrid; i++) {
	  if (cutflags[i]==1) {
	    ngrid1++;
	    }
	  }
	xgrid1=dvector(1,ngrid1);
	ygrid1=dvector(1,ngrid1);
	zgrid1=dvector(1,ngrid1);
	ii=1;
	for (i=1; i<=ngrid; i++) {
	  if (cutflags[i]==1) {
	    xgrid1[ii] = xgrid[i];
	    ygrid1[ii] = ygrid[i];
	    zgrid1[ii] = zgrid[i];
	    ii++;
	    }
	  }
	free_dvector(xgrid,1,ngrid);
	free_dvector(ygrid,1,ngrid);
	free_dvector(zgrid,1,ngrid);
	free_ivector(cutflags,1,ngrid);
	xgrid=xgrid1;
	ygrid=ygrid1;
	zgrid=zgrid1;
	ngrid=ii-1;
	for (i=1; i<=natoms; i++) {
	  xatom[i] += x_sphere;
	  yatom[i] += y_sphere;
	  zatom[i] += z_sphere;
	  }
	for (i=1; i<=ngrid; i++) {
	  xgrid[i] += x_sphere;
	  ygrid[i] += y_sphere;
	  zgrid[i] += z_sphere;
	  }
	fprintf(stderr,"done.\n");
	return;
}





static void dice_data()
{
	int i,gridnr;
	double sigma;

	fprintf(stderr,"Gauss-sampling convex volume ... ");
        ndata=NR_OF_DICE_DATA;
	sigma=water_shell_thickness;
	xdata = dvector(1,ndata);
	ydata = dvector(1,ndata);
	zdata = dvector(1,ndata);
	for (i=1;i<=ndata;i++) {
	  gridnr=1+floor(ngrid*ran2(&seed));
	  if ((gridnr<1)||(gridnr>ngrid)) {
	    fprintf(stderr,"\n\nERROR: This should not happen. Sorry.\n");
	    terminate();
	    }
	  xdata[i]=xgrid[gridnr]+sigma*rang();
	  ydata[i]=ygrid[gridnr]+sigma*rang();
	  zdata[i]=zgrid[gridnr]+sigma*rang();
	  }
	free_dvector(xgrid,1,ngrid);
	free_dvector(ygrid,1,ngrid);
	free_dvector(zgrid,1,ngrid);
	fprintf(stderr,"done.\n");
	return;
}


static void do_pca()
{
	int	i,j;
        double  *e;

	cov = dmatrix(1,3,1,3);
	gauss_xm = dvector(1,3);
	gauss_y = dmatrix(1,3,1,ncodebook);
	gauss_w = dcube(1,3,1,3,1,ncodebook);
	gauss_rho2 = dmatrix(1,3,1,ncodebook);
	p = dvector(1,ncodebook);
	ev = dmatrix(1,3,1,3);
	lambda = dvector(1,3);
	g=dvector(1,ncodebook);
	a=dvector(1,ncodebook);
	gauss_dist=dvector(1,ncodebook);
	xx=dvector(1,3);
	xt=dvector(1,3);
	xtt=dvector(1,3);

	/*** calculate covariance matrix: ***/

	for (i=1;i<=3;i++) {
	  gauss_xm[i]=0.0;
	  }
	for (i=1;i<=ndata;i++) {
	  gauss_xm[1] += xdata[i];
	  gauss_xm[2] += ydata[i];
	  gauss_xm[3] += zdata[i];
	  }
	for (i=1;i<=3;i++) {
	  gauss_xm[i] /= (1.0*ndata);
	  }
	for (i=1;i<=3;i++) {
	  for (j=1;j<=3;j++) {
	    cov[i][j]=0.0;
	    }
	  }
	for (i=1;i<=ndata;i++) {
	  cov[1][1] += (xdata[i]-gauss_xm[1])*(xdata[i]-gauss_xm[1]);
	  cov[1][2] += (xdata[i]-gauss_xm[1])*(ydata[i]-gauss_xm[2]);
	  cov[1][3] += (xdata[i]-gauss_xm[1])*(zdata[i]-gauss_xm[3]);
	  cov[2][1] += (ydata[i]-gauss_xm[2])*(xdata[i]-gauss_xm[1]);
	  cov[2][2] += (ydata[i]-gauss_xm[2])*(ydata[i]-gauss_xm[2]);
	  cov[2][3] += (ydata[i]-gauss_xm[2])*(zdata[i]-gauss_xm[3]);
	  cov[3][1] += (zdata[i]-gauss_xm[3])*(xdata[i]-gauss_xm[1]);
	  cov[3][2] += (zdata[i]-gauss_xm[3])*(ydata[i]-gauss_xm[2]);
	  cov[3][3] += (zdata[i]-gauss_xm[3])*(zdata[i]-gauss_xm[3]);
	  }
	for (i=1;i<=3;i++) {
	  for (j=1;j<=3;j++) {
	    cov[i][j] /= (1.0*ndata);
	    }
	  }

        /*** diagonalize covariance-matrix: ***/

	e = dvector(1,3);
	for (i=1;i<=3;i++) {
	  for (j=1;j<=3;j++) {
	    ev[i][j] = cov[i][j];
	    }
	  }
	tred2(ev,3,lambda,e);
	tqli(lambda,e,3,ev);  /* the eigenvectors j are now in ev[][j] */
	free_dvector(e,1,3);

	/* ugly sort */
	if (lambda[1]<lambda[2]) {
	  temp=lambda[1]; lambda[1]=lambda[2]; lambda[2]=temp;
	  temp=ev[1][1]; ev[1][1]=ev[1][2]; ev[1][2]=temp;
	  temp=ev[2][1]; ev[2][1]=ev[2][2]; ev[2][2]=temp;
	  temp=ev[3][1]; ev[3][1]=ev[3][2]; ev[3][2]=temp;
	  }
	if (lambda[1]<lambda[3]) {
	  temp=lambda[1]; lambda[1]=lambda[3]; lambda[3]=temp;
	  temp=ev[1][1]; ev[1][1]=ev[1][3]; ev[1][3]=temp;
	  temp=ev[2][1]; ev[2][1]=ev[2][3]; ev[2][3]=temp;
	  temp=ev[3][1]; ev[3][1]=ev[3][3]; ev[3][3]=temp;
	  }
	if (lambda[2]<lambda[3]) {
	  temp=lambda[2]; lambda[2]=lambda[3]; lambda[3]=temp;
	  temp=ev[1][2]; ev[1][2]=ev[1][3]; ev[1][3]=temp;
	  temp=ev[2][2]; ev[2][2]=ev[2][3]; ev[2][3]=temp;
	  temp=ev[3][2]; ev[3][2]=ev[3][3]; ev[3][3]=temp;
	  }
        return;
}


static void init_codebook_positions()
{
	int	i,r;

	for (r=1;r<=ncodebook;r++) {
	  for (i=1;i<=3;i++) {
	    gauss_y[i][r]=gauss_xm[i];
	    }
	  }
        return;
}



static void init_codebook_directions()
{
	int	i,j,r;

	for (r=1;r<=ncodebook;r++) {
	  for (i=1;i<=3;i++) {
	    for (j=1;j<=3;j++) {
	      gauss_w[i][j][r]=ev[i][j];
	      }
	    }
	  }
        return;
}



static void init_codebook_variances()
{
	int	i,r;

	rho0=sqrt(lambda[1]);
	for (r=1;r<=ncodebook;r++) {
	  for (i=1;i<=3;i++) {
	    gauss_rho2[i][r]=FACTOR_FOR_INITIAL_LAMBDA_FOR_LEARNS_GAUSSIANS*lambda[i];
	    }
	  }
        return;
}



static void init_parameters()
{
	int r;

	entropy=1.0e10;
	min_entropy=1.0e10;
	min_rhoaver2=1.0e10;
	enttfac=ENTROPY_AVERAGE_LENGTH/(1.0*timax);

        /*** init gauss weights p[r] ***/

	for (r=1;r<=ncodebook;r++) {
	  p[r]=1.0/(1.0*ncodebook);
	  }
	return;
}



static void compute_activities()
{
        int r,i,j;

	for (r=1;r<=ncodebook;r++) {
	  for (i=1;i<=3;i++) {
	    xt[i]=0.0;
	    }
	  for (i=1;i<=3;i++) {
	    for (j=1;j<=3;j++) {
	      xt[i] += (gauss_w[j][i][r]*(xx[j]-gauss_y[j][r]));   /* XT=W^T*(x-y_r) */
	      }
            }
	  for (i=1;i<=3;i++) {
            xt[i] /= gauss_rho2[i][r];                         /* XT=GT*XT */
            }
	  for (i=1;i<=3;i++) {
	    xtt[i]=0.0;
	    }
	  for (i=1;i<=3;i++) {
	    for (j=1;j<=3;j++) {
	      xtt[i] += (gauss_w[i][j][r]*xt[j]);           /* XTT=W*XT */
	      }
            }
	  gauss_dist[r]=0.0;
	  for (j=1;j<=3;j++) {
	    gauss_dist[r] += ((xx[j]-gauss_y[j][r])*xtt[j]);    /* dist_r=(x-y_r)^T*XTT */
	    }
	  gauss_dist[r] *= 0.5;
	  }
	distmin=1.0e200;
	for (r=1;r<=ncodebook;r++) {
	  if (gauss_dist[r] <= distmin) distmin=gauss_dist[r];
	  }
	for (r=1;r<=ncodebook;r++) {
	  gauss_dist[r] -= distmin;
	  }
	for (r=1;r<=ncodebook;r++) {
	  if (gauss_dist[r]>100.0) {
	    g[r]=0.0;
	    }
          else {
	    g[r]=exp(-gauss_dist[r]);
	    }
	  vol2=pow(2.0*PI,3.0/2.0);
	  for (i=1;i<=3;i++) {
	    vol2 *= sqrt(gauss_rho2[i][r]);
	    }
	  g[r] /= vol2;
	  }
	atot=0.0;
	for (r=1;r<=ncodebook;r++) {
	  atot += (p[r]*g[r]);
	  }
	for (r=1;r<=ncodebook;r++) {
	  a[r] = (p[r]*g[r])/atot;
	  }

	/* sum moving-averaged kulbak-entropy: */
	entropy = (1.0-enttfac)*entropy - enttfac*(log(atot)-distmin);
	if (entropy < min_entropy) {
	  min_entropy=entropy;
	  min_rhoaver2=rhoaver2;
	  }
        return;
}


static void learn_codebook_positions(double eps)
{
	  int i,r;

	  /* learn codebook vectors y[][r]  */
	  /* y_r := y_r + eps*[a_r*(x-y_r) + noise]  */

	  for (r=1;r<=ncodebook;r++) {
	    for (i=1;i<=3;i++) {
	      gauss_y[i][r] += (eps*(a[r]*(xx[i]-gauss_y[i][r])
			       +NOISE_FOR_LEARN_CODEBOOK_POSITIONS*rang()/sqrt(gauss_rho2[1][r])));
	      }
	    }
          return;
}




static void learn_codebook_directions(double eps)
{
	  int i,j,k,r;

	  /* learn codebook directions w[][j][r]  */
	  /* w_j = w_j + eps*a_r*(w_j*(x-y))*(x-y) */

	  for (r=1;r<=ncodebook;r++) {
	    for (j=1;j<=3;j++) {
	      prod=0.0;
	      for (i=1;i<=3;i++) {
		prod += (gauss_w[i][j][r]*(xx[i]-gauss_y[i][r]));
		}
	      for (i=1;i<=3;i++) {
		gauss_w[i][j][r] += ((eps*a[r]*prod*(xx[i]-gauss_y[i][r]))/(rho0*rho0));
		}
	      }
	    }

	  /* ortho-normalize codebook directions w[][j][r] (schmidt) */

	  for (r=1;r<=ncodebook;r++) {
	    for (k=2;k<=3;k++) {
	      for (j=1;j<k;j++) {
		prod1=0.0;
		prod2=0.0;
	        for (i=1;i<=3;i++) {
                  prod1 += (gauss_w[i][k][r]*gauss_w[i][j][r]);
                  prod2 += (gauss_w[i][j][r]*gauss_w[i][j][r]);
		  }
	        for (i=1;i<=3;i++) {
		  gauss_w[i][k][r] -= (prod1/prod2*gauss_w[i][j][r]);
		  }
		}
   	      }
	    for (j=1;j<=3;j++) {
	      prod=0.0;
	      for (i=1;i<=3;i++) {
		prod += (gauss_w[i][j][r]*gauss_w[i][j][r]);
		}
              prod=1.0/sqrt(prod);
	      for (i=1;i<=3;i++) {
		gauss_w[i][j][r] *= prod;
		}
              }   
            }
          return;
}



static void learn_codebook_variances_soft_individual(double eps, double mu)
{
	  int r,i,j;

	  /* soft volume (2): individual constraint to rhoaver */
	  for (r=1;r<=ncodebook;r++) {
	    for (j=1;j<=3;j++) {
	      prod=0.0;
	      for (i=1;i<=3;i++) {
		prod += (gauss_w[i][j][r]*(xx[i]-gauss_y[i][r]));
		}
	      prod2=prod*prod;
	      gauss_rho2[j][r] += (eps*a[r]*((prod2-gauss_rho2[j][r])-mu*(gauss_rho2[j][r]-rhoaver2)));
	      }
	    }
          return;
}



static void learn_codebook_variances_soft_global(double eps, double mu)
{
	  int r,i,j;
	  double mufactor;

	  /* soft volume (1): global */
	  for (r=1;r<=ncodebook;r++) {
	    rhoav2=1.0;
	    for (j=1;j<=3;j++) {
	      rhoav2 *= (gauss_rho2[j][r]);
	      }
	    rhoav2=pow(rhoav2,1.0/3.0);
	    mufactor=mu*(1.0-rhoaver2/rhoav2);
	    for (j=1;j<=3;j++) {
	      prod=0.0;
	      for (i=1;i<=3;i++) {
		prod += (gauss_w[i][j][r]*(xx[i]-gauss_y[i][r]));
		}
	      prod2=prod*prod;
	      gauss_rho2[j][r] += (eps*a[r]*((prod2-gauss_rho2[j][r])-mufactor*gauss_rho2[j][r]));
	      }
	    }
          return;
}


static void equalize_codebook_variances()
{
	  int r,j;

	  /* this routine is only used if natoms==1, so that we have *really* a sphere */
	  for (r=1;r<=ncodebook;r++) {
	    rhoav2=1.0;
	    for (j=1;j<=3;j++) {
	      rhoav2 *= (gauss_rho2[j][r]);
	      }
	    rhoav2=pow(rhoav2,1.0/3.0);
	    for (j=1;j<=3;j++) {
	      gauss_rho2[j][r] = rhoav2;
	      }
	    }
          return;
}



static void learn_codebook_variances_hard(double eps)
{
	  int r,i,j;
	  double mufactor;

	  /* hard volume */
	  for (r=1;r<=ncodebook;r++) {
	    for (j=1;j<=3;j++) {
	      prod=0.0;
	      for (i=1;i<=3;i++) {
		prod += (gauss_w[i][j][r]*(xx[i]-gauss_y[i][r]));
		}
	      prod2=prod*prod;
	      gauss_rho2[j][r] += eps*a[r]*(prod2-gauss_rho2[j][r]);
	      }
	    rhoav2=1.0;
	    for (j=1;j<=3;j++) {
	      rhoav2 *= (gauss_rho2[j][r]);
	      }
	    rhoav2=pow(rhoav2,1.0/3.0);
	    mufactor=rhoaver2/rhoav2;
	    for (j=1;j<=3;j++) {
	      gauss_rho2[j][r] *= mufactor;
	      }
	    }
          return;
}


static void outputs()
{
	  int r,i,j;

	  /*** output y[][r], sqrt(rho2[j][r]), w[i][j][r] ***/

	  printf("%f\n",sqrt(rhoaver2));
	  printf("%f\n",entropy);
	  for (r=1;r<=ncodebook;r++) {
	    for (j=1;j<=3;j++) {
	      printf("%f\n",gauss_y[j][r]);
	      }
	    for (j=1;j<=3;j++) {
	      printf("%f\n",sqrt(gauss_rho2[j][r]));
	      }
	    for (j=1;j<=3;j++) {
	      for (i=1;i<=3;i++) {
		printf("%f ",gauss_w[i][j][r]);
	        }
	      printf("\n");
	      }
	    }
	  return;

}



static double get_rhoaver2(double r1, double r2, long int tmax, long int t)
{
	  double rav;

	  rav = r1*pow((r2/r1),((1.0*t)/(1.0*tmax)));
	  return (rav*rav);
}


static double get_eps(double e1, double e2, long int tmax, long int t)
{
	  return (e1*pow((e2/e1),((1.0*t)/(1.0*tmax))));
}





static void learn_gaussians()
{
	if (ncodebook==1) {
	  fprintf(stderr,"Approximating volume with 1 gaussian: Phase 1/4...");
	  }
	else {
	  fprintf(stderr,"Approximating volume with %d gaussians: Phase 1/4...",ncodebook);
	  }
	do_pca();
	timax=TIMAX_FACTOR*ndata;
	init_parameters();
	init_codebook_positions();
	init_codebook_directions();
	init_codebook_variances();
        rhoaver2 = pow((gauss_rho2[1][1]*gauss_rho2[2][1]*gauss_rho2[3][1]),1.0/(1.0*3.0));
	rhoaverstart=sqrt(rhoaver2)*RHOAVERSTART_FACTOR;
	rhoaverend=sqrt(rhoaver2)*RHOAVEREND_FACTOR;
	for (ti=0;ti<timax;ti++) {
	  xx[1]=xdata[1+ti%ndata];
	  xx[2]=ydata[1+ti%ndata];
	  xx[3]=zdata[1+ti%ndata];
	  rhoaver2 = get_rhoaver2(rhoaverstart,rhoaverend,timax,ti);
	  compute_activities();
          learn_codebook_positions(EPS_CONST_FOR_LEARN_POSITIONS);
          learn_codebook_directions(EPS_CONST_FOR_LEARN_DIRECTIONS);
          learn_codebook_variances_hard(EPS_CONST_FOR_LEARN_VARIANCES);
	  }

	fprintf(stderr,"2/4...");
	rhoaverend = sqrt(min_rhoaver2);
	timax=TIMAX_FACTOR*ndata;
	init_parameters();
	init_codebook_positions();
	init_codebook_directions();
	init_codebook_variances();
        rhoaver2 = pow((gauss_rho2[1][1]*gauss_rho2[2][1]*gauss_rho2[3][1]),1.0/(1.0*3.0));
	rhoaverstart=sqrt(rhoaver2)*RHOAVERSTART_FACTOR;
	for (ti=0;ti<timax;ti++) {
	  xx[1]=xdata[1+ti%ndata];
	  xx[2]=ydata[1+ti%ndata];
	  xx[3]=zdata[1+ti%ndata];
	  rhoaver2 = get_rhoaver2(rhoaverstart,rhoaverend,timax,ti);
	  compute_activities();
          learn_codebook_positions(EPS_CONST_FOR_LEARN_POSITIONS);
          learn_codebook_directions(EPS_CONST_FOR_LEARN_DIRECTIONS);
          learn_codebook_variances_hard(EPS_CONST_FOR_LEARN_VARIANCES);
	  }

	fprintf(stderr,"3/4...");
	timax=TIMAX_FACTOR*ndata;
	init_parameters();
        rhoaver2 = rhoaverend*rhoaverend;
	epsstart=EPS_START_FOR_LEARN_DIRECTIONS;
	epsend  =EPS_END_FOR_LEARN_DIRECTIONS;
	for (ti=0;ti<timax;ti++) {
	  xx[1]=xdata[1+ti%ndata];
	  xx[2]=ydata[1+ti%ndata];
	  xx[3]=zdata[1+ti%ndata];
	  compute_activities();
	  eps = get_eps(epsstart,epsend,timax,ti);
          learn_codebook_positions(eps*EPS_CONST_FOR_LEARN_POSITIONS/EPS_CONST_FOR_LEARN_DIRECTIONS);
          learn_codebook_directions(eps);
          learn_codebook_variances_hard(eps);
	  }

	fprintf(stderr,"4/4...");
	timax=TIMAX_FACTOR*ndata;
	init_parameters();
        rhoaver2 = rhoaverend*rhoaverend;
	epsstart=EPS_START_FOR_LEARN_DIRECTIONS;
	epsend  =EPS_END_FOR_LEARN_DIRECTIONS;
	for (ti=0;ti<timax;ti++) {
	  xx[1]=xdata[1+ti%ndata];
	  xx[2]=ydata[1+ti%ndata];
	  xx[3]=zdata[1+ti%ndata];
	  compute_activities();
	  eps = get_eps(epsstart,epsend,timax,ti);
          learn_codebook_directions(eps);
          learn_codebook_variances_soft_individual(eps,MU_FOR_LEARN_VARIANCES_INDIVIDUAL);
	  }
	if (natoms<=1) {
	  equalize_codebook_variances();
	  }
	free_dvector(xdata,1,ndata);
	free_dvector(ydata,1,ndata);
	free_dvector(zdata,1,ndata);
	free_dvector(gauss_xm,1,3);
	free_dmatrix(cov,1,3,1,3);
	free_dvector(p,1,ncodebook);
	free_dmatrix(ev,1,3,1,3);
	free_dvector(lambda,1,3);
	free_dvector(g,1,ncodebook);
	free_dvector(a,1,ncodebook);
	free_dvector(gauss_dist,1,ncodebook);
	free_dvector(xx,1,3);
	free_dvector(xt,1,3);
	free_dvector(xtt,1,3);
	fprintf(stderr,"done.\n");
        return;
}





static double density_varrho(double x, double y, double z, double rhofac)
{
	int r;
        double sum,rhoq,di;
	double w11,w12,w13;
	double w21,w22,w23;
	double w31,w32,w33;
	double xt1,xt2,xt3,xtt1,xtt2,xtt3;
	double dx,dy,dz;

	sum=0.0;
	rhoq=2.0*rhofac*rhofac;
	for (r=1; r<=ncodebook; r++) {
          dx=x-gauss_y[1][r];
          dy=y-gauss_y[2][r];
          dz=z-gauss_y[3][r];
	  w11=gauss_w[1][1][r];
	  w12=gauss_w[1][2][r];
	  w13=gauss_w[1][3][r];
	  w21=gauss_w[2][1][r];
	  w22=gauss_w[2][2][r];
	  w23=gauss_w[2][3][r];
	  w31=gauss_w[3][1][r];
	  w32=gauss_w[3][2][r];
	  w33=gauss_w[3][3][r];
	  xt1 = (w11*dx + w21*dy + w31*dz)/gauss_rho2[1][r];
	  xt2 = (w12*dx + w22*dy + w32*dz)/gauss_rho2[2][r];
	  xt3 = (w13*dx + w23*dy + w33*dz)/gauss_rho2[3][r];
	  xtt1 = (w11*xt1) + (w12*xt2) + (w13*xt3);
	  xtt2 = (w21*xt1) + (w22*xt2) + (w23*xt3);
	  xtt3 = (w31*xt1) + (w32*xt2) + (w33*xt3);
	  di=(dx*xtt1 + dy*xtt2 + dz*xtt3)/rhoq;
	  if (di<100.0) sum += (exp(-di)/gauss_vol[r]);
	  }
	return (sum);
}




static double density(double x, double y, double z)
{
	int r;
        double sum,di;
	double w11,w12,w13;
	double w21,w22,w23;
	double w31,w32,w33;
	double xt1,xt2,xt3,xtt1,xtt2,xtt3;
	double dx,dy,dz;

	sum=0.0;
	for (r=1; r<=ncodebook; r++) {
          dx=x-gauss_y[1][r];
          dy=y-gauss_y[2][r];
          dz=z-gauss_y[3][r];
	  w11=gauss_w[1][1][r];
	  w12=gauss_w[1][2][r];
	  w13=gauss_w[1][3][r];
	  w21=gauss_w[2][1][r];
	  w22=gauss_w[2][2][r];
	  w23=gauss_w[2][3][r];
	  w31=gauss_w[3][1][r];
	  w32=gauss_w[3][2][r];
	  w33=gauss_w[3][3][r];
	  xt1 = (w11*dx + w21*dy + w31*dz)/gauss_rho2[1][r];
	  xt2 = (w12*dx + w22*dy + w32*dz)/gauss_rho2[2][r];
	  xt3 = (w13*dx + w23*dy + w33*dz)/gauss_rho2[3][r];
	  xtt1 = (w11*xt1) + (w12*xt2) + (w13*xt3);
	  xtt2 = (w21*xt1) + (w22*xt2) + (w23*xt3);
	  xtt3 = (w31*xt1) + (w32*xt2) + (w33*xt3);
	  di=0.5*(dx*xtt1 + dy*xtt2 + dz*xtt3);
	  if (di<100.0) sum += (gauss_height[r]*exp(-di));
	  }
	return (sum);
}


static double density_grad_curve(double x, double y, double z, double *gx, double *gy, double *gz, double *g, double *cr, double *curvature)
{
	int r;
        double sum,prod;
	double w11,w12,w13;
	double w21,w22,w23;
	double w31,w32,w33;
	double t11,t12,t13;
	double t21,t22,t23;
	double t31,t32,t33;
	double a11,a12,a13;
	double a21,a22,a23;
	double a31,a32,a33;
	double m11,m12,m13;
	double m21,m22,m23;
	double m31,m32,m33;
	double cxx,cxy,cxz,cyx,cyy,cyz,czx,czy,czz;
	double dx,dy,dz,ex;
	double g0x,g0y,g0z;
	double v1,v2,v3;
	double temp;

	sum=0.0;
	*gx = 0.0;
	*gy = 0.0;
	*gz = 0.0;
	cxx = 0.0;
	cxy = 0.0;
	cxz = 0.0;
	cyx = 0.0;
	cyy = 0.0;
	cyz = 0.0;
	czx = 0.0;
	czy = 0.0;
	czz = 0.0;
	for (r=1; r<=ncodebook; r++) {
          dx=x-gauss_y[1][r];
          dy=y-gauss_y[2][r];
          dz=z-gauss_y[3][r];
	  w11=gauss_w[1][1][r];
	  w12=gauss_w[1][2][r];
	  w13=gauss_w[1][3][r];
	  w21=gauss_w[2][1][r];
	  w22=gauss_w[2][2][r];
	  w23=gauss_w[2][3][r];
	  w31=gauss_w[3][1][r];
	  w32=gauss_w[3][2][r];
	  w33=gauss_w[3][3][r];
	  t11 = w11/gauss_rho2[1][r];        /* T := rho2^(-1) * W^T */
	  t12 = w21/gauss_rho2[1][r];
	  t13 = w31/gauss_rho2[1][r];
	  t21 = w12/gauss_rho2[2][r];
	  t22 = w22/gauss_rho2[2][r];
	  t23 = w32/gauss_rho2[2][r];
	  t31 = w13/gauss_rho2[3][r];
	  t32 = w23/gauss_rho2[3][r];
	  t33 = w33/gauss_rho2[3][r];
          a11 = w11*t11 + w12*t21 + w13*t31;  /* A := W * T    ( == W * rho2^(-1) * W^T)  */
          a12 = w11*t12 + w12*t22 + w13*t32;
          a13 = w11*t13 + w12*t23 + w13*t33;
          a21 = w21*t11 + w22*t21 + w23*t31;
          a22 = w21*t12 + w22*t22 + w23*t32;
          a23 = w21*t13 + w22*t23 + w23*t33;
          a31 = w31*t11 + w32*t21 + w33*t31;
          a32 = w31*t12 + w32*t22 + w33*t32;
          a33 = w31*t13 + w32*t23 + w33*t33;
          v1 = a11*dx + a12*dy + a13*dz;      /* v := A*x */
          v2 = a21*dx + a22*dy + a23*dz;
          v3 = a31*dx + a32*dy + a33*dz;
	  prod = dx*v1 + dy*v2 + dz*v3;       /* prod := x^T * v   ( == x^T * A * x )  */
	  if (prod<100.0) {
	    m11 = -a11 + v1*v1;               /* M := -A + (Ax)(Ax)^T */
	    m12 = -a12 + v1*v2;
	    m13 = -a13 + v1*v3;
	    m21 = -a21 + v2*v1;
	    m22 = -a22 + v2*v2;
	    m23 = -a23 + v2*v3;
	    m31 = -a31 + v3*v1;
	    m32 = -a32 + v3*v2;
	    m33 = -a33 + v3*v3;
	    ex = gauss_height[r]*exp(-prod/2.0);
	    sum += ex;
	    *gx -= (v1*ex);
	    *gy -= (v2*ex);
	    *gz -= (v3*ex);
	    cxx += (m11*ex);
	    cxy += (m12*ex);
	    cxz += (m13*ex);
	    cyx += (m21*ex);
	    cyy += (m22*ex);
	    cyz += (m23*ex);
	    czx += (m31*ex);
	    czy += (m32*ex);
	    czz += (m33*ex);
	    }
	  }
	*g =sqrt((*gx)*(*gx) + (*gy)*(*gy) + (*gz)*(*gz));
	if ((*g)!=0.0) {
	  g0x = (*gx)/(*g);
	  g0y = (*gy)/(*g);
	  g0z = (*gz)/(*g);
	  }
	else {
	  g0x = 0.0;
	  g0y = 0.0;
	  g0z = 0.0;
	  *gx = 1.0e-50;
	  *gy = 0.0;
	  *gz = 0.0;
	  *g  = 1.0e-50;
	  }
        v1 = cxx*g0x + cxy*g0y + cxz*g0z;
        v2 = cyx*g0x + cyy*g0y + cyz*g0z;
        v3 = czx*g0x + czy*g0y + czz*g0z;
	*cr = v1*g0x + v2*g0y + v3*g0z;
	temp=sqrt((*gx)*(*gx)+(*gy)*(*gy)+(*gz)*(*gz));
	if (temp!=0.0) {
	  *curvature=-0.5*((*gx)*(*gx)*(cyy+czz) + (*gy)*(*gy)*(cxx+czz) + (*gz)*(*gz)*(cxx+cyy)
	  	        -2.0*((*gx)*(*gy)*cxy + (*gx)*(*gz)*cxz + (*gy)*(*gz)*cyz))/(temp*temp*temp);
	  }
	else {
	  *curvature = 0.0;
	  }
	return (sum);
}


static double density_gradient(double x, double y, double z, double *gx, double *gy, double *gz)
{
	int r;
        double sum,di;
	double w11,w12,w13;
	double w21,w22,w23;
	double w31,w32,w33;
	double xt1,xt2,xt3,xtt1,xtt2,xtt3;
	double dx,dy,dz,ex;

	sum=0.0;
	*gx = 0.0;
	*gy = 0.0;
	*gz = 0.0;
	for (r=1; r<=ncodebook; r++) {
          dx=x-gauss_y[1][r];
          dy=y-gauss_y[2][r];
          dz=z-gauss_y[3][r];
	  w11=gauss_w[1][1][r];
	  w12=gauss_w[1][2][r];
	  w13=gauss_w[1][3][r];
	  w21=gauss_w[2][1][r];
	  w22=gauss_w[2][2][r];
	  w23=gauss_w[2][3][r];
	  w31=gauss_w[3][1][r];
	  w32=gauss_w[3][2][r];
	  w33=gauss_w[3][3][r];
	  xt1 = (w11*dx + w21*dy + w31*dz)/gauss_rho2[1][r];
	  xt2 = (w12*dx + w22*dy + w32*dz)/gauss_rho2[2][r];
	  xt3 = (w13*dx + w23*dy + w33*dz)/gauss_rho2[3][r];
	  xtt1 = (w11*xt1) + (w12*xt2) + (w13*xt3);
	  xtt2 = (w21*xt1) + (w22*xt2) + (w23*xt3);
	  xtt3 = (w31*xt1) + (w32*xt2) + (w33*xt3);
	  di=0.5*(dx*xtt1 + dy*xtt2 + dz*xtt3);
	  if (di<100.0) {
	    ex = gauss_height[r]*exp(-di);
	    sum += ex;
	    *gx -= xtt1*ex;
	    *gy -= xtt2*ex;
	    *gz -= xtt3*ex;
	    }
	  }
	return (sum);
}




static double get_min_distance(double rhofac, double level)
{
	int i;
	double dx,dy,dz,rr,rrr;
	double x0,y0,z0;
	double x1,y1,z1;
	double x2,y2,z2;
	double x3,y3,z3;
	double dens1,dens2,dens3;
	double dmin;
	double theta,thetastep;
	double phi,phistep;
	double cosphi,sinphi;

	dmin=1.0e20;
        thetastep=2.0*PI/(1.0+floor(2.0*PI/
		  atan(RESOLUTION_FOR_GET_MIN_DISTANCE/r_sphere)));
	for (theta=thetastep/2.0; theta<PI-thetastep/2.0+thetastep*0.1; theta+=thetastep) {
	  x0=sin(theta)*SCALE_FACTOR_FOR_GET_MIN_DISTANCE*(r_sphere+2.0*water_shell_thickness);
	  y0=0.0;
	  z0=cos(theta)*SCALE_FACTOR_FOR_GET_MIN_DISTANCE*(r_sphere+2.0*water_shell_thickness);
	  rrr=r_sphere*sin(theta);
          phistep=2.0*PI/(1.0+floor(2.0*PI/atan(RESOLUTION_FOR_GET_MIN_DISTANCE/rrr)));
	  for (phi=0.0; phi<2.0*PI-phistep*0.1; phi+=phistep) {
	    cosphi=cos(phi);
	    sinphi=sin(phi);
	    x1=x_sphere + cosphi*x0 - sinphi*y0;
	    y1=y_sphere + sinphi*x0 + cosphi*y0;
	    z1=z_sphere + z0;
	    x3=x_sphere;
	    y3=y_sphere;
	    z3=z_sphere;
	    dens1=density_varrho(x1,y1,z1,rhofac);
	    dens3=density_varrho(x3,y3,z3,rhofac);
	    if (dens1>level) {
	      fprintf(stderr,"\n\nERROR: dens1>level (%f > %f)\n",dens1,level);
	      fprintf(stderr,"       Increasing SCALE_FACTOR_FOR_GET_MIN_DISTANCE should help.\n");
	      terminate();
	      }
	    if (dens3<level) {
	      return(0.0);   /* to enforce an increase of rho_factor */
	      }
	    else {
	      do {
                dx=x3-x1;
                dy=y3-y1;
                dz=z3-z1;
	        rr=dx*dx+dy*dy+dz*dz;
	        x2=(x1+x3)/2.0;
	        y2=(y1+y3)/2.0;
	        z2=(z1+z3)/2.0;
	        dens2=density_varrho(x2,y2,z2,rhofac);
	        if (dens2>level) {
	          dens3=dens2;
	          x3=x2;
	          y3=y2;
	          z3=z2;
	          }
                else {
	          dens1=dens2;
	          x1=x2;
	          y1=y2;
	          z1=z2;
	          }
	        } while (rr>ACCURACY_OF_GET_MIN_DISTANCE*ACCURACY_OF_GET_MIN_DISTANCE);
	      x2=(x1+x3)/2.0;
	      y2=(y1+y3)/2.0;
	      z2=(z1+z3)/2.0;
	      for (i=1;i<=natoms;i++) {
                dx=x2-xatom[i];
                dy=y2-yatom[i];
                dz=z2-zatom[i];
	        rr=dx*dx+dy*dy+dz*dz;
	        if (rr<dmin) dmin=rr;
	        }
	      }
	    }
	  }
        return(sqrt(dmin));
}



static void adjust_surface_distance(void)
{
	int i;
	double volaver,rhofac,rhofac_try,rhofac_shrink_factor;
	double d_now,d_try;

	fprintf(stderr,"\nAdjusting boundary distance by varying gaussian widths ...          ");
	volaver=0.0;
	gauss_vol=dvector(1,ncodebook);
	gauss_height=dvector(1,ncodebook);
	for (i=1; i<=ncodebook; i++) {
	  gauss_vol[i]=pow((gauss_rho2[1][i]*gauss_rho2[2][i]*gauss_rho2[3][i]),1.0/(2.0*3.0));
	  volaver += gauss_vol[i];
	  }
	volaver /= (1.0*ncodebook);
	rhofac=INITIAL_RHOFAC;
        d_now = get_min_distance(rhofac,SURFACE_LEVEL/volaver);
	if (water_shell_thickness > d_now) {
	  fprintf(stderr,"\n\nERROR: water_shell_thickness > d_now (%f > %f)\n",water_shell_thickness,d_now);
	  fprintf(stderr,"       Increasing INITIAL_RHOFAC should help.\n");
	  terminate();
	  }
	rhofac_shrink_factor = INITIAL_RHOFAC_SHRINK_FACTOR;
	do {
	  rhofac_try=rhofac*rhofac_shrink_factor;
	  d_try=get_min_distance(rhofac_try,SURFACE_LEVEL/volaver);
	  if (water_shell_thickness<d_try) {
	    rhofac=rhofac_try;
	    d_now=d_try;
	    }
          else {
            rhofac_shrink_factor = (5.0+3.0*rhofac_shrink_factor)/(5.0+3.0);
	    }
	  fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b%8.3f A",d_try);
	  } while (fabs(d_now-water_shell_thickness)>ACCURACY_OF_ADJUST_SURFACE_DISTANCE);
	fprintf(stderr,"\n");
	fprintf(stderr,"Gaussian widths are being scaled by a factor of %f\n",rhofac);
	for (i=1; i<=ncodebook; i++) {
	  gauss_rho2[1][i] *= (rhofac*rhofac);
	  gauss_rho2[2][i] *= (rhofac*rhofac);
	  gauss_rho2[3][i] *= (rhofac*rhofac);
	  gauss_height[i]   = volaver/(gauss_vol[i]);
	  }
        return;
}


static void save_gaussian_parameters()
{
	int r;
	FILE	*fp;

	if (ncodebook==1) fprintf(stderr,"Saving parameter for one gaussian to file gaussians.lis ... ");
	if (ncodebook!=1) fprintf(stderr,"Saving parameter for %d gaussians to file gaussians.lis ... ",ncodebook);
	if ((fp=fopen("gaussians.lis","w"))==NULL) {
	  fprintf(stderr,"\n\nERROR: Could not write to file gaussians.lis\n");
	  terminate();
	  }
	fprintf(fp,"[ENCLOSING SPHERE (X,Y,Z,R)]\n");
	fprintf(fp,"%20e\n",x_sphere);
	fprintf(fp,"%20e\n",y_sphere);
	fprintf(fp,"%20e\n",z_sphere);
	fprintf(fp,"%20e\n",r_sphere);
	fprintf(fp,"[NR OF GAUSSIANS]\n");
	fprintf(fp,"%d\n",ncodebook);
	fprintf(fp,"[SURFACE LEVEL]\n");
	fprintf(fp,"%20e\n",SURFACE_LEVEL);
	fprintf(fp,"[GAUSSIAN CENTERS]\n");
	for (r=1; r<=ncodebook; r++) {
	  fprintf(fp,"%20e %20e %20e\n",gauss_y[1][r],gauss_y[2][r],gauss_y[3][r]);
	  }
	fprintf(fp,"[GAUSSIAN VARIANCES]\n");
	for (r=1; r<=ncodebook; r++) {
	  fprintf(fp,"%20e %20e %20e\n",gauss_rho2[1][r],gauss_rho2[2][r],gauss_rho2[3][r]);
	  }
	fprintf(fp,"[GAUSSIAN ORIENTATIONS]\n");
	for (r=1; r<=ncodebook; r++) {
	  fprintf(fp,"%20e %20e %20e\n",gauss_w[1][1][r],gauss_w[1][2][r],gauss_w[1][3][r]);
	  fprintf(fp,"%20e %20e %20e\n",gauss_w[2][1][r],gauss_w[2][2][r],gauss_w[2][3][r]);
	  fprintf(fp,"%20e %20e %20e\n",gauss_w[3][1][r],gauss_w[3][2][r],gauss_w[3][3][r]);
	  }
	fclose(fp);
	fprintf(stderr,"done.\n");
	return;
}


static void read_gaussian_parameters()
{
	int r;
	float fdummy;
	FILE	*fp;

	fprintf(stderr,"Reading parameter for gaussians from file gaussians.lis ... ");
	if ((fp=fopen("gaussians.lis","r"))==NULL) {
	  fprintf(stderr,"\n\nERROR: Could not read from file gaussians.lis\n");
	  terminate();
	  }
	skip_line(fp);
	fscanf(fp,"%f",&fdummy);
	x_sphere=fdummy;
	fscanf(fp,"%f",&fdummy);
	y_sphere=fdummy;
	fscanf(fp,"%f",&fdummy);
	z_sphere=fdummy;
	fscanf(fp,"%f",&fdummy);
	r_sphere=fdummy;
	skip_line(fp);
	fscanf(fp,"%d",&ncodebook);
	if (ncodebook>10000) {
	  fprintf(stderr,"\n\nERROR: Found too many gaussians (%d).\n",ncodebook);
	  fprintf(stderr,"       Probably file gaussians.lis is corrupted.\n");
	  terminate();
	  }
	if (ncodebook <1) {
	  fprintf(stderr,"\n\nERROR: Found too few gaussians (%d).\n",ncodebook);
	  fprintf(stderr,"       Probably file gaussians.lis is corrupted.\n");
	  terminate();
	  }
	skip_line(fp);
	fscanf(fp,"%f",&fdummy);
	if (fabs(fdummy-SURFACE_LEVEL)>1.0e-5) {
	  fprintf(stderr,"\n\nERROR: The surface level found in gaussians.lis (%f)\n",fdummy);
	  fprintf(stderr,"       differs from the one used here (%f).\n",SURFACE_LEVEL);
	  fprintf(stderr,"       Probably file gaussians.lis is corrupted.\n");
	  terminate();
	  }
	skip_line(fp);
	gauss_y = dmatrix(1,3,1,ncodebook);
	gauss_w = dcube(1,3,1,3,1,ncodebook);
	gauss_rho2 = dmatrix(1,3,1,ncodebook);
	for (r=1; r<=ncodebook; r++) {
	  fscanf(fp,"%f",&fdummy);
	  gauss_y[1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_y[2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_y[3][r]=fdummy;
	  }
	skip_line(fp);
	for (r=1; r<=ncodebook; r++) {
	  fscanf(fp,"%f",&fdummy);
	  gauss_rho2[1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_rho2[2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_rho2[3][r]=fdummy;
	  }
	skip_line(fp);
	for (r=1; r<=ncodebook; r++) {
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[1][1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[1][2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[1][3][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[2][1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[2][2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[2][3][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[3][1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[3][2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[3][3][r]=fdummy;
	  }
	fclose(fp);
	if (ncodebook==1) fprintf(stderr,"found %d gaussian.\n",ncodebook);
	if (ncodebook!=1) fprintf(stderr,"found %d gaussians.\n",ncodebook);
	return;
}



static void read_boundary_parameters()
{
	int r;
	float fdummy;
	FILE	*fp;

	fprintf(stderr,"Reading boundary parameters from file boundary.lis ... ");
	if ((fp=fopen("boundary.lis","r"))==NULL) {
	  fprintf(stderr,"\n\nERROR: Could not read from file boundary.lis\n");
	  terminate();
	  }
	skip_line(fp);
	fscanf(fp,"%f",&fdummy);
	water_shell_thickness=fdummy;
	skip_line(fp);
	fscanf(fp,"%f",&fdummy);
	x_sphere=fdummy;
	fscanf(fp,"%f",&fdummy);
	y_sphere=fdummy;
	fscanf(fp,"%f",&fdummy);
	z_sphere=fdummy;
	fscanf(fp,"%f",&fdummy);
	r_sphere=fdummy;
	skip_line(fp);
	fscanf(fp,"%d",&ncodebook);
	if (ncodebook>10000) {
	  fprintf(stderr,"\n\nERROR: Found too many gaussians (%d).\n",ncodebook);
	  fprintf(stderr,"       Probably file boundary.lis is corrupted.\n");
	  terminate();
	  }
	if (ncodebook <1) {
	  fprintf(stderr,"\n\nERROR: Found too few gaussians (%d).\n",ncodebook);
	  fprintf(stderr,"       Probably file boundary.lis is corrupted.\n");
	  terminate();
	  }
	skip_line(fp);
	fscanf(fp,"%f",&fdummy);
	if (fabs(fdummy-SURFACE_LEVEL)>1.0e-5) {
	  fprintf(stderr,"\n\nERROR: The surface level found in boundary.lis (%f)\n",fdummy);
	  fprintf(stderr,"       differs from the one used here (%f).\n",SURFACE_LEVEL);
	  fprintf(stderr,"       Probably file boundary.lis is corrupted.\n");
	  terminate();
	  }
	skip_line(fp);
	gauss_y = dmatrix(1,3,1,ncodebook);
	gauss_w = dcube(1,3,1,3,1,ncodebook);
	gauss_rho2 = dmatrix(1,3,1,ncodebook);
	gauss_height = dvector(1,ncodebook);
	for (r=1; r<=ncodebook; r++) {
	  fscanf(fp,"%f",&fdummy);
	  gauss_height[r]=fdummy;
	  }
	skip_line(fp);
	for (r=1; r<=ncodebook; r++) {
	  fscanf(fp,"%f",&fdummy);
	  gauss_y[1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_y[2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_y[3][r]=fdummy;
	  }
	skip_line(fp);
	for (r=1; r<=ncodebook; r++) {
	  fscanf(fp,"%f",&fdummy);
	  gauss_rho2[1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_rho2[2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_rho2[3][r]=fdummy;
	  }
	skip_line(fp);
	for (r=1; r<=ncodebook; r++) {
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[1][1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[1][2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[1][3][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[2][1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[2][2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[2][3][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[3][1][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[3][2][r]=fdummy;
	  fscanf(fp,"%f",&fdummy);
	  gauss_w[3][3][r]=fdummy;
	  }
	fclose(fp);
	if (ncodebook==1) fprintf(stderr,"\n        found %d gaussian;\n",ncodebook);
	if (ncodebook!=1) fprintf(stderr,"\n        found %d gaussians;\n",ncodebook);
	fprintf(stderr,"        minimum water shell thickness is %5.2f A.\n",water_shell_thickness);
	return;
}




static void save_boundary_parameters()
{
	int r;
	FILE	*fp;

	fprintf(stderr,"Saving boundary parameters to file boundary.lis ... ");
	if ((fp=fopen("boundary.lis","w"))==NULL) {
	  fprintf(stderr,"\n\nERROR: Could not write to file boundary.lis\n");
	  terminate();
	  }
        fprintf(fp,"[MINIMUM WATER SHELL THICKNESS]\n");
        fprintf(fp,"%20e\n",water_shell_thickness);
        fprintf(fp,"[ENCLOSING SPHERE (X,Y,Z,R)]\n");
        fprintf(fp,"%20e\n",x_sphere);
        fprintf(fp,"%20e\n",y_sphere);
        fprintf(fp,"%20e\n",z_sphere);
        fprintf(fp,"%20e\n",r_sphere);
	fprintf(fp,"[NR OF GAUSSIANS]\n");
	fprintf(fp,"%d\n",ncodebook);
	fprintf(fp,"[SURFACE LEVEL]\n");
	fprintf(fp,"%20e\n",SURFACE_LEVEL);
	fprintf(fp,"[GAUSSIAN HEIGHTS]\n");
	for (r=1; r<=ncodebook; r++) {
	  fprintf(fp,"%20e\n",gauss_height[r]);
	  }
	fprintf(fp,"[GAUSSIAN CENTERS]\n");
	for (r=1; r<=ncodebook; r++) {
	  fprintf(fp,"%20e %20e %20e\n",gauss_y[1][r],gauss_y[2][r],gauss_y[3][r]);
	  }
	fprintf(fp,"[GAUSSIAN VARIANCES]\n");
	for (r=1; r<=ncodebook; r++) {
	  fprintf(fp,"%20e %20e %20e\n",gauss_rho2[1][r],gauss_rho2[2][r],gauss_rho2[3][r]);
	  }
	fprintf(fp,"[GAUSSIAN ORIENTATIONS]\n");
	for (r=1; r<=ncodebook; r++) {
	  fprintf(fp,"%20e %20e %20e\n",gauss_w[1][1][r],gauss_w[1][2][r],gauss_w[1][3][r]);
	  fprintf(fp,"%20e %20e %20e\n",gauss_w[2][1][r],gauss_w[2][2][r],gauss_w[2][3][r]);
	  fprintf(fp,"%20e %20e %20e\n",gauss_w[3][1][r],gauss_w[3][2][r],gauss_w[3][3][r]);
	  }
	fclose(fp);
	fprintf(stderr,"done.\n");
	return;
}






static void create_volume(void)
{
	long int ii,i;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double rr,rr2,dmin,dd;
	double x,y,z;
	double dx,dy,dz,ddd;

	fprintf(stderr,"\nCreating volume within boundary ... ");
        xmin=floor(x_sphere-r_sphere-water_shell_thickness-ADD_LAYER_FOR_CREATE_VOLUME);
        xmax=floor(x_sphere+r_sphere+water_shell_thickness+ADD_LAYER_FOR_CREATE_VOLUME+1.0);
        ymin=floor(y_sphere-r_sphere-water_shell_thickness-ADD_LAYER_FOR_CREATE_VOLUME);
        ymax=floor(y_sphere+r_sphere+water_shell_thickness+ADD_LAYER_FOR_CREATE_VOLUME+1.0);
        zmin=floor(z_sphere-r_sphere-water_shell_thickness-ADD_LAYER_FOR_CREATE_VOLUME);
        zmax=floor(z_sphere+r_sphere+water_shell_thickness+ADD_LAYER_FOR_CREATE_VOLUME+1.0);
	rr=r_sphere+water_shell_thickness+ADD_LAYER_FOR_CREATE_VOLUME+1.0;
	rr2=rr*rr;
	ngrid=0;
	for (x=xmin; x<xmax+0.1; x++) {
	  for (y=ymin; y<ymax+0.1; y++) {
	    for (z=zmin; z<zmax+0.1; z++) {
	      dx=x-x_sphere;
	      dy=y-y_sphere;
	      dz=z-z_sphere;
	      if (dx*dx+dy*dy+dz*dz<=rr2) {
		if (density(x,y,z)>=SURFACE_LEVEL) ngrid++;
		}
              }
	    }
	  }
	xgrid=dvector(1,ngrid);
	ygrid=dvector(1,ngrid);
	zgrid=dvector(1,ngrid);
	dgrid=dvector(1,ngrid);
	densgrid=dvector(1,ngrid);
	ii=1;
	for (x=xmin; x<xmax+0.1; x++) {
	  for (y=ymin; y<ymax+0.1; y++) {
	    for (z=zmin; z<zmax+0.1; z++) {
	      dx=x-x_sphere;
	      dy=y-y_sphere;
	      dz=z-z_sphere;
	      if (dx*dx+dy*dy+dz*dz<=rr2) {
		ddd=density(x,y,z);
		if (ddd>=SURFACE_LEVEL) {
		  xgrid[ii]=x;
		  ygrid[ii]=y;
		  zgrid[ii]=z;
		  densgrid[ii]=ddd;
		  dmin=1.0e20;
		  for (i=1; i<=natoms; i++) {
		    dx=x-xatom[i];
		    dy=y-yatom[i];
		    dz=z-zatom[i];
		    dd=dx*dx+dy*dy+dz*dz;
		    if (dd<dmin) dmin=dd;
		    }
		  dgrid[ii]=sqrt(dmin);
		  ii++;
		  }
		}
              }
	    }
	  }
	fprintf(stderr,"done.\n");
	fprintf(stderr,"Sorting grid points according to distance from solute ... ");
	sort5(ngrid,dgrid,xgrid,ygrid,zgrid,densgrid);
	fprintf(stderr,"done.\n");
	return;
}





static void surface_statistics(void)
{
	int count,countmax;
	long int i;
	double gx,gy,gz,gx0,gy0,gz0,dd,dd0,dx,dy,dz,x0,y0,z0,x1,y1,z1,r;
	double initial_stepsize,stepsize;
	double xtry,ytry,ztry;
	double gxnew,gynew,gznew,ddnew;
	double estimated_distance,true_distance;
	double gr,cur,s2,kk,xx0,yy0,yy1,yy2,err;
	double ymin1,ymin2,ymin3,ymin4,ymin5,ymin6,ymin7;
	double dmin1,dmin2,dmin3,dmin4,dmin5,dmin6,dmin7;
	double curvature;
	double dist_0,dist_1,dens_min,dd_0;
	FILE *fp;

	fprintf(stderr,"Performing distance approximation statistics with %ld grid points ...     ",ngrid);
	sdistgrid=dvector(1,ngrid);
	distgrid=dvector(1,ngrid);
	gradgrid=dvector(1,ngrid);
	curvegrid=dvector(1,ngrid);
	curvaturegrid=dvector(1,ngrid);

	/* (for benchmark) */
	/*
	fprintf(stderr,"\n\nstart benchmark.\n\n");
	system("date");
	for (i=1; i<=ngrid; i++) {
	  x0=xgrid[i];
	  y0=ygrid[i];
	  z0=zgrid[i];
	  dd=density_grad_curve(x0,y0,z0,&gx,&gy,&gz,&gr,&cur,&curvature);
	  yy0 = dd;
	  yy1 = gr;
	  yy2 = cur;
	  s2  = yy0*yy0/(yy1*yy1-yy0*yy2);
	  kk  = yy0*exp(0.5*yy1*yy1/(yy1*yy1-yy0*yy2));
	  xx0 = yy0*yy1/(yy1*yy1-yy0*yy2);
	  estimated_distance = fabs(xx0 - sqrt(2.0*s2*log(kk/SURFACE_LEVEL)));
	  }
	system("date");
	fprintf(stderr,"\nend benchmark.\n\n");
	*/

	for (i=1; i<=ngrid; i++) {
	  if (floor(((i-1.0)*100.0)/(1.0*ngrid))!=floor((i*100.0)/(1.0*ngrid))) {
	    fprintf(stderr,"\b\b\b\b%3d%%",(int)(floor((i*100.0)/(1.0*ngrid))));
	    }
	  x0=xgrid[i];
	  y0=ygrid[i];
	  z0=zgrid[i];
	  dd=density_grad_curve(x0,y0,z0,&gx,&gy,&gz,&gr,&cur,&curvature);
	  dd0=dd;
	  gx0=gx;
	  gy0=gy;
	  gz0=gz;
          initial_stepsize=INITIAL_STEPSIZE_FOR_SURFACE_STATISTICS;
          stepsize=initial_stepsize;
          countmax=2*floor(water_shell_thickness/initial_stepsize);
          while (stepsize>ACCURACY_FOR_SURFACE_STATISTICS) {
            count=0;
            do {
              count++;
              r=sqrt(gx*gx+gy*gy+gz*gz);
              xtry = x0 - (stepsize*gx/r);
              ytry = y0 - (stepsize*gy/r);
              ztry = z0 - (stepsize*gz/r);
              ddnew = density_gradient(xtry,ytry,ztry,&gxnew,&gynew,&gznew);
              if (ddnew >= SURFACE_LEVEL) {
                x0=xtry;
                y0=ytry;
                z0=ztry;
                gx=gxnew;
                gy=gynew;
                gz=gznew;
                dd=ddnew;
                }
              } while ((ddnew >= SURFACE_LEVEL) && (count<countmax));
            stepsize /= 2.0;
	    countmax=3;
            }
	  dx=x0-xgrid[i];
	  dy=y0-ygrid[i];
	  dz=z0-zgrid[i];
	  densgrid[i]=dd0;
	  gradgrid[i]=sqrt(gx0*gx0+gy0*gy0+gz0*gz0);
	  curvegrid[i]=cur;
	  curvaturegrid[i]=curvature;
	  sdistgrid[i]=sqrt(dx*dx+dy*dy+dz*dz);
	  if ((sdistgrid[i]>DISTANCE_TIMES_WATER_SHELL_THICKNESS*water_shell_thickness) || (fabs(ddnew-SURFACE_LEVEL)>SURFACE_LEVEL/1000.0)) {
	    sdistgrid[i]=DISTANCE_TIMES_WATER_SHELL_THICKNESS*water_shell_thickness;
	    }
	  }
	fprintf(stderr,"\n");
	if (write_surface_statistics_flag==1) {
	  if ((fp=fopen("surface_stat.lis","w"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not write to file surface_stat.lis\n");
	    terminate();
	    }
	  fprintf(stderr,"Writing boundary coordinates to file surface_stat.lis ... ");
	  fprintf(fp,"[SURFACE COORDINATES (x,y,z,curvature)]\n");
	  for (i=1; i<=ngrid; i++) {
	    if (densgrid[i]<(1.0+BOUNDARY_SURFACE_THICKNESS_FACTOR)*SURFACE_LEVEL) {
              fprintf(fp,"%f %f %f  %f\n",xgrid[i],ygrid[i],zgrid[i],curvaturegrid[i]);
	      }
	    }
	  fclose(fp);
	  fprintf(stderr,"done.\n");
	  }
	/* Now let's check whether our distance computation via  */
	/*   densgrid[],gradgrid[], and curvegrid[] is correct: */
	fprintf(stderr,"Checking approximation for distance from boundary ...     ");
	ymin1=1.0e20;
	ymin2=1.0e20;
	ymin3=1.0e20;
	ymin4=1.0e20;
	ymin5=1.0e20;
	ymin6=1.0e20;
	ymin7=1.0e20;
	dens_min=1.0e10;
	for (i=1; i<=ngrid; i++) {
	  if (floor(((i-1.0)*100.0)/(1.0*ngrid))!=floor((i*100.0)/(1.0*ngrid))) {
	    fprintf(stderr,"\b\b\b\b%3d%%",(int)(floor((i*100.0)/(1.0*ngrid))));
	    }
	  x0=xgrid[i];
	  y0=ygrid[i];
	  z0=zgrid[i];
	  dd=density_grad_curve(x0,y0,z0,&gx,&gy,&gz,&gr,&cur,&curvature);
	  yy0 = dd;
	  dd_0 = dd;
	  yy1 = sqrt(gx*gx+gy*gy+gz*gz);
	  yy2 = cur;
	  s2  = yy0*yy0/(yy1*yy1-yy0*yy2);
	  kk  = yy0*exp(0.5*yy1*yy1/(yy1*yy1-yy0*yy2));
	  xx0 = yy0*yy1/(yy1*yy1-yy0*yy2);
	  dist_0 = fabs(xx0 - sqrt(2.0*s2*log(kk/SURFACE_LEVEL)));
	  if (sdistgrid[i]<DISTANCE_TIMES_WATER_SHELL_THICKNESS*water_shell_thickness) {
            x1 = x0 - dist_0*gx/yy1;
            y1 = y0 - dist_0*gy/yy1;
            z1 = z0 - dist_0*gz/yy1;
	    dd=density_grad_curve(x1,y1,z1,&gx,&gy,&gz,&gr,&cur,&curvature);
	    curvaturegrid[i]=curvature;
	    yy0 = dd;
	    yy1 = sqrt(gx*gx+gy*gy+gz*gz);
	    yy2 = cur;
	    s2  = yy0*yy0/(yy1*yy1-yy0*yy2);
	    kk  = yy0*exp(0.5*yy1*yy1/(yy1*yy1-yy0*yy2));
	    xx0 = yy0*yy1/(yy1*yy1-yy0*yy2);
	    dist_1 = fabs(xx0 - sqrt(2.0*s2*log(kk/SURFACE_LEVEL)));
	    estimated_distance = dist_0 + dist_1;
	    distgrid[i]=estimated_distance;
	    true_distance = sdistgrid[i];
	    err=fabs(true_distance-estimated_distance);
	    if (err>0.01) {
	      if (dd_0<ymin1) {
	        ymin1=dd_0;
	        dmin1=true_distance;
	        }
	      }
	    if (err>0.02) {
	      if (dd_0<ymin2) {
	        ymin2=dd_0;
	        dmin2=true_distance;
	        }
	      }
	    if (err>0.05) {
	      if (dd_0<ymin3) {
	        ymin3=dd_0;
	        dmin3=true_distance;
	        }
	      }
	    if (err>0.1) {
	      if (dd_0<ymin4) {
	        ymin4=dd_0;
	        dmin4=true_distance;
	        }
	      }
	    if (err>0.2) {
	      if (dd_0<ymin5) {
	        ymin5=dd_0;
	        dmin5=true_distance;
	        }
	      }
	    if (err>0.5) {
	      if (dd_0<ymin6) {
	        ymin6=dd_0;
	        dmin6=true_distance;
	        }
	      }
	    if (err>1.0) {
	      if (dd_0<ymin7) {
	        ymin7=dd_0;
	        dmin7=true_distance;
		}
	      }
	    }
	  else {
	    curvaturegrid[i]=0.0;
            distgrid[i]=10.0;
	    if (dens_min > dd) dens_min=dd;
	    }
	  }
	fprintf(stderr,"\n");
	if (write_surface_statistics_flag==1) {
	  if ((fp=fopen("surface_stat.lis","a"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not append to file surface_stat.lis\n");
	    terminate();
	    }
	  fprintf(fp,"[MINIMUM INVALID DENSITY]\n");
	  fprintf(fp,"%f\n",dens_min);
	  fprintf(fp,"[DISTANCE ERROR STATISTICS (ABS. ERR / DENSITY / DISTANCE)]\n");
	  fprintf(fp,"0.01  %f %f\n",ymin1,dmin1);
	  fprintf(fp,"0.02  %f %f\n",ymin2,dmin2);
	  fprintf(fp,"0.05  %f %f\n",ymin3,dmin3);
	  fprintf(fp,"0.10  %f %f\n",ymin4,dmin4);
	  fprintf(fp,"0.20  %f %f\n",ymin5,dmin5);
	  fprintf(fp,"0.50  %f %f\n",ymin6,dmin6);
	  fprintf(fp,"1.00  %f %f\n",ymin7,dmin7);
	  fclose(fp);
	  }
	if (write_volume_statistics_flag==1) {
	  fprintf(stderr,"Writing grid points within boundary to file volume_stat.lis ... ");
	  if ((fp=fopen("volume_stat.lis","w"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not write to file volume_stat.lis\n");
	    terminate();
	    }
	  for (i=1; i<=ngrid; i++) {
            fprintf(fp,"%f %f %f %f %f %f %f\n",xgrid[i],ygrid[i],zgrid[i],densgrid[i],sdistgrid[i],distgrid[i],curvaturegrid[i]);
	    }
	  fclose(fp);
	  fprintf(stderr,"done.\n");
	  }
	free_dvector(sdistgrid,1,ngrid);
	free_dvector(gradgrid,1,ngrid);
	free_dvector(curvegrid,1,ngrid);
	free_dvector(curvaturegrid,1,ngrid);
	return;
}




static double vdw_energy_force(double x, double y, double z, double *fx, double *fy, double *fz)
{
    double evdw,dx,dy,dz,dd,fvdw,r_2,r_6;
    int i;

    *fx = 0.0;
    *fy = 0.0;
    *fz = 0.0;
    evdw=0.0;
    for (i=1; i<=natoms; i++) {
      dx=x-xatom[i];
      dy=y-yatom[i];
      dz=z-zatom[i];
      dd=dx*dx+dy*dy+dz*dz;
      if (dd<=100.0) {
	r_2=1.0/dd;
        r_6=r_2*r_2*r_2;
        evdw += (r_6* (r_6*a_vdw_atom[i]-b_vdw_atom[i]));
        fvdw = 6.0 * r_2 * r_6* (2.0*r_6*a_vdw_atom[i]-b_vdw_atom[i]);
        *fx += (fvdw * dx);
        *fy += (fvdw * dy);
        *fz += (fvdw * dz);
	}
      }
    if (nwater>0) {
      for (i=1; i<=nwater; i++) {
        dx=x-xwater[i];
        dy=y-ywater[i];
        dz=z-zwater[i];
        dd=dx*dx+dy*dy+dz*dz;
        if (dd<=100.0) {
  	  r_2=1.0/dd;
          r_6=r_2*r_2*r_2;
          evdw += (r_6* (r_6*a_vdw-b_vdw));
          fvdw = 6.0 * r_2 * r_6* (2.0*r_6*a_vdw-b_vdw);
          *fx += (fvdw * dx);
          *fy += (fvdw * dy);
          *fz += (fvdw * dz);
  	  }
        }
      }
    return(evdw);
}


static void place_water_molecules(void)
{
        long int ii,jj,i;
	double sig,eps;
	double sig_6,sig_12;
	double dx,dy,dz;
	double xx,yy,zz;
	double r_min,r_aver,fx,fy,fz,evdw,stepsize,dd;
	double xxnew,yynew,zznew,evdwnew,fxnew,fynew,fznew,dummy;

	sig=SIG_H2O;
	eps=0.15;
        sig_6   = sig*sig*sig*sig*sig*sig;
        sig_12  = sig_6*sig_6;
        a_vdw   = 4.0*sig_12*eps;
        b_vdw   = 4.0*sig_6*eps;
	r_min=sig-0.5;   /* 0.5: half grid spacing */
	r_aver=sig*pow(2.0,1.0/6.0);
	a_vdw_atom=dvector(1,natoms);
	b_vdw_atom=dvector(1,natoms);
        for (i=1; i<=natoms; i++) {
	  sig = 0.5*(sigatom[i]+SIG_H2O);
	  eps = sqrt(0.15*epsatom[i]);
	  sig_6 = sig*sig*sig*sig*sig*sig;
	  sig_12 = sig_6*sig_6;
	  a_vdw_atom[i] = 4.0*sig_12*eps;
	  b_vdw_atom[i] = 4.0*sig_6*eps;
	  }
	/* estimate number of water molecules: */
	maxnwater=0;
	for (ii=1; ii<=ngrid; ii++) {
	  if (dgrid[ii]>r_min) maxnwater++;
	  }
	maxnwater=floor((MAXNWATER_FACTOR*maxnwater)/(4.0/3.0*PI*r_aver*r_aver*r_aver/(2.0*2.0*2.0)))+1;
	fprintf(stderr,"\nEstimated number of water molecules: %ld\n",maxnwater);
	fprintf(stderr,"Placing water molecules ...        ");

	xwater=dvector(1,maxnwater);
	ywater=dvector(1,maxnwater);
	zwater=dvector(1,maxnwater);
	replace_water_by_ion=ivector(1,maxnwater);
	for (ii=1; ii<=maxnwater; ii++) replace_water_by_ion[ii]=0;
	nwater=0;

	for (ii=1; ii<=ngrid; ii++) {
	  if (dgrid[ii]>r_min) {
            xx=xgrid[ii];
            yy=ygrid[ii];
            zz=zgrid[ii];
            stepsize=INITIAL_STEPSIZE_FOR_PLACE_WATER;
            evdw = vdw_energy_force(xx, yy, zz, &fx, &fy, &fz);
            dd=sqrt(fx*fx+fy*fy+fz*fz);
	    fx /= dd;
	    fy /= dd;
	    fz /= dd;
	    xxnew=xx+stepsize*fx;
	    yynew=yy+stepsize*fy;
	    zznew=zz+stepsize*fz;
	    do {
              evdwnew = vdw_energy_force(xxnew, yynew, zznew, &fxnew, &fynew, &fznew);
	      if (evdwnew<evdw) {
		evdw=evdwnew;
                dd=sqrt(fxnew*fxnew+fynew*fynew+fznew*fznew);
		fx=fxnew/dd;
		fy=fynew/dd;
		fz=fznew/dd;
		xx=xxnew;
		yy=yynew;
		zz=zznew;
		}
	      else {
		stepsize /= 2.0;
		}
	      xxnew=xx+stepsize*fx;
	      yynew=yy+stepsize*fy;
	      zznew=zz+stepsize*fz;
	      } while (stepsize>FINAL_STEPSIZE_FOR_PLACE_WATER);
	    if (evdw<=0.0) {
              dd = density_gradient(xx,yy,zz,&dummy,&dummy,&dummy);
              if (dd >= SURFACE_LEVEL) {
	        nwater++;
	        fprintf(stderr,"\b\b\b\b\b\b\b%-7ld",nwater);
	        if (nwater>=maxnwater) {
		  fprintf(stderr,"\n\nERROR: Too many water molecules (%ld)\n",nwater);
		  fprintf(stderr,"       Increase MAXNWATER_FACTOR!\n");
		  terminate();
		  }
	        xwater[nwater] = xx;
	        ywater[nwater] = yy;
	        zwater[nwater] = zz;
	        /* mask grid points around new water molecule: */
	        if (ii<ngrid) {
	          for (jj=ii+1; jj<=ngrid; jj++) {
		    dx=xgrid[jj]-xx;
		    dy=ygrid[jj]-yy;
		    dz=zgrid[jj]-zz;
		    if (dx*dx+dy*dy+dz*dz<=r_min*r_min) dgrid[jj]=0.0;
		    }
		  }
		}
	      }
	    }
	  }
	fprintf(stderr,"\b\b\b\b\b\b\bfound %ld\n",nwater);
	free_dvector(xgrid,1,ngrid);
	free_dvector(ygrid,1,ngrid);
	free_dvector(zgrid,1,ngrid);
	free_dvector(dgrid,1,ngrid);
	free_dvector(densgrid,1,ngrid);
	return;
}


/*


                  EPS      SIGMA

 NONBonded  C       0.1200   3.7418
 NONBonded  CH1E    0.0486   4.2140
 NONBonded  CH2E    0.1142   3.9823
 NONBonded  CH3E    0.1811   3.8576
 NONBonded  N       0.2384   2.8509
 NONBonded  O       0.1591   2.8509
 NONBonded  OC      0.6469   2.8509
 NONBonded  S       0.0430   3.3676


sig_ij  =(sigma[i]+sigma[j])/2.0;
eps_ij  =sqrt(eps[i]*eps[j]);

                  a_vdw        b_vdw
nbfix ot  ot  581980.4948  595.0436396
nbfix ht  ht  3.085665E-06 7.533363E-04
nbfix ht  ot  327.8404792  10.47230620

sig_6   = sig*sig*sig*sig*sig*sig;
sig_12  = sig_6*sig_6;
a_vdw   = 4.0*sig_12*eps;
b_vdw   = 4.0*sig_6*eps;
r_6     = r_2*r_2*r_2;
vdw_loc_energy = r_6 * (r_6 * a_vdw - b_vdw) / 2.0;
fvdw    = r_2 * (24.0*vdw_loc_energy + 6.0*r_6*b_vdw);
fx_h = fvdw * rx;
fy_h = fvdw * ry;
fz_h = fvdw * rz;


*/


static void group_water_molecules(void)
{
	long int i,j,count,n1,n2,n3;
	int foundflag,finishflag;
	double dx,dy,dz,x1,y1,z1,r_aver,rcontact2;
	FILE *fp;

	fprintf(stderr,"Grouping water molecules ...     ");
	r_aver=SIG_H2O*pow(2.0,1.0/6.0);
	rcontact2=r_aver*r_aver*1.2*1.2;
	groupflags = ivector(1,nwater);
	nr_of_water_groups=1;   /* 1: bulk water group */
	for (i=1; i<=nwater; i++) groupflags[i]=0;
	groupflags[nwater]=nr_of_water_groups;  /* (the outermost water defines group nr.one) */
	count=1;
	do {
	  do {
	    foundflag=0;
	    for (i=1; i<=nwater; i++) {
	      if (groupflags[i]==nr_of_water_groups) {
		x1=xwater[i];
		y1=ywater[i];
		z1=zwater[i];
	        for (j=1; j<=nwater; j++) {
	          if (groupflags[j]==0) {
		    dx = x1 - xwater[j];
		    dy = y1 - ywater[j];
		    dz = z1 - zwater[j];
		    if ((dx*dx + dy*dy + dz*dz)<=rcontact2) {
		      groupflags[j]=nr_of_water_groups;
                      foundflag=1;
		      count++;
	              if (floor(((count+1)*100.0+0.5)/(1.0*nwater))!=floor((count*100.0+0.5)/(1.0*nwater))) {
	                fprintf(stderr,"\b\b\b\b%3d%%",(int)(floor((count*100.0+1.0)/(1.0*nwater))));
	                }
		      }
		    }
		  }
		}
	      }

	    } while (foundflag==1);
	  finishflag=1;
          for (i=nwater; i>=1; i--) {
	    if (groupflags[i]==0) {
	      nr_of_water_groups++;
	      groupflags[i]=nr_of_water_groups;
	      finishflag=0;
	      i=0;
	      }
	    }
	  } while (finishflag==0);
	fprintf(stderr,"\b\b\b\b100%%\n\n");
	n1=0;
	n2=0;
	n3=0;
	if (nr_of_water_groups==1) {
	  fprintf(stderr,"The water molecules form one group.\n");
	  }
        else {
	  fprintf(stderr,"The water molecules form %d groups.\n",nr_of_water_groups);
	  for (i=1; i<=nr_of_water_groups; i++) {
	    count=0;
	    for (j=1; j<=nwater; j++) {
	      if (groupflags[j]==i) count++;
	      }
	    if (count==1) fprintf(stderr,"Group %ld contains 1 water molecule\n",i);
	    if (count!=1) fprintf(stderr,"Group %ld contains %ld water molecules\n",i,count);
	    if (count==1) n1++;
	    if (count==2) n2++;
	    if (count==3) n3++;
	    }
	  if (n1==1) fprintf(stderr,"So we have one single, ");
	  if (n1!=1) fprintf(stderr,"So we have %ld singles, ",n1);
	  if (n2==1) fprintf(stderr,"one twin, and ");
	  if (n2!=1) fprintf(stderr,"%ld twins, and ",n2);
	  if (n3==1) fprintf(stderr,"one triplet\n");
	  if (n3!=1) fprintf(stderr,"%ld triplets\n",n3);
	  }
	fprintf(stderr,"\n");
	if (write_surface_statistics_flag==1) {
	  fprintf(stderr,"Appending grouping statistics to file surface_stat.lis ... ");
	  if ((fp=fopen("surface_stat.lis","a"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not append to file surface_stat.lis\n");
	    terminate();
	    }
	  fprintf(fp,"[NUMBER OF GROUPS]\n");
	  fprintf(fp,"%d\n",nr_of_water_groups);
	  fprintf(fp,"[NUMBER OF SINGLES]\n");
	  fprintf(fp,"%ld\n",n1);
	  fprintf(fp,"[NUMBER OF TWINS]\n");
	  fprintf(fp,"%ld\n",n2);
	  fprintf(fp,"[NUMBER OF TRIPLETS]\n");
	  fprintf(fp,"%ld\n",n3);
	  fprintf(fp,"[GROUP OCCUPANCY]\n");
	  for (i=1; i<=nr_of_water_groups; i++) {
	    count=0;
	    for (j=1; j<=nwater; j++) {
	      if (groupflags[j]==i) count++;
	      }
	    fprintf(fp,"%ld %ld\n",i,count);
	    }
	  fclose(fp);
	  fprintf(stderr,"done.\n");
	  }
	return;
}       


static void place_ions()
{
	int number_cl, number_na;
	int number_ions;
	long int i0,i1,i2,i,j,j0;
	double xx,yy,zz,q,sum,dx,dy,dz,dd,r,delta_e,ex;
	double total_positive_charge,total_negative_charge;
    double total_positive_charge2,total_negative_charge2;
	double sum_na,sum_cl;
	double *xion;
	double *yion;
	double *zion;
	double *dens_na;
	double *dens_cl;
        double *prob_na;
        double *prob_cl;
        double *prob_na_total;
        double *prob_cl_total;
	double *protein_e_field;
	double *debye_hueckel_potential;
	int    *positions_na;
	int    *positions_cl;
        int    *ion_is_here;
	int    *positions_ion;
	int    *charge_ion;

	fprintf(stderr,"Placing ions ... ");
	replace_water_by_ion=ivector(1,nwater);
	nbulkwater=0;
	for (i=1; i<=nwater; i++) {
          if (groupflags[i]==1) nbulkwater++;
	  }
        xion=dvector(1,nbulkwater);
        yion=dvector(1,nbulkwater);
        zion=dvector(1,nbulkwater);
        debye_hueckel_potential=dvector(1,nbulkwater);
        protein_e_field=dvector(1,nbulkwater);
        dens_na=dvector(1,nbulkwater);
        dens_cl=dvector(1,nbulkwater);
        prob_na=dvector(1,nbulkwater);
        prob_cl=dvector(1,nbulkwater);
        prob_na_total=dvector(1,nbulkwater);
        prob_cl_total=dvector(1,nbulkwater);
	j=1;
	for (i=1; i<=nwater; i++) {
          if (groupflags[i]==1) {
	    xion[j] = xwater[i];
	    yion[j] = ywater[i];
	    zion[j] = zwater[i];
	    j++;
	    }
	  }
	for (j=1; j<=nbulkwater; j++) {
          protein_e_field[j] = 0.0;
          debye_hueckel_potential[j] = 0.0;
	  prob_na_total[j]=0.0;
	  prob_cl_total[j]=0.0;
	  }
	for (i=1; i<=natoms; i++) {
	  q=qatom[i];
	  if (q!=0.0) {
	    xx=xatom[i];
	    yy=yatom[i];
	    zz=zatom[i];
	    for (j=1; j<=nbulkwater; j++) {
	      dx=xion[j]-xx;
	      dy=yion[j]-yy;
	      dz=zion[j]-zz;
	      dd=sqrt(dx*dx+dy*dy+dz*dz);
	      protein_e_field[j] += q/dd;
	      debye_hueckel_potential[j] = (q*exp(-dd/DEBYE_HUECKEL_LENGTH)/dd);
	      }
	    sum_na = 0.0;
	    sum_cl = 0.0;
	    for (j=1; j<=nbulkwater; j++) {
	      prob_na[j]=-(ANGSTROM_PER_METER*ELEMENTARY_CHARGE*ELEMENTARY_CHARGE*debye_hueckel_potential[j]/KT/FOUR_PI_EPS_EPSR);
	      if (prob_na[j]<-1.0) prob_na[j]=-1.0;
	      prob_cl[j]= (ANGSTROM_PER_METER*ELEMENTARY_CHARGE*ELEMENTARY_CHARGE*debye_hueckel_potential[j]/KT/FOUR_PI_EPS_EPSR);
	      if (prob_cl[j]<-1.0) prob_cl[j]=-1.0;
	      sum_na += prob_na[j];
	      sum_cl += prob_cl[j];
	      }
	    for (j=1; j<=nbulkwater; j++) {
	      prob_na_total[j] -= (prob_na[j]*q/(2.0*sum_na));
	      prob_cl_total[j] += (prob_cl[j]*q/(2.0*sum_cl));
	      /* (why 1/2? because, analytically, the integral
	      over ion density change caused by one unit charge
	      is equal to 1/2 unit charges for sodium and also
	      1/2 (opposite sign) for chloride. that makes a total
	      of one unit charge, which shields the ion.) */
	      }
            }
          }
	for (j=1; j<=nbulkwater; j++) {
	  dens_na[j] = 1.0/ISOTONIC_WATER_PER_ION + prob_na_total[j];
	  if (dens_na[j]<0.0) dens_na[j]=0.0;
	  dens_cl[j] = 1.0/ISOTONIC_WATER_PER_ION + prob_cl_total[j];
	  if (dens_cl[j]<0.0) dens_cl[j]=0.0;
          }

	total_positive_charge=0.0;
	total_negative_charge=0.0;
	for (j=1; j<=nbulkwater; j++) {
	  total_positive_charge += dens_na[j];
	  total_negative_charge += dens_cl[j];
	  }
	number_na = floor(total_positive_charge+0.5);
	number_cl = floor(total_negative_charge+0.5);

	/* The user might have requested a total charge with the -q flag.
	 * If so, we adjust the number of ions to place here.
	 *
     * Do not modify total_positive_charge and total_negative_charge here,
     * since then the ion placement in the next code block cannot work properly */
	if (enforced_total_charge_flag==1)
	{
	    total_positive_charge2 = total_positive_charge;
	    total_negative_charge2 = total_negative_charge;
	    if (number_cl-number_na>floor(total_charge-enforced_total_charge+0.5)) {
	        do {
	            total_positive_charge2 += 0.001;
	            total_negative_charge2 -= 0.001;
	            number_na = floor(total_positive_charge2+0.5);
	            number_cl = floor(total_negative_charge2+0.5);
	        } while (number_cl-number_na>floor(total_charge-enforced_total_charge+0.5));
	    }
	    if (number_cl-number_na<floor(total_charge-enforced_total_charge+0.5)) {
	        do {
	            total_positive_charge2 -= 0.001;
	            total_negative_charge2 += 0.001;
	            number_na = floor(total_positive_charge2+0.5);
	            number_cl = floor(total_negative_charge2+0.5);
	        } while (number_cl-number_na<floor(total_charge-enforced_total_charge+0.5));
	    }
	    if (number_cl < 0 || number_na < 0)
	    {
	        fprintf(stderr, "\nERROR: Cannot match the requested total charge of %d.\n",
	                         enforced_total_charge);
	        exit(1);
	    }

	}

	if (number_na>0) positions_na = ivector(1,number_na);
	if (number_cl>0) positions_cl = ivector(1,number_cl);
        ion_is_here=ivector(1,nbulkwater);
	for (j=1; j<=nbulkwater; j++) {
	  ion_is_here[j]=0;
	  }

	/* place debye-hueckel chloride ions: */
	if (number_cl>0) {
	  for (i=1; i<=number_cl; i++) {
	    do {
              r=total_negative_charge*ran2(&seed);
	      sum=0.0;
	      for (j=1; sum<r; j++) {
	        i0=j;
	        if (dens_cl[j]>0.0) sum += dens_cl[j];
	        }
	      } while (ion_is_here[i0]==1);
            ion_is_here[i0]=1;
	    positions_cl[i]=i0;
	    }
	  }

	/* place debye-hueckel sodium ions: */
	if (number_na>0) {
	  for (i=1; i<=number_na; i++) {
	    do {
              r=total_positive_charge*ran2(&seed);
	      sum=0.0;
	      for (j=1; sum<r; j++) {
	        i0=j;
	        if (dens_na[j]>0.0) sum += dens_na[j];
	        }
	      } while (ion_is_here[i0]==1);
            ion_is_here[i0]=1;
	    positions_na[i]=i0;
	    }
	  }
        
        if ((number_na>0) && (number_cl>0)) {
	  fprintf(stderr,"placed %d sodium ion",number_na);
	  if (number_na!=1) fprintf(stderr,"s");
	  fprintf(stderr," and %d chloride ion",number_cl);
	  if (number_cl!=1) fprintf(stderr,"s");
	  fprintf(stderr,".\n");
	  }
        if ((number_na>0) && (number_cl==0)) {
	  fprintf(stderr,"placed %d sodium ion",number_na);
	  if (number_na!=1) fprintf(stderr,"s");
	  fprintf(stderr,".\n");
	  }
        if ((number_na==0) && (number_cl>0)) {
	  fprintf(stderr,"placed %d chloride ion",number_cl);
	  if (number_cl!=1) fprintf(stderr,"s");
	  fprintf(stderr,".\n");
	  }
        if ((number_na==0) && (number_cl==0)) {
	  fprintf(stderr,"placed no ions.\n");
	  }

	/* make a series of monte-carlo moves to equilibrate the ion cloud: */
	number_ions=0;
	if (number_na+number_cl>0) {
	  fprintf(stderr,"Equilibrating ions (%ld Monte Carlo moves) ...     ",(long int)NR_OF_DEBYE_HUECKEL_MONTE_CARLO_MOVES);
	  number_ions=number_na+number_cl;
	  positions_ion = ivector(1,number_ions);
	  charge_ion = ivector(1,number_ions);
	  if (number_na>0) {
	    for (i=1; i<=number_na; i++) {
	      positions_ion[i] = positions_na[i];
	      charge_ion[i]=1;
	      }
	    }
	  if (number_cl>0) {
	    for (i=1; i<=number_cl; i++) {
	      positions_ion[i+number_na] = positions_cl[i];
	      charge_ion[i+number_na]=-1;
	      }
	    }

  	  for (i=1; i<=NR_OF_DEBYE_HUECKEL_MONTE_CARLO_MOVES; i++) {
	    if (floor(((i-1.0)*100.0)/(1.0*NR_OF_DEBYE_HUECKEL_MONTE_CARLO_MOVES))!=floor((i*100.0)/(1.0*NR_OF_DEBYE_HUECKEL_MONTE_CARLO_MOVES))) {
	      fprintf(stderr,"\b\b\b\b%3d%%",(int)(floor((i*100.0)/(1.0*NR_OF_DEBYE_HUECKEL_MONTE_CARLO_MOVES))));
	      }
	    i0=1+floor(((float)number_ions)*ran2(&seed));
	    i1=positions_ion[i0];
	    do {
	      i2=1+floor(((float)nbulkwater)*ran2(&seed));
	      } while (ion_is_here[i2]==1);
	    delta_e = charge_ion[i0]*(protein_e_field[i2] - protein_e_field[i1]);
	    for (j=1; j<=number_ions; j++) {
	      j0=positions_ion[j];
	      if (j!=i0) {
		dx = xion[i2] - xion[j0];
		dy = yion[i2] - yion[j0];
		dz = zion[i2] - zion[j0];
		dd=sqrt(dx*dx+dy*dy+dz*dz);
		delta_e += (charge_ion[i0]*charge_ion[j]/dd);
		}
	      if (j0!=i1) {
		dx = xion[i1] - xion[j0];
		dy = yion[i1] - yion[j0];
		dz = zion[i1] - zion[j0];
		dd=sqrt(dx*dx+dy*dy+dz*dz);
		delta_e -= (charge_ion[i0]*charge_ion[j]/dd);
		}
              }
	    if (delta_e<0.0) {
	      ion_is_here[i1]=0;
	      positions_ion[i0]=i2;
	      ion_is_here[i2]=1;
	      }
	    else {
	      ex=(ANGSTROM_PER_METER*ELEMENTARY_CHARGE*ELEMENTARY_CHARGE/(KT*FOUR_PI_EPS_EPSR))*delta_e;
	      if (ex<100.0) {
	        if (exp(-ex)>ran2(&seed)) {
	          ion_is_here[i1]=0;
	          positions_ion[i0]=i2;
		  ion_is_here[i2]=1;
		  }
		}
	      }
            }
	  }
	fprintf(stderr,"\n\n");

	j=1;
	for (i=1; i<=nwater; i++) replace_water_by_ion[i]=0;
	if (number_na+number_cl>0) {
	  for (i=1; i<=nwater; i++) {
            if (groupflags[i]==1) {
	      if (ion_is_here[j]) {
	        for (i0=1; i0<=number_ions; i0++) {
	          if (positions_ion[i0]==j) replace_water_by_ion[i]=charge_ion[i0];
		  }
	        }
	      j++;
	      }
	    }
	  }
        free_dvector(xion,1,nbulkwater);
        free_dvector(yion,1,nbulkwater);
        free_dvector(zion,1,nbulkwater);
        free_dvector(dens_na,1,nbulkwater);
        free_dvector(dens_cl,1,nbulkwater);
        free_dvector(prob_na,1,nbulkwater);
        free_dvector(prob_cl,1,nbulkwater);
        free_dvector(prob_na_total,1,nbulkwater);
        free_dvector(prob_cl_total,1,nbulkwater);
        free_dvector(debye_hueckel_potential,1,nbulkwater);
        free_dvector(protein_e_field,1,nbulkwater);
        free_ivector(ion_is_here,1,nbulkwater);
	if (number_na+number_cl>0) {
	  if (number_na>0) free_ivector(positions_na,1,number_na);
	  if (number_cl>0) free_ivector(positions_cl,1,number_cl);
	  free_ivector(positions_ion,1,number_ions);
	  free_ivector(charge_ion,1,number_ions);
	  }
	nions=number_ions;
	return;
}





static void rotate_x(double *x1, double *y1, double *z1, double w)
{
	double dummy;

	dummy = (*y1)*cos(w) - (*z1)*sin(w);
	*z1 = (*y1)*sin(w) + (*z1)*cos(w);
	*y1 = dummy;
	return;
}



static void rotate_y(double *x1, double *y1, double *z1, double w)
{
	double dummy;

	dummy = (*z1)*cos(w) - (*x1)*sin(w);
	*x1 = (*z1)*sin(w) + (*x1)*cos(w);
	*z1 = dummy;
	return;
}



static void rotate_z(double *x1, double *y1, double *z1, double w)
{
	double dummy;

	dummy = (*x1)*cos(w) - (*y1)*sin(w);
	*y1 = (*x1)*sin(w) + (*y1)*cos(w);
	*x1 = dummy;
	return;
}





static void dice_hydrogens(double x0, double y0, double z0, double *x1, double *y1, double *z1, double *x2, double *y2, double *z2)
{
	double alpha,theta,phi;

	*x1 = BOND_LENGTH_H2O * sin(BOND_ANGLE_H2O/2.0);
	*x2 = -(*x1);
	*y1 = 0.0;
	*y2 = 0.0;
	*z1 = BOND_LENGTH_H2O * cos(BOND_ANGLE_H2O/2.0);
	*z2 = *z1;
	alpha = 2.0*PI*ran2(&seed);
	rotate_z(x1,y1,z1,alpha);
	rotate_z(x2,y2,z2,alpha);
	theta = acos(2.0*ran2(&seed)-1.0);
	rotate_y(x1,y1,z1,theta);
	rotate_y(x2,y2,z2,theta);
	phi = 2.0*PI*ran2(&seed);
	rotate_z(x1,y1,z1,phi);
	rotate_z(x2,y2,z2,phi);
	*x1 += x0;
	*y1 += y0;
	*z1 += z0;
	*x2 += x0;
	*y2 += y0;
	*z2 += z0;
	return;
}       



static void check_max_atom_nr(long int nr)
{
	if (nr>MAX_NR_OF_PDB_ATOMS-2) {
	  fprintf(stderr,"\n\nERROR: Too many atoms for pdb-file. Maximum is %ld.\n",(long int)MAX_NR_OF_PDB_ATOMS);
	  terminate();
	  }
	return;
}       


static void check_max_segment_nr(int nr)
{
	if (nr>MAX_NR_OF_PDB_SEGMENTS) {
	  fprintf(stderr,"\n\nERROR: Too many segments for pdb-file. Maximum is %d.\n",MAX_NR_OF_PDB_SEGMENTS);
	  fprintf(stderr,"       You seem to have encountered the unlikely situation\n");
	  fprintf(stderr,"       that there are more than 900 isolated water groups,\n");
	  fprintf(stderr,"       each of which is given an individual segment-id.\n");
	  fprintf(stderr,"       Look for lines 'res_id++;' in the source and remove them!\n");
	  terminate();
	  }
	return;
}       


static void write_xplor_segment(FILE *fp, int segid, long int nr, int diff)
{
	fprintf(fp,"segment\n");
	fprintf(fp,"  name=W%3d\n",segid);
	fprintf(fp,"  molecule\n");
	fprintf(fp,"    name=TIP3\n");
	fprintf(fp,"    number=%ld\n",nr);
	fprintf(fp,"  end\n");
	fprintf(fp,"end\n");
	fprintf(fp,"vector do (resid = encode(decode(resid) + %d))   (segid W%3d)\n\n",diff,segid);
	return;
}       

static void write_xplor_segment_na(FILE *fp, long int nr)
{
	fprintf(fp,"segment\n");
	fprintf(fp,"  name=NA\n");
	fprintf(fp,"  molecule\n");
	fprintf(fp,"    name=INA\n");
	fprintf(fp,"    number=%ld\n",nr);
	fprintf(fp,"  end\n");
	fprintf(fp,"end\n");
	return;
}       

static void write_xplor_segment_cl(FILE *fp, long int nr)
{
	fprintf(fp,"segment\n");
	fprintf(fp,"  name=CL\n");
	fprintf(fp,"  molecule\n");
	fprintf(fp,"    name=ICL\n");
	fprintf(fp,"    number=%ld\n",nr);
	fprintf(fp,"  end\n");
	fprintf(fp,"end\n");
	return;
}       


static void write_water_molecules(char *infile, char *outfile)
{
	long int i,count;
	int j,prev_res_id,res_id,seg_id;
	double x0,y0,z0,x1,y1,z1,x2,y2,z2;
	char   infile1[MAX_STRING_SIZE+1];
	char   outfile1[MAX_STRING_SIZE+1];
	char   xplorfile[MAX_STRING_SIZE+1];
	char   psfinfile[MAX_STRING_SIZE+1];
	char   psfoutfile[MAX_STRING_SIZE+1];
	FILE *fp_in, *fp_out, *fp_xplor;

	strcpy(infile1,infile);
	strcat(infile1,".pdb");
	strcpy(outfile1,outfile);
	strcat(outfile1,".pdb");
	strcpy(psfinfile,infile);
	strcat(psfinfile,".psf");
	strcpy(psfoutfile,outfile);
	strcat(psfoutfile,".psf");
	strcpy(xplorfile,XPLOR_SCRIPT_NAME);
	if (output_bulk_water_only_flag==0) {
	  if ((output_no_solute_flag==0) && (no_solute_present_flag==0)) {
	    if (add_ions_flag==0) {
	      fprintf(stderr,"Saving solute and %ld water molecules to file %s ... ",nwater,outfile1);
	      }
	    else {
	      fprintf(stderr,"Saving solute, %ld water molecules, and %ld ions to file %s ... ",nwater-nions,nions,outfile1);
	      }
	    }
	  else {
	    if (add_ions_flag==0) {
	      fprintf(stderr,"Saving %ld water molecules to file %s ... ",nwater,outfile1);
	      }
	    else {
	      fprintf(stderr,"Saving %ld water molecules and %ld ions to file %s ... ",nwater-nions,nions,outfile1);
	      }
	    }
	  }
	else {
	  count=0;
	  for (i=1; i<=nwater; i++) {if (groupflags[i]==1) count++;}
	  if ((output_no_solute_flag==0) && (no_solute_present_flag==0)) {
	    if (add_ions_flag==0) {
	      fprintf(stderr,"Saving solute and %ld bulk water molecules to file %s ... ",count,outfile1);
	      }
	    else {
	      fprintf(stderr,"Saving solute, %ld bulk water molecules, and %ld ions to file %s ... ",count-nions,nions,outfile1);
	      }
	    }
	  else {
	    if (add_ions_flag==0) {
	      fprintf(stderr,"Saving %ld bulk water molecules to file %s ... ",count,outfile1);
	      }
	    else {
	      fprintf(stderr,"Saving %ld bulk water molecules and %ld ions to file %s ... ",count-nions,nions,outfile1);
	      }
	    }
	  }
	if ((fp_out=fopen(outfile1,"w"))==NULL) {
	  fprintf(stderr,"\n\nERROR: Could not write to file %s .\n",outfile1);
	  terminate();
	  }
	fprintf(fp_out,"REMARK Created by SOLVATE 1.0.1 (C) 1996-2010 Helmut Grubmuller\n");
	fprintf(fp_out,"REMARK http://www.imo.physik.uni-muenchen.de/~grubi\n");
	fprintf(fp_out,"REMARK Helmut.Grubmueller@Physik.uni-muenchen.de\n");
	fprintf(fp_out,"REMARK \n");
	if ((output_no_solute_flag==0) && (no_solute_present_flag==0)) {
	  if ((fp_in=fopen(infile1,"r"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not read from file %s .\n",outfile1);
	    fprintf(stderr,"       (Yes, we still need the file as a template!)\n");
	    terminate();
	    }
	  for (i=1; i<=natoms; i++) {
            patch_pdb_coords(fp_in,fp_out,xatom[i],yatom[i],zatom[i]);
            }
          }
	else {
	  max_atom_nr_in_pdb=0;
	  }
	seg_id=100;
	res_id=1;
	if ((output_bulk_water_only_flag==0) && (nr_of_water_groups>1)) {
	  for (j=2; j<=nr_of_water_groups; j++) {
	    for (i=1; i<=nwater; i++) {
	      if (groupflags[i]==j) {
	        max_atom_nr_in_pdb++;
		check_max_atom_nr(max_atom_nr_in_pdb);
	        x0=xwater[i];
	        y0=ywater[i];
	        z0=zwater[i];
	        dice_hydrogens(x0,y0,z0,&x1,&y1,&z1,&x2,&y2,&z2);
                add_pdb_tip(fp_out,0,x0,y0,z0,max_atom_nr_in_pdb,res_id,seg_id);
	        max_atom_nr_in_pdb++;
                add_pdb_tip(fp_out,1,x1,y1,z1,max_atom_nr_in_pdb,res_id,seg_id);
	        max_atom_nr_in_pdb++;
                add_pdb_tip(fp_out,2,x2,y2,z2,max_atom_nr_in_pdb,res_id,seg_id);
		res_id++;
	        if (res_id>MAX_NR_OF_PDB_GROUPS) {
	          res_id=1;
	          seg_id++;
		  check_max_segment_nr(seg_id);
	          }
		}
	      }
	    if (res_id!=1) {
	      seg_id++;
	      check_max_segment_nr(seg_id);
	      }
	    }
	  }
	for (i=1; i<=nwater; i++) {
	  if ((groupflags[i]==1) && (replace_water_by_ion[i]==0)) {
	    max_atom_nr_in_pdb++;
	    check_max_atom_nr(max_atom_nr_in_pdb);
	    x0=xwater[i];
	    y0=ywater[i];
	    z0=zwater[i];
	    dice_hydrogens(x0,y0,z0,&x1,&y1,&z1,&x2,&y2,&z2);
            add_pdb_tip(fp_out,0,x0,y0,z0,max_atom_nr_in_pdb,res_id,seg_id);
	    max_atom_nr_in_pdb++;
            add_pdb_tip(fp_out,1,x1,y1,z1,max_atom_nr_in_pdb,res_id,seg_id);
	    max_atom_nr_in_pdb++;
            add_pdb_tip(fp_out,2,x2,y2,z2,max_atom_nr_in_pdb,res_id,seg_id);
	    res_id++;
	    if (res_id>MAX_NR_OF_PDB_GROUPS) {
	      res_id=1;
	      seg_id++;
	      check_max_segment_nr(seg_id);
	      }
	    }
	  }
	res_id=1;
	for (i=1; i<=nwater; i++) {
	  if (replace_water_by_ion[i]>0) {
	    max_atom_nr_in_pdb++;
	    check_max_atom_nr(max_atom_nr_in_pdb);
	    x0=xwater[i];
	    y0=ywater[i];
	    z0=zwater[i];
            add_pdb_ion(fp_out,replace_water_by_ion[i],x0,y0,z0,max_atom_nr_in_pdb,res_id);
	    res_id++;
	    if (res_id>MAX_NR_OF_PDB_GROUPS) {   /* very unlikely */
	      res_id=1;
	      }
	    }
	  }
	res_id=1;
	for (i=1; i<=nwater; i++) {
	  if (replace_water_by_ion[i]<0) {
	    max_atom_nr_in_pdb++;
	    check_max_atom_nr(max_atom_nr_in_pdb);
	    x0=xwater[i];
	    y0=ywater[i];
	    z0=zwater[i];
            add_pdb_ion(fp_out,replace_water_by_ion[i],x0,y0,z0,max_atom_nr_in_pdb,res_id);
	    res_id++;
	    if (res_id>MAX_NR_OF_PDB_GROUPS) {   /* very unlikely */
	      res_id=1;
	      }
	    }
	  }
	if ((output_no_solute_flag==0) && (no_solute_present_flag==0)) {
          fclose(fp_in);
	  }
	fprintf(fp_out,"END\n");
        fclose(fp_out);
	fprintf(stderr,"done.\n");

	/* output xplor-script */

	if (output_mkpsf_script==1) {
	  fprintf(stderr,"Writing XPLOR-script %s ... ",xplorfile);
	  if ((fp_xplor=fopen(xplorfile,"w"))==NULL) {
	    fprintf(stderr,"\n\nERROR: Could not write to file %s .\n",xplorfile);
	    terminate();
	    }
	  fprintf(fp_xplor,"!***************************************************************\n");
	  fprintf(fp_xplor,"!* Make structure file for water shell                         *\n");
	  fprintf(fp_xplor,"!* Created by SOLVATE 1.0.1 (C) 1996-2010 Helmut Grubmuller    *\n");
	  fprintf(fp_xplor,"!* http://www.mpibpc.mpg.de/home/grubmueller/downloads/solvate *\n");
	  fprintf(fp_xplor,"!* hgrubmu@gwdg.de                                             *\n");
	  fprintf(fp_xplor,"!***************************************************************\n\n\n");
	  fprintf(fp_xplor,"struct @%s end\n\n",psfinfile);
	  fprintf(fp_xplor,"topo @%s end\n",XPLOR_TOPOLOGY_FILE_WATER);
	  fprintf(fp_xplor,"topo @%s end\n",XPLOR_TOPOLOGY_FILE_IONS);
	  fprintf(fp_xplor,"param @%s end\n\n",XPLOR_PARAMETER_FILE_WATER);
	  fprintf(fp_xplor,"param @%s end\n\n",XPLOR_PARAMETER_FILE_IONS);
	  seg_id=100;
	  res_id=1;
	  count=0;
	  if ((output_bulk_water_only_flag==0) && (nr_of_water_groups>1)) {
	    for (j=2; j<=nr_of_water_groups; j++) {
	      prev_res_id=res_id;
	      for (i=1; i<=nwater; i++) {
	        if (groupflags[i]==j) {
		  count++;
		  res_id++;
	          if (res_id>MAX_NR_OF_PDB_GROUPS) {
		    write_xplor_segment(fp_xplor,seg_id,count,prev_res_id-1);
	            res_id=1;
	            prev_res_id=res_id;
	            seg_id++;
	            count=0;
	            }
		  }
	        }
	      if (res_id!=1) {
		write_xplor_segment(fp_xplor,seg_id,count,prev_res_id-1);
	        seg_id++;
	        count=0;
	        }
	      }
	    }
	  prev_res_id=res_id;
	  for (i=1; i<=nwater; i++) {
	    if ((groupflags[i]==1) && (replace_water_by_ion[i]==0)) {
	      count++;
	      res_id++;
	      if (res_id>MAX_NR_OF_PDB_GROUPS) {
		write_xplor_segment(fp_xplor,seg_id,count,prev_res_id-1);
	        res_id=1;
	        prev_res_id=res_id;
	        seg_id++;
		count=0;
	        }
	      }
	    }
	  write_xplor_segment(fp_xplor,seg_id,count,prev_res_id-1);
	  count=0;
	  for (i=1; i<=nwater; i++) {
	    if (replace_water_by_ion[i]==1) count++;
	    }
          if (count>0) {
	    write_xplor_segment_na(fp_xplor,count);
	    }
	  count=0;
	  for (i=1; i<=nwater; i++) {
	    if (replace_water_by_ion[i]==-1) count++;
	    }
          if (count>0) {
	    write_xplor_segment_cl(fp_xplor,count);
	    }

	  fprintf(fp_xplor,"coord @%s\n\n",outfile1);
	  fprintf(fp_xplor,"vector do (q = 0.0) (resname TIP3)\n");
	  fprintf(fp_xplor,"vector do (b = 0.0) (resname TIP3)\n\n");
	  fprintf(fp_xplor,"vector do (b = 0.0) (resname INA)\n\n");
	  fprintf(fp_xplor,"vector do (b = 0.0) (resname ICL)\n\n");
	  fprintf(fp_xplor,"write struct output=%s end\n",psfoutfile);
	  fprintf(fp_xplor,"write coors  output=%s end\n\n",outfile1);
	  fprintf(fp_xplor,"stop\n");
          fclose(fp_xplor);
	  fprintf(stderr,"done.\n");
	  fprintf(stderr,"(Use `xplor < %s' to create structure fiele %s.)\n",xplorfile,psfoutfile);
	  }
	free_ivector(groupflags,1,nwater);
	free_dvector(xwater,1,maxnwater);
	free_dvector(ywater,1,maxnwater);
	free_dvector(zwater,1,maxnwater);
	free_ivector(replace_water_by_ion,1,maxnwater);
	return;
}       






/*****************************************/
/*** came from numerical recipies ********/
/*****************************************/

static double	***dcube(int nrl, int nrh, int ncl, int nch, int npl, int nph)
{
	int i,j;
	double ***m;

	m=(double ***) alloc_ptrs(nrh-nrl+1);
	if (!m) {
	  fprintf(stderr,"\n\nERROR: allocation failure 1 in dmatrix()\n");
	  terminate();
	  }
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]= (double **) alloc_ptrs(nch-ncl+1);
		if (!m[i]) {
		  fprintf(stderr,"\n\nERROR: allocation failure 2 in dmatrix()\n");
		  terminate();
		  }
		m[i] -= ncl;
	}

	for(i=nrl;i<=nrh;i++) {
	  for(j=ncl;j<=nch;j++) {
		m[i][j]= alloc_doubles(nph-npl+1);
		if (!m[i][j]) {
		  fprintf(stderr,"\n\nERROR: allocation failure 3 in dmatrix()\n");
		  terminate();
		  }
		m[i][j] -= npl;
	  }
	}
	return m;
}



static void	free_dcube(double ***m, int nrl, int nrh, int ncl, int nch, int npl, int nph)
{
	int i=nph; /*** no warning/complete confusion? ***/
	int j=nch; /*** no warning/complete confusion? ***/

	for(i=nrh;i>=nrl;i--) {
	  for(j=nch;j>=ncl;j--) {
	    free((char*) (m[i][j]+npl));
	    }
	  }
	for(i=nrh;i>=nrl;i--) {
	  free((char*) (m[i]+ncl));
	  }
	free((char*) (m+nrl));
}


static double	**dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) alloc_ptrs(nrh-nrl+1);
	if (!m) {
	  fprintf(stderr,"\n\nERROR: allocation failure 1 in dmatrix()\n");
	  terminate();
	  }
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]= alloc_doubles(nch-ncl+1);
		if (!m[i]) {
		  fprintf(stderr,"\n\nERROR: allocation failure 2 in dmatrix()\n");
		  terminate();
		  }
		m[i] -= ncl;
	}
	return m;
}

static void	free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i=nch; /*** no warning ***/

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


/********************************************/


static int *ivector(int nl,int nh)
{
        int *v;

	v=alloc_ints(nh-nl+1);
	if (!v) {
	  fprintf(stderr,"\n\nERROR: allocation failure in ivector()\n");
	  terminate();
	  }
        return v-nl;
}

static void free_ivector(int *v, int nl, int nh)
{
	int i = nh; /*** no warning ***/
	i++;
	free((char*) (v+nl));
}

/********************************************/


static double *dvector(int nl,int nh)
{
        double *v;

	v=alloc_doubles(nh-nl+1);
	if (!v) {
	  fprintf(stderr,"\n\nERROR: allocation failure in dvector()\n");
	  terminate();
	  }
        return v-nl;
}

static void	free_dvector(double *v, int nl, int nh)
{
	int i = nh; /*** no warning ***/
	i++;
	free((char*) (v+nl));
}

/*********************************************/

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

static void	tred2(double **a,int n,double *d,double *e)
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = f>0 ? -sqrt(h) : sqrt(h);
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
			/* Next statement can be omitted if eigenvectors not wanted */
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	/* Next statement can be omitted if eigenvectors not wanted */
	d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}

/****************************************************************/
static void	tqli(double *d,double *e,int n,double **z)
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) {
				  fprintf(stderr,"\n\nERROR: Too many iterations in TQLI\n");
				  terminate();
				  }
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=sqrt((g*g)+1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					if (fabs(f) >= fabs(g)) {
						c=g/f;
						r=sqrt((c*c)+1.0);
						e[i+1]=f*r;
						c *= (s=1.0/r);
					} else {
						s=f/g;
						r=sqrt((s*s)+1.0);
						e[i+1]=g*r;
						s *= (c=1.0/r);
					}
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					p=s*r;
					d[i+1]=g+p;
					g=c*r-b;
				/* Next loop can be omitted if eigenvectors not wanted */
					for (k=1;k<=n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				d[l]=d[l]-p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

/***************************************************************************************/

static void sort(int n, double *ra)
{
        int l,j,ir,i;
        double rra;

        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1) {
                        rra=ra[--l];
                } else {
                        rra=ra[ir];
                        ra[ir]=ra[1];
                        if (--ir == 1) {
                                ra[1]=rra;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && ra[j] < ra[j+1]) ++j;
                        if (rra < ra[j]) {
                                ra[i]=ra[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                ra[i]=rra;
        }
}


static void sort2(int n, double *ra, double *rb)
{
        int l,j,ir,i;
        double rra;
        int rrb;

        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1) {
                        rra=ra[--l];
                        rrb=rb[l];
                } else {
                        rra=ra[ir];
                        rrb=rb[ir];
                        ra[ir]=ra[1];
                        rb[ir]=rb[1];
                        if (--ir == 1) {
                                ra[1]=rra;
                                rb[1]=rrb;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && ra[j] < ra[j+1]) ++j;
                        if (rra < ra[j]) {
                                ra[i]=ra[j];
                                rb[i]=rb[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                ra[i]=rra;
                rb[i]=rrb;
        }
}


static void sort5(int n, double *ra, double *rb, double *rc, double *rd, double *re)
{
	int l,j,ir,i;
	float rre,rrd,rrc,rrb,rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
			rrb=rb[l];
			rrc=rc[l];
			rrd=rd[l];
			rre=re[l];
		} else {
			rra=ra[ir];
			rrb=rb[ir];
			rrc=rc[ir];
			rrd=rd[ir];
			rre=re[ir];
			ra[ir]=ra[1];
			rb[ir]=rb[1];
			rc[ir]=rc[1];
			rd[ir]=rd[1];
			re[ir]=re[1];
			if (--ir == 1) {
				ra[1]=rra;
				rb[1]=rrb;
				rc[1]=rrc;
				rd[1]=rrd;
				re[1]=rre;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir)	{
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				rb[i]=rb[j];
				rc[i]=rc[j];
				rd[i]=rd[j];
				re[i]=re[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
		rb[i]=rrb;
		rc[i]=rrc;
		rd[i]=rrd;
		re[i]=rre;
	}
	return;
}




static int *alloc_ints(int len)
{
        int *ptr;
        size_t  nobytes;

        if (len<1) len = 1;
        nobytes = (size_t) (len * sizeof(int));
        ptr = (int *)malloc(nobytes);
        return (ptr);
}


static double *alloc_doubles(int len)
{
        double *ptr;
        size_t  nobytes;

        if (len<1) len = 1;
        nobytes = (size_t) (len * sizeof(double));
        ptr = (double *)malloc(nobytes);
        return (ptr);
}


static int   **alloc_ptrs(int len)
{
        int   **ptr;
        size_t  nobytes;

        if (len < 1)    len = 1;
        nobytes = (size_t) (len * sizeof(int*));
        ptr = (int **)malloc(nobytes);
        return (ptr);
}



#define M 714025
#define IA 1366
#define IC 150889

static double ran2(long *idum)
{
        static long iy,ir[98];
        static int iff=0;
        int j;

        if (*idum < 0 || iff == 0) {
                iff=1;
                if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
                for (j=1;j<=97;j++) {
                        *idum=(IA*(*idum)+IC) % M;
                        ir[j]=(*idum);
                }
                *idum=(IA*(*idum)+IC) % M;
                iy=(*idum);
        }
        j=1 + 97.0*iy/M;
        if (j > 97 || j < 1) {
	  fprintf(stderr,"\n\nERROR: This should really not happen. Sorry.\n");
	  terminate();
	  }
        iy=ir[j];
        *idum=(IA*(*idum)+IC) % M;
        ir[j]=(*idum);
        return (double) iy/M;
}

#undef M
#undef IA
#undef IC



static double rang(void)
/* Gibt Zufallsvariable mit in 8.Ordnung angenaeherter Gaussverteilung
   um 0.0 und Varianz 1.0 */
{
  int i;
  double r;

  r=0.0;
  for (i=1; i<=8; i++) {
    r=r+ran2(&seed);
    }
/*
        Zentrieren:
        v=v-8*1/2
        Mitteln:
        v=v/8
        1/(2*sqrt(3)) ist Varianz der Gleichverteilung in [0,1]:
        v=v*2*sqrt(3)
        Varianz der Summe von n Zufallszahlen geht mit 1/sqrt(n):
        v=v*sqrt(8)
        daher:
        v=(v-4)/sqrt(3/2)
*/
  return (r-4.0)*1.2247448713;
}




static void bad_input(void)
{
	  fprintf(stderr,"usage: solvate [-t thick] [-r radius] [-n ngauss] [-ug] [-ub]\n");
	  fprintf(stderr,"               [-s] [-v] [-bulk] [-w] [-ion] [-q charge] [-psf] [infile] outfile\n");
	  fprintf(stderr,"\n");
	  fprintf(stderr,"  thick:   minimum water shell thickness in Angstrom (default: 10.0)\n");
	  fprintf(stderr,"  radius:  maximum boundary curvature radius in Angstrom (default: 100000.0)\n");
	  fprintf(stderr,"  ngauss:  use ngauss gaussians to define solvent boundary (default: 1)\n");
	  fprintf(stderr,"  -s       write surface statistics to file surface_stat.lis\n");
	  fprintf(stderr,"  -v       write volume statistics to file volume_stat.lis\n");
	  fprintf(stderr,"  -ug      use gaussians from file gaussians.lis\n");
	  fprintf(stderr,"           (skip gaussian learning phase)\n");
	  fprintf(stderr,"  -ub      use boundary definition from file boundary.lis\n");
	  fprintf(stderr,"           (skip gaussian learning phase and boundary adjustment)\n");
	  fprintf(stderr,"  -bulk    output only bulk water (i.e., no buried water molecules)\n");
	  fprintf(stderr,"  -w       output only water molecules (i.e., not the solute)\n");
	  fprintf(stderr,"  -ion     add ions (na+/cl-) in isotonic/Debye-Hueckel concentration\n");
	  fprintf(stderr,"  charge:  force total charge of output system after ion placement\n");
	  fprintf(stderr,"           (Debye-Hueckel ion placement in small systems may lead to a net charge)\n");
	  fprintf(stderr,"  -psf     output XPLOR-script %s for generating a psf-file\n",XPLOR_SCRIPT_NAME);
	  fprintf(stderr,"  infile:  pdb-file (no extension) of solute\n");
	  fprintf(stderr,"           (if omitted, a water droplet is created)\n");
	  fprintf(stderr,"  outfile: pdb-file (no extension) of water shell\n");
	  fprintf(stderr,"\n");
          exit(1);
}



/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


int main(int argc,char *argv[])
{
	int i;
	double dummy;
        char infile[MAX_STRING_SIZE+1];
        char outfile[MAX_STRING_SIZE+1];

	fprintf(stderr,"\n\n\n");
	fprintf(stderr,"  .:::::..    SOLVATE 1.0.1 (C) 1996-2010 Helmut Grubmuller\n");
	fprintf(stderr," ::oOOOoo::        Purpose: make a water shell around solutes\n");
	fprintf(stderr," ::OOO::::        http://www.mpibpc.mpg.de/home/grubmueller/downloads/solvate\n");
	fprintf(stderr,"  :::::         hgrubmu@gwdg.de\n\n");
	if ((argc<2) ||(argc>16)) bad_input();

	water_shell_thickness = DEFAULT_WATER_SHELL_THICKNESS;
	boundary_radius = DEFAULT_BOUNDARY_RADIUS;
	ncodebook = DEFAULT_NR_OF_GAUSSIANS;
	write_surface_statistics_flag = 0;
	write_volume_statistics_flag = 0;
	use_gaussians_flag = 0;
	use_boundary_flag = 0;
	output_bulk_water_only_flag = 0;
	no_solute_present_flag = 0;
	output_no_solute_flag = 0;
	add_ions_flag = 0;
	output_mkpsf_script = 0;
	water_shell_thickness_flag = 0;
	ncodebook_flag = 0;
	sprintf(infile," ");
	sprintf(outfile," ");

	for (i=1; i<argc; i++) {
          if (strcmp(argv[i],"-t") == 0) {
	    if (i==argc-1) {
	      fprintf(stderr,"\n\nERROR: Missing water shell thickness.\n");
              exit(1);
	      }
	    water_shell_thickness=atof(argv[i+1]);
	    i++;
	    if (water_shell_thickness<1.0) {
	      fprintf(stderr,"\n\nERROR: Water shell thickness must be at least 1 A.\n");
              exit(1);
	      }
	    if (water_shell_thickness>1000.0) {
	      fprintf(stderr,"\n\nERROR: Water shell thickness must be < 1000 A.\n");
              exit(1);
	      }
	    water_shell_thickness_flag=1;
	    }
          else if (strcmp(argv[i],"-r") == 0) {
	    if (i==argc-1) {
	      fprintf(stderr,"\n\nERROR: Missing boundary curvature radius.\n");
              exit(1);
	      }
	    boundary_radius=atof(argv[i+1]);
	    i++;
	    if (boundary_radius<=0.0) {
	      fprintf(stderr,"\n\nERROR: A positive boundary curvature radius is recommended.\n");
              exit(1);
	      }
	    }
          else if (strcmp(argv[i],"-n") == 0) {
	    if (i==argc-1) {
	      fprintf(stderr,"\n\nERROR: Missing number of gaussians.\n");
              exit(1);
	      }
	    dummy=atof(argv[i+1]);
	    ncodebook=(int)(floor(dummy));
	    i++;
	    if (((double)ncodebook)!=dummy) {
	      fprintf(stderr,"\n\nERROR: Have you ever seen %f gaussians ?\n",dummy-ncodebook);
              exit(1);
	      }
	    if (ncodebook<1) {
	      fprintf(stderr,"\n\nERROR: You should at least use *one* gaussian.\n");
              exit(1);
	      }
	    if (ncodebook>1000) {
	      fprintf(stderr,"\n\nERROR: Use of more than 1000 gaussians is nonsense.\n");
              exit(1);
	      }
	    ncodebook_flag=1;
	    }
          else if (strcmp(argv[i],"-s") == 0) {
	    write_surface_statistics_flag = 1;
	    }
          else if (strcmp(argv[i],"-v") == 0) {
	    write_volume_statistics_flag = 1;
	    }
          else if (strcmp(argv[i],"-ug") == 0) {
	    use_gaussians_flag = 1;
	    }
          else if (strcmp(argv[i],"-ub") == 0) {
	    use_boundary_flag = 1;
	    }
          else if (strcmp(argv[i],"-bulk") == 0) {
	    output_bulk_water_only_flag = 1;
	    }
          else if (strcmp(argv[i],"-w") == 0) {
	    output_no_solute_flag = 1;
	    }
          else if (strcmp(argv[i],"-ion") == 0) {
	    add_ions_flag = 1;
	    }
          else if (strcmp(argv[i],"-q") == 0) {
              if (i==argc-1) {
                  fprintf(stderr,"\n\nERROR: Missing value for the total charge after ion placement.\n");
                      exit(1);
                  }
              enforced_total_charge=atof(argv[i+1]);
              i++;
              enforced_total_charge_flag=1;
          }
          else if (strcmp(argv[i],"-psf") == 0) {
	    output_mkpsf_script = 1;
	    }
	  else {
	    if (strncmp(argv[i],"-",1) == 0) {
	      fprintf(stderr,"\n\nERROR: Illegal option %s\n\n",argv[i]);
              bad_input();
	      }
	    else {
	      if (strcmp(infile," ") == 0) {
	        strcpy(infile,argv[i]);
		}
	      else {
	        if (strcmp(outfile," ") == 0) {
		  strcpy(outfile,argv[i]);
		  }
		else bad_input();
		}
	      }
	    }
	  }
	if (strcmp(infile," ") == 0) {
	  fprintf(stderr,"\n\nERROR: Missing output filename\n\n");
	  bad_input();
	  }
	if (strcmp(outfile," ") == 0) {
	  strcpy(outfile,infile);
	  strcpy(infile," ");
	  no_solute_present_flag=1;
	  }
	if ((use_boundary_flag == 1) && (water_shell_thickness_flag==1)) {
	  fprintf(stderr,"\n\nWARNING: Given water shell thickness is ignored,\n");
	  fprintf(stderr,"         since it is beeing read from file boundary.lis\n\n");
	  }
	if (((use_boundary_flag == 1) || (use_gaussians_flag==1)) && (ncodebook_flag==1)) {
	  fprintf(stderr,"\n\nWARNING: Given number of gaussians is ignored,\n");
	  fprintf(stderr,"         since it is beeing read from file.\n\n");
	  }
	if ((output_mkpsf_script == 1) && (output_no_solute_flag==1)) {
	  fprintf(stderr,"\n\nWARNING: You should not use both -w and -psf.\n");
	  fprintf(stderr,"         (XPLOR-script will not work if solute is not output.)\n\n");
	  }
	if (enforced_total_charge_flag==1 && add_ions_flag==0) {
	  fprintf(stderr,"\n\nWARNING: -q will only take effect together with -ion.\n\n");
	  }


	/*
	printf("\n");
	printf("water_shell_thickness:         %f\n",water_shell_thickness);
	printf("boundary_radius:               %f\n",boundary_radius);
	printf("ncodebook:                     %d\n",ncodebook);
	printf("write_surface_statistics_flag: %d\n",write_surface_statistics_flag);
	printf("write_volume_statistics_flag:  %d\n",write_volume_statistics_flag);
	printf("use_gaussians_flag:            %d\n",use_gaussians_flag);
	printf("use_boundary_flag:             %d\n",use_boundary_flag);
	printf("output_bulk_water_only_flag:   %d\n",output_bulk_water_only_flag);
	printf("no_solute_present_flag:        %d\n",no_solute_present_flag);
	printf("output_no_solute_flag:         %d\n",output_no_solute_flag);
	printf("output_mkpsf_script:           %d\n",output_mkpsf_script);
	printf("infile:                       |%s|\n",infile);
	printf("outfile:                      |%s|\n",outfile);
	printf("\n");
	*/
        seed=123;
	read_pdb_file(infile);
	if (add_ions_flag==1) read_atom_charges_from_psf_file(infile);
	if ((use_gaussians_flag==0) && (use_boundary_flag==0)) {
          calc_bounding_sphere();
          create_convex_volume();
	  dice_data();
	  learn_gaussians();
          save_gaussian_parameters();
	  adjust_surface_distance();
	  save_boundary_parameters();
	  }
	if ((use_gaussians_flag==1) && (use_boundary_flag==0)) {
	  read_gaussian_parameters();
	  adjust_surface_distance();
	  save_boundary_parameters();
	  }
	if ((use_gaussians_flag==0) && (use_boundary_flag==1)) {
	  read_boundary_parameters();
	  }
	if ((use_gaussians_flag==1) && (use_boundary_flag==1)) {
	  fprintf(stderr,"ERROR: You cannot use both  -ug and -ub\n\n");
	  exit(1);
	  }
	create_volume();
        surface_statistics();
        place_water_molecules();
        group_water_molecules();
	if (add_ions_flag==1) place_ions();
        write_water_molecules(infile,outfile);

	free_dmatrix(gauss_y,1,3,1,ncodebook);
	free_dcube(gauss_w,1,3,1,3,1,ncodebook);
	free_dmatrix(gauss_rho2,1,3,1,ncodebook);
	free_dvector(xatom,1,natoms);
	free_dvector(yatom,1,natoms);
	free_dvector(zatom,1,natoms);
	free_dvector(qatom,1,natoms);
	free_dvector(sigatom,1,natoms);
	free_dvector(epsatom,1,natoms);
	free_dvector(a_vdw_atom,1,natoms);
	free_dvector(b_vdw_atom,1,natoms);

	fprintf(stderr,"\nSolvate ready. No errors.\n\n");
	return 0;
}
