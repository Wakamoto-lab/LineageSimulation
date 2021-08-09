/************************
<additionalfunctions.h>
 ***********************/

/************* PRE-SET VARIABLES *****************/

/************************
    MERSENNE TWISTER
 ***********************/
/* NECESSARY FOR RANDOM VALUE GENERATION */
/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
//#define MAX 32768 //Not used in this header file

/************************
 GNUPLOT
 ***********************/
#define GNUPLOT_PATH "gnuplot -persist"
#define GNUPLOT_PATH_NoPersist "gnuplot"

/************* GLOBAL VARIABLES ***************/
/************************
    MERSENNE TWISTER
 ***********************/
static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/************* DECLARATION OF FUNCTIONS ***************/
/************************
    FILE NAME HANDLER
 ***********************/
void AddNumberExtToFileName(char fname[], int rank, char extention[],
			    long number);
void AddTwoNumberPlusExtToFileName(char fname[], int rank1, long number1, 
				   int rank2, long number2, char extention[]);
// The following functions are used in 'AddNumberExtToFileName()'
void NumberStringReturn(char numberstring[], int rank, long number);
void NumberOrder(long number, int *order);
int NumberOrderReturn(long number);

/************************
    FILE HANDLING
 ***********************/
void FileCopy(char fname_r[], char fname_w[]);
void FileContentShow(char fname[]);
long FileLineNumberReturn(char fname[]);

/************************
    MERSENNE TWISTER
 ***********************/
void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void); /* Random No in [0,1]*/
double genrand_real2(void); /* Random No in [0,1)*/
double genrand_real3(void); /* Random No in (0,1)*/
double genrand_res53(void);

/************************
RANDOM NUMBER GENERATORS
 ***********************/
double rand_range(double min, double max);
double gauss_rand(double mean, double SD);
double lognormal_rand(double meantake, double SDtake);
long binomial_rand(long trial, double prob);
long poisson(long average);
double gamma_rand_ave(double average , double k);
double gamma_rand(double k, double theta);
double NextValueGamma(double shape, double scale, double autocorr, double value);
// The following functions are used in 'gamma_rand()'
double exprand (void);
double gamrand (double a);

/************************
SPECIAL FUNCTIONS
 ************************/
double gammafunc(double x);
double loggamma(double x);
double igamma(double k, double x); //Incomplete gamma function
double p_gamma(double a, double x, double loggamma_a);
double q_gamma(double a, double x, double loggamma_a);
int comb(int n, int k);
long comb_long(long n, long k);

/************************
STATISTICS AND COLUMN HANDLING
 ***********************/
void MeanVarAdd(double *mean, double *var, double value);
void MeanVarFinal(double *mean, double *var, double counter);
void MeanVarCalcArray(double array[], long k,
		      double *mean, double *var);
double Correl(double mean_x, double mean_y, double SD_x, double SD_y, double xy_mean);
double CorrelFromFile(char fname[]);
double CorrelFromArray(double x[], double y[], long clNum);
double CorrelFromFile_ColumnSpecified(char fname[], int xcl, int ycl, int maxcl);
double CorrelFromFile_GrowthInverted(char fname[], int xindex, int yindex);
double CovarFromFile(char fname[]);
void ColumnAdd(double column[], long clnum, double value, double interval);
void ColumnInitialize_long(long column[], long clnum);
void ColumnInitialize_int(int column[], long clnum);
void ColumnInitialize_double(double column[], long clnum);
void ColumnRatio(double column[], long clnum, double counter);
void ColumnExport(double column[], long clnum, double interval, char fname[]);
void MakingColumnDataFromFile(char inputfname[], char outputfname[],
			      double interval);
void MakingColumnDataFromFile_ColumnSpecified(char inputfname[], 
					      char outputfname[],
					      double interval,
					      int target, int allNum);
void MakingColumnDataFromFile_Positive(char inputfname[], char outputfname[],
				       double interval);
void ExtractingNonZeroValueFromColumn(char fname[]);
void MeanVarCalcFile(char fname[], double *mean, double *var);
void MeanVarCalcFile_ColumnSpecified(char fname[], double *mean, double *var,
				     int target, int clnum);
void MeanVarCalcFile_GrowthConverted(char fname[], double *mean, double *var);
void RemoveLargeOutlierFromFile(char fname[], double maxValue);
void RemoveSmallOutlierFromFile(char fname[], double minValue);
void ExtractRegionDataFromFile(char fname[], double xmin, double xmax,
			       double ymin, double ymax);
void MakingCumulativeDataFile(char fname[], char fname_out[]);
void MakingSurvivalFuncFile(char fname[], char fname_out[]);
void MakingSurvivalFuncFromRawData(char fname[], char fname_out[],
				   double plotMin);
void MakingCumulativeDataFileFromRawData(char fname[], char fname_out[],
					 double plotMin);
void MakingCumulativeDataFile_wXInterval(char fname[], char fname_out[],
					 double factor);
void ColumnInitialize_double_increment(double column[], long clnum,
				       double initial_value, double increment);
void ColumnInitialize_long_increment(long column[], long clnum,
				     long initial_value, long increment);
void ColumnInitialize_int_increment(int column[], long clnum,
				    int initial_value, int increment);
double LeastSquares(double x[], double y[], long clnum, int index);
double LeastSquaresFromFile(char fname[], int index);
double LogLeastSquares(double x[], double y[], long clnum, int index);
double LeastSquaresYInterceptFixed(double x[], double y[], long clnum, double yintercept);
double LeastSquaresFromFileYInterceptFixed(char fname[], double yintercept);
void AcorrFromFile(char fname[], char fname_out[]);
void AcorrFromFile_Ver2(char fname[], char fname_out[], long plotIndex);
double MomentCalcFromFile(char fname[], int order);
double CumulantCalcFromFile(char fname[], int order);
double ArraySum_double(double array[], long k);
long ArraySum_long(long array[], long k);
int ArraySum_int(int array[], long k);
void GetWeibullParFile(char fname[], double *shape, double *scale);
void GetWeibullParameter(double value[], long k,
			 double *shape, double *scale);
void RandomDifferenceDataGenerator(char fname[], char fname_out[], long dataNum);
void FrechetParameterAssign(double *s, double *a,
			    double mean, double var);

/************************
 GNUPLOT
 ***********************/
void plot_linespoints(char filename[], char plottype[], char xlabel[], 
		      char ylabel[], int xlogscale, int ylogscale, 
		      int eps_output, int datanum);
void plot_linespoints_cl(char filename[], char plottype[], char xlabel[], 
			 char ylabel[], int xlogscale, int ylogscale, 
			 int eps_output, int clpos);
void plot_histogram(char filename[], char xlabel[], char ylabel[], 
		    float fill_level, int eps_output_index, int datanum);
double returnBoxWidth(char fname[], int dataNum);
void plot_histogram_width_specified(char filename[], char xlabel[], char ylabel[], 
				    float fill_level, int eps_output_index, int datanum,
				    double binSize, double binRatio);
void plot_histogram_wGaussFit(char filename[], char xlabel[], char ylabel[], 
			      float fill_level, int eps_output_index, 
			      int datanum,
			      double mean, double SD, double interval);
void plot_histogram_wGammaFit(char filename[], char xlabel[], char ylabel[], 
			      float fill_level, int eps_output_index, 
			      int datanum,
			      double mean, double SD, double interval);
void plot_histogram_wWeibullFit(char filename[], char xlabel[], char ylabel[], 
				float fill_level, int eps_output_index, 
				int datanum,
				double shape, double scale, double interval);
//gnuplot command interface
void gpclose(FILE *gp);
FILE *gpopen(void);
FILE *gpopen_nopersist(void);
void gpplot(FILE *gp, char fname[], char argument[]);
void gpcommand(FILE *gp, char argument[]);
void gpeps(FILE *gp, char fname[]);
void gpStandardFormat(FILE *gp);
void gpsave(FILE *gp, char fname[]);
void eps_output(FILE *gp, char fname[]);
void gpAxisLabels(FILE *gp, char xlabel[], char ylabel[]);
void ChangeExtentionTo(char fname[], char ext[]);
void gpLabels(FILE *gp, char xlabel[], char ylabel[]);

long FileLineNumberReturn(char fname[]);
void FileContentShow(char fname[]);
void FileCopy(char fname_r[], char fname_w[]);
int FileExistCheck(char fname[]);

void DataSort_LowToHigh(double value[], long k);
void DataSort_LowToHigh_DoubleArray(double value[], double value2[], long k);
void DataSort_LowToHigh_Long(long value[], long k);
void DataSort_LowToHigh_Int(int value[], long k);
void DataSort_LowToHigh_File(char fname[]);
void DataSort_HighToLow(double value[], long k);
void DataSort_HighToLow_Long(long value[], long k);
void DataSort_HighToLow_File(char fname[]);
double PercentileValueReturn(double c[], long k, double percentile);
double PercentileValueReturnFile(char fname[], double percentile);
double MaxValueReturnFile(char fname[]);
double MinValueReturnFile(char fname[]);

void InputPositiveDoubleValue(double *value, char sentence[]);
void InputNonNegativeDoubleValue(double *value, char sentence[]);
void InputPositiveLongValue(long *value, char sentence[]);
void InputNonNegativeLongValue(long *value, char sentence[]);
void InputPositiveIntValue(int *value, char sentence[]);
void InputNonNegativeIntValue(int *value, char sentence[]);

/************************
    STRING OPERATION
 ***********************/
void ExtractStringFromEnd(char result[], char original[], int sampleNo);
void latterStringReturn(char name[], char namesub[], int charnum);
int numberReturnFromString(char sentence[], int charnum);

/************* FUNCTION CONTENTS ***************/

/************************
    FILE NAME HANDLER
 ***********************/
void AddNumberExtToFileName(char fname[], int rank, char extention[],
			    long number)
{
  int i;
  int order=1;
  char fname_temp[100];
  char numberstring[100];

  NumberOrder(number, &order);//Determination of order
  if(rank<order){
    printf("Number is too large in 'AddNumberExtToFileName'\n");
    exit(EXIT_FAILURE);
  }

  NumberStringReturn(numberstring, rank, number);
  sprintf(fname_temp, "%s%s%s", fname, numberstring, extention);
  strcpy(fname, fname_temp);
}

void AddTwoNumberPlusExtToFileName(char fname[], int rank1, long number1, 
				   int rank2, long number2, char extention[])
{
  int i;
  int order=1;
  char fname_temp[100];
  char fname_temp2[100];
  char numberstring[100];

  NumberOrder(number1, &order);//Determination of order
  if(rank1<order){
    printf("Number is too large in 'AddTwoNumberPlusExtToFileName'\n");
    exit(EXIT_FAILURE);
  }
  NumberStringReturn(numberstring, rank1, number1);
  sprintf(fname_temp, "%s%s%s", fname, numberstring, "_");

  order = 1;
  NumberOrder(number2, &order);//Determination of order
  if(rank2<order){
    printf("Number is too large in 'AddTwoNumberPlusExtToFileName'\n");
    exit(EXIT_FAILURE);
  }
  NumberStringReturn(numberstring, rank2, number2);
  sprintf(fname_temp2, "%s%s%s", fname_temp, numberstring, extention);

  strcpy(fname, fname_temp2);
}

void NumberOrder(long number, int *order)
{
  if(number/10>=1){
    (*order)=(*order)+1;
    NumberOrder(number/10, order);
  }
}

int NumberOrderReturn(long number)
{
  int order = 1;
  NumberOrder(number, &order);
  return order;
}

void NumberStringReturn(char numberstring[], int rank, long number)
{
  int order=1;
  int i;
  char temp[100];

  NumberOrder(number, &order);//Determination of order
  if(rank<order){
    printf("Number is too large in 'NumberStringReturn'\n");
    exit(EXIT_FAILURE);
  }

  if(rank==order){
    sprintf(numberstring, "%ld", number);
  }else{
    for(i=0;i<rank-order;i++){
      temp[i]='0';
    }
    temp[rank-order]='\0';
    sprintf(numberstring, "%s%ld", temp, number);
  }
}

/************************
    MERSENNE TWISTER
 ***********************/
/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/***************************
  RANDOM NUMBER GENERATORS
 ***************************/
//This function returns a random number in the range (min, max)
double rand_range(double min, double max)
{

  if(min >= max){
    printf("ERROR: min is larger than max in rand_range()\n");
    exit(EXIT_FAILURE);
  }

  return min+(max-min)*genrand_real3();

}

//THIS FUNCTION RETURNS IN a and b GAUSS RANDOM NO WITH MEAN AND SD
double gauss_rand(double mean, double SD){

  double x, y;

  if(SD<0.0){
    printf("Minus SD value was passed in gauss_rand!!\n");
    exit(EXIT_FAILURE);
  }

  x=1.0-genrand_real1();
  y=1.0-genrand_real1();

  return pow(-2*log(x), 0.5)*sin(2*3.14159265*y)*SD+mean;
  
}

//THIS FUNCTION RETERNS LOG-NORMAL RANDOM NO WITH mean=meantake and sd=SDtake
double lognormal_rand(double meantake, double SDtake){

  double mean;
  double SD;

  if(SDtake<0.0 || meantake<0.0){
    printf("Minus mean or SD value was passed in lognormal_rand!!\n");
    exit(EXIT_FAILURE);
  }

  SD=pow(log(1+pow(SDtake,2.0)/pow(meantake,2.0)),0.5);
  mean=log(meantake)-log(1.0+(SDtake*SDtake)/(meantake*meantake))/2.0;

  return exp(gauss_rand(mean, SD));

}


//THIS FUNCTION RETURNS BINOMIAL RANDOM NO WITH N=trial and p=prob
long binomial_rand(long trial, double prob){

  long i, counter=0;

  if(prob<0 || prob>1){
    printf("Impossible probability in binomial!!\n");
    exit(EXIT_FAILURE);
  }

  if(trial>0){  
    for(i=0; i<trial; i++){
      if(genrand_real1()<=prob){
	counter++;
      }
    }
  }
	
  return counter;
}

long poisson(long average)
{
	double prob;
	long i, counter = 0, repeatnum;
	int threshold1 = 100;
	int threshold2 = 1000;
	
	if(average <= threshold1){
		repeatnum = 10000;
	}else if(average > threshold1 && average < threshold2){
		repeatnum = 100000;
	}else{
		repeatnum = 1000000;
	}
	
	prob = (double)average/(double)repeatnum;
	for(i=0; i<repeatnum; i++){
		if(genrand_real1()<=prob){
			counter++;
		}
	}	
	return counter;
}

//THIS FUNCTION RETURNS GAMMA RANDOM NUMBER WITH SHAPE
//PARAMETER k AND AVERAGE average
double gamma_rand_ave(double average , double k)
{
  double theta;

  theta=(double)average/k;

  return gamrand(k)*theta; //return gamma(k, theta);
}
//THIS FUNCTION RETURNS GAMMA RANDOM NUMBER WITH
//THE INPUT OF SHAPE AND SCALE PARAMETERS
double gamma_rand(double k, double theta)
{
  return gamrand(k)*theta;
}

double shifted_gamma_rand(double k, double theta, double shift)
{
  return gamma_rand(k, theta) + shift;
}

double gamrand (double a)
{
  double x, y, z;
  double u, v, w, b, c, e;
  double check;
  int accept = 0;
  if (a < 1)
    {
      /* Johnk's generator. Devroye (1986) p.418 */
      e = exprand();
      do {
	while(1){
	  check=genrand_real1();
	  if(check>0 && check<1){
	    x = pow(check, 1 / a);
	    break;
	  }
	}
	while(1){
	  check=genrand_real1();
	  if(check>0 && check<1){
	    y = pow(check, 1 / (1 - a));
	    break;
	  }
	}	
      } while (x + y > 1);
      return (e * x / (x + y));
    }else{
      /* Best's rejection algorithm. Devroye (1986) p.410 */
      b = a - 1;

      if(b<=0.0){
	return exprand();
      }else{
	c = 3 * a - 0.75;
	do {
	  /* generate */
	  while(1){
	    u = genrand_real1();
	    if(u>0 && u<1){
	      break;
	    }
	  }
	  while(1){
	    v = genrand_real1();
	    if(v>0 && v<1){
	      break;
	    }
	  }
	  w = u * (1 - u);
	  y = sqrt(c / w) * (u - 0.5);
	  x = b + y;
	  if (x >= 0)
	    {
	      z = 64 * w * w * w * v * v;
	      if (z <= 1 - (2 * y * y) / x)
		{
		  accept = 1;
		} else {
		  if (log(z) < 2 * (b * log(x / b) - y))
		    accept = 1;
		}
	    }
	} while (accept != 1);
	return x;
      }
    }
}



double NextValueGamma(double shape, double scale, double autocorr, double value)
{
  double alpha, beta;

  alpha = (value*autocorr/scale+shape*(1.0-autocorr))/(1.0-pow(autocorr, 2.0));
  beta = (1.0-pow(autocorr, 2.0))*scale;

  return gamma_rand(alpha, beta);

}


double exprand (void)
{
  double check;

  while(1){
    check=genrand_real1();
    if(check>0){
      break;
    }
  }
  return (- log(check));
}

/************************
SPECIAL FUNCTIONS
 ************************/
double gammafunc(double x)
{
  double PI = 3.14159265358979324;

    if (x < 0)
        return PI / (sin(PI * x) * exp(loggamma(1 - x)));
    return exp(loggamma(x));
}

//This function returns the logarithm of gamma function \Log(\Gamma(x))
double loggamma(double x)
{
  double v, w;
  double B0 =  1.0;
  double B1 = (-1.0 / 2.0);
  double B2 =  ( 1.0 / 6.0);
  double B4 =  (-1.0 / 30.0);
  double B6 =  ( 1.0 / 42.0);
  double B8 =  (-1.0 / 30.0);
  double B10 = ( 5.0 / 66.0);
  double B12 = (-691.0 / 2730.0);
  double B14 = ( 7.0 / 6.0);
  double B16 = (-3617.0 / 510.0);
  double threshold = 8.0;
  double LOG_2PI = 1.83787706640934548;

  v = 1;
  while (x < threshold) {  v *= x;  x++;  }
  w = 1 / (x * x);
  return ((((((((B16 / (16 * 15))  * w + (B14 / (14 * 13))) * w
	       + (B12 / (12 * 11))) * w + (B10 / (10 *  9))) * w
	     + (B8  / ( 8 *  7))) * w + (B6  / ( 6 *  5))) * w
	   + (B4  / ( 4 *  3))) * w + (B2  / ( 2 *  1))) / x
    + 0.5 * LOG_2PI - log(v) - x + (x - 0.5) * log(x);
}

//This function returns incomplete gamma function (\int_0^a x^(a-1)exp(-x)dx)/(\Gamma(a))
double igamma(double k, double x)
{
  return p_gamma(k, x, loggamma(k));
}

//This function returns incomplete gamma function (\int_0^a x^(a-1)exp(-x)dx)/(\Gamma(a))
//This function need an input of log(\Gamma(a))
double p_gamma(double a, double x, double loggamma_a)
{
  int k;
  double result, term, previous;
  
  if (x >= 1 + a) return 1 - q_gamma(a, x, loggamma_a);
  if (x == 0)     return 0;
  result = term = exp(a * log(x) - x - loggamma_a) / a;
  for (k = 1; k < 1000; k++) {
    term *= x / (a + k);
    previous = result;  result += term;
    if (result == previous) return result;
  }
  printf("p_gamma(): not converged!\n");
  return result;
}

//This function is 1-p_gamma()
double q_gamma(double a, double x, double loggamma_a)
{
  int k;
  double result, w, temp, previous;
  double la = 1, lb = 1 + x - a;
  
  if (x < 1 + a) return 1 - p_gamma(a, x, loggamma_a);
  w = exp(a * log(x) - x - loggamma_a);
  result = w / lb;
  for (k = 2; k < 1000; k++) {
    temp = ((k - 1 - a) * (lb - la) + (k + x) * lb) / k;
    la = lb;  lb = temp;
    w *= (k - 1 - a) / k;
    temp = w / (la * lb);
    previous = result;  result += temp;
    if (result == previous) return result;
  }
  printf("q_gamma(): Not converged!\n");
  return result;
}

int comb(int n, int k)
{
  if(k==0 || k == n) return 1;
  return comb(n-1, k-1)+comb(n-1, k);
}

long comb_long(long n, long k)
{
  if(k==0 || k == n) return 1;
  return comb(n-1, k-1)+comb(n-1, k);
}



/************************
STATISTICS AND COLUMN HANDLING
 ***********************/
void MeanVarAdd(double *mean, double *var, double value)
{
  (*mean)=(*mean)+value;
  (*var)=(*var)+pow(value, 2.0);
}

void MeanVarFinal(double *mean, double *var, double counter)
{
  (*mean)=(*mean)/counter;
  (*var)=(*var)/counter-pow((*mean), 2.0);
}

void MeanVarCalcFile(char fname[], double *mean, double *var)
{
  double value;
  double counter=0.0;
  FILE *fin=fopen(fname, "r");

  (*mean)=0.0;
  (*var)=0.0;

  while(fscanf(fin, "%lf", &value)!=EOF){
    MeanVarAdd(mean, var, value);
    counter=counter+1.0;
  }
  MeanVarFinal(mean, var, counter);
  fclose(fin);
}

void MeanVarCalcFile_ColumnSpecified(char fname[], double *mean, double *var,
				     int target, int clnum)
{
  double value;
  double counter=0.0;
  long lineNum;
  long i, j;
  FILE *fin=fopen(fname, "r");

  (*mean)=0.0;
  (*var)=0.0;

  lineNum = FileLineNumberReturn(fname);

  for(i=0; i<lineNum; i++){
    for(j=0; j<clnum; j++){
      fscanf(fin, "%lf", &value);
      if(j == target - 1){
	MeanVarAdd(mean, var, value);
	counter=counter+1.0;
      }
    }
  }

  MeanVarFinal(mean, var, counter);
  fclose(fin);
}

void MeanVarCalcFile_GrowthConverted(char fname[], double *mean, double *var)
{
  double value;
  double counter=0.0;
  FILE *fin=fopen(fname, "r");

  (*mean)=0.0;
  (*var)=0.0;

  while(fscanf(fin, "%lf", &value)!=EOF){
    if(value > 0.0){
      MeanVarAdd(mean, var, log(2.0)/value);
      counter=counter+1.0;
    }
  }
  MeanVarFinal(mean, var, counter);
  fclose(fin);
}

void ExtractRegionDataFromFile(char fname[], double xmin, double xmax,
			       double ymin, double ymax)
{
  FILE *fin, *fout;
  long dataNum;
  long i;
  double *datax, *datay;

  dataNum = FileLineNumberReturn(fname);
  datax = (double *)malloc(sizeof(double)*dataNum);
  datay = (double *)malloc(sizeof(double)*dataNum);
  if(datax == NULL || datay == NULL){
    printf("'datax' or 'datay' pointer couldn't be assigned!\n");
    exit(1);
  }

  fin = fopen(fname, "r");
  for(i=0; i<dataNum; i++){
    fscanf(fin, "%lf\t%lf\n", &datax[i], &datay[i]);
  }
  fclose(fin);

  fout = fopen(fname, "w");
  for(i=0; i<dataNum; i++){
    if(datax[i] >= xmin && datax[i] <= xmax &&
       datay[i] >= ymin && datay[i] <= ymax){
      fprintf(fout, "%lf\t%lf\n", datax[i], datay[i]);
    }
  }
  fclose(fout);

  printf("%ld lines were removed from %s\n",
	 dataNum-FileLineNumberReturn(fname), fname);
  
  free(datax);
  free(datay);
}

void RemoveLargeOutlierFromFile(char fname[], double maxValue)
{
  FILE *fin, *fout;
  long dataNum;
  long i;
  double *data;
  long counter = 0;

  dataNum = FileLineNumberReturn(fname);
  data = (double *)malloc(sizeof(double)*dataNum);
  if(data == NULL){
    printf("'data' pointer couldn't be assigned in 'RemoveLargeOutlierFromFile()'\n");
    exit(EXIT_FAILURE);
  }

  fin = fopen(fname, "r");
  for(i=0; i<dataNum; i++){
    fscanf(fin, "%lf", &data[i]);
  }
  fclose(fin);

  fout = fopen(fname, "w");
  for(i=0; i<dataNum; i++){
    if(data[i]<=maxValue){
      fprintf(fout, "%lf\n", data[i]);
    }else{
      printf("Outlier data %lf was removed from '%s'\n", data[i], fname);
      counter++;
    }
  }
  if(counter >0) printf("Total %ld data points were removed as outliers\n", counter);
  fclose(fout);
  free(data);
}


void RemoveSmallOutlierFromFile(char fname[], double minValue)
{
  FILE *fin, *fout;
  long dataNum;
  long i;
  double *data;
  long counter = 0;

  dataNum = FileLineNumberReturn(fname);
  data = (double *)malloc(sizeof(double)*dataNum);
  if(data == NULL){
    printf("'data' pointer couldn't be assigned in 'RemoveLargeOutlierFromFile()'\n");
    exit(EXIT_FAILURE);
  }

  fin = fopen(fname, "r");
  for(i=0; i<dataNum; i++){
    fscanf(fin, "%lf", &data[i]);
  }
  fclose(fin);

  fout = fopen(fname, "w");
  for(i=0; i<dataNum; i++){
    if(data[i]>=minValue){
      fprintf(fout, "%lf\n", data[i]);
    }else{
      printf("Outlier data %lf was removed from '%s'\n", data[i], fname);
      counter++;
    }
  }

  if(counter >0) printf("Total %ld data points were removed as outliers\n", counter);
  fclose(fout);
  free(data);
}

void MakingCumulativeDataFile(char fname[], char fname_out[])
{
  double x, y;
  double sum = 0.0;
  FILE *fin, *fout;

  fin = fopen(fname, "r");
  fout = fopen(fname_out, "w");

  while(fscanf(fin, "%lf", &x)!=EOF){
    fscanf(fin, "%lf", &y);
    sum = sum + y;
    fprintf(fout, "%lf\t%lf\n", x, sum);
  }
  fclose(fin);
  fclose(fout);
}


void MakingSurvivalFuncFile(char fname[], char fname_out[])
{
  double x, y;
  double sum = 0.0;
  FILE *fin, *fout;

  fin = fopen(fname, "r");
  fout = fopen(fname_out, "w");

  while(fscanf(fin, "%lf", &x)!=EOF){
    fscanf(fin, "%lf", &y);
    sum = sum + y;
    fprintf(fout, "%lf\t%lf\n", x, 1.0-sum);
  }
  fclose(fin);
  fclose(fout);
}

void MakingSurvivalFuncFromRawData(char fname[], char fname_out[],
				   double plotMin)
{
  long lineNum;
  double *value;
  long i;
  FILE *fin, *fout;

  lineNum = FileLineNumberReturn(fname);
  value = (double *)malloc(sizeof(double)*lineNum);
  if(value == NULL){
    printf("value pointer couldn't be assigned in MakingSurvivalFuncFromRawData()\n");
    exit(1);
  }

  fin = fopen(fname, "r");
  fout = fopen(fname_out, "w");
  fprintf(fout, "%lf\t1.0\n", plotMin);

  for(i=0; i<lineNum; i++){
    fscanf(fin, "%lf", &value[i]);
  }
  fclose(fin);

  DataSort_LowToHigh(value, lineNum);
  for(i=0; i<lineNum; i++){
    fprintf(fout, "%lf\t%lf\n", value[i], 1.0-((double)(i+1))/((double)lineNum));
  }
  fclose(fout);

  free(value);
}

void MakingCumulativeDataFileFromRawData(char fname[], char fname_out[],
					 double plotMin)
{
  long lineNum;
  double *value;
  long i;
  FILE *fin, *fout;

  lineNum = FileLineNumberReturn(fname);
  value = (double *)malloc(sizeof(double)*lineNum);
  if(value == NULL){
    printf("value pointer couldn't be assigned in MakingSurvivalFuncFromRawData()\n");
    exit(1);
  }

  fin = fopen(fname, "r");
  fout = fopen(fname_out, "w");
  fprintf(fout, "%lf\t0.0\n", plotMin);

  for(i=0; i<lineNum; i++){
    fscanf(fin, "%lf", &value[i]);
  }
  fclose(fin);

  DataSort_LowToHigh(value, lineNum);
  for(i=0; i<lineNum; i++){
    fprintf(fout, "%lf\t%lf\n", value[i], ((double)(i+1))/((double)lineNum));
  }
  fclose(fout);

  free(value);
}

void MakingCumulativeDataFile_wXInterval(char fname[], char fname_out[],
					 double factor)
{
  double x, y;
  double sum = 0.0;
  FILE *fin, *fout;

  fin = fopen(fname, "r");
  fout = fopen(fname_out, "w");

  while(fscanf(fin, "%lf", &x)!=EOF){
    fscanf(fin, "%lf", &y);
    sum = sum + y*factor;
    fprintf(fout, "%lf\t%lf\n", x, sum);
  }
  fclose(fin);
  fclose(fout);
}


void MeanVarCalcArray(double array[], long k,
		      double *mean, double *var)
{
  long i;
  double meanCalc = 0.0, varCalc = 0.0;

  for(i=0; i<k; i++){
    MeanVarAdd(&meanCalc, &varCalc, array[i]);
  }
  MeanVarFinal(&meanCalc, &varCalc, (double)k);
  (*mean) = meanCalc;
  (*var) = varCalc;
}

/******************************/
/* THIS FUNCTION RETURNS CORRELATION COEFFICIENT */
/******************************/
double Correl(double mean_x, double mean_y, double SD_x, double SD_y,
	      double xy_mean)
{
  return (xy_mean-mean_x*mean_y)/SD_x/SD_y;
}

double CorrelFromFile(char fname[])
{
  double counter=0.0;
  double x, y;
  double mean_x=0.0, mean_y=0.0, SD_x=0.0, SD_y=0.0, xy_mean=0.0;
  FILE *fin;

  fin = fopen(fname, "r");

  while(fscanf(fin, "%lf", &x)!=EOF){
    fscanf(fin, "%lf", &y);
    MeanVarAdd(&mean_x, &SD_x, x);
    MeanVarAdd(&mean_y, &SD_y, y);
    xy_mean = xy_mean + x*y;
    counter = counter + 1.0;
  }

  MeanVarFinal(&mean_x, &SD_x, counter);
  MeanVarFinal(&mean_y, &SD_y, counter);
  xy_mean = xy_mean/counter;

  fclose(fin);

  return Correl(mean_x, mean_y, pow(SD_x, 0.5), pow(SD_y, 0.5), xy_mean);
}

double CorrelFromArray(double x[], double y[], long clNum)
{
  double meanx, meany, varx, vary, meanxy;
  long i;

  MeanVarCalcArray(x, clNum, &meanx, &varx);
  MeanVarCalcArray(y, clNum, &meany, &vary);

  meanxy = 0.0;
  for(i=0; i<clNum; i++){
    meanxy = meanxy + x[i]*y[i];
  }
  meanxy = meanxy / clNum;

  return Correl(meanx, meany, pow(varx, 0.5), pow(vary, 0.5), meanxy);
  
}

double CorrelFromFile_ColumnSpecified(char fname[], int xcl, int ycl, int maxcl)
{
  double counter=0.0;
  double x, y;
  double mean_x=0.0, mean_y=0.0, SD_x=0.0, SD_y=0.0, xy_mean=0.0;
  int i;
  double *value;
  FILE *fin;

  value = (double *)malloc(sizeof(double)*maxcl);

  fin = fopen(fname, "r");

  while(fscanf(fin, "%lf", &x)!=EOF){
    value[0] = x;
    for(i=1; i<maxcl; i++){
      fscanf(fin, "%lf", &value[i]);
    }
    MeanVarAdd(&mean_x, &SD_x, value[xcl-1]);
    MeanVarAdd(&mean_y, &SD_y, value[ycl-1]);
    xy_mean = xy_mean + value[xcl-1]*value[ycl-1];
    counter = counter + 1.0;
  }

  MeanVarFinal(&mean_x, &SD_x, counter);
  MeanVarFinal(&mean_y, &SD_y, counter);
  xy_mean = xy_mean/counter;

  fclose(fin);
  free(value);

  return Correl(mean_x, mean_y, pow(SD_x, 0.5), pow(SD_y, 0.5), xy_mean);
}

//xindex = 1: log(2.0)/x; yindex = 1: log(2.0)/y
double CorrelFromFile_GrowthInverted(char fname[], int xindex, int yindex)
{
  double counter=0.0;
  double x, y;
  double mean_x=0.0, mean_y=0.0, SD_x=0.0, SD_y=0.0, xy_mean=0.0;
  FILE *fin;

  fin = fopen(fname, "r");

  while(fscanf(fin, "%lf", &x)!=EOF){
    fscanf(fin, "%lf", &y);
    if(x > 0.0 && y>0.0){
      if(xindex == 1) x = log(2.0)/x;
      if(yindex == 1) y = log(2.0)/y;
      MeanVarAdd(&mean_x, &SD_x, x);
      MeanVarAdd(&mean_y, &SD_y, y);
      xy_mean = xy_mean + x*y;
      counter = counter + 1.0;
    }
  }

  MeanVarFinal(&mean_x, &SD_x, counter);
  MeanVarFinal(&mean_y, &SD_y, counter);
  xy_mean = xy_mean/counter;

  fclose(fin);

  return Correl(mean_x, mean_y, pow(SD_x, 0.5), pow(SD_y, 0.5), xy_mean);
}

double CovarFromFile(char fname[])
{
  double counter=0.0;
  double x, y;
  double mean_x=0.0, mean_y=0.0, SD_x=0.0, SD_y=0.0, xy_mean=0.0;
  FILE *fin;

  fin = fopen(fname, "r");

  while(fscanf(fin, "%lf", &x)!=EOF){
    fscanf(fin, "%lf", &y);
    MeanVarAdd(&mean_x, &SD_x, x);
    MeanVarAdd(&mean_y, &SD_y, y);
    xy_mean = xy_mean + x*y;
    counter = counter + 1.0;
  }

  MeanVarFinal(&mean_x, &SD_x, counter);
  MeanVarFinal(&mean_y, &SD_y, counter);
  xy_mean = xy_mean/counter;

  fclose(fin);

  //return Correl(mean_x, mean_y, pow(SD_x, 0.5), pow(SD_y, 0.5), xy_mean);
  return xy_mean - mean_x*mean_y;
}

/******************************/
/* This function extract data from file and
   export column data*/
/******************************/
void MakingColumnDataFromFile(char inputfname[], char outputfname[],
			      double interval)
{
  double value;
  double column[1000];
  double counter=0.0;
  double max=0.0;
  long clnum;

  FILE* fin;

  fin=fopen(inputfname, "r");
  if(fin==NULL){
    printf("%s cannot be opened\n", inputfname);
    exit(EXIT_FAILURE);
  }
  ColumnInitialize_double(column, 1000);
  while(fscanf(fin, "%lf", &value)!=EOF){
    ColumnAdd(column, 1000, value, interval);
    if(max<value){
      max = value;
    }
    counter=counter+1.0;
  }
  fclose(fin);
  clnum = (long)(max/interval)+1;
  //added on 2010/04/26
  if(clnum == 1){
    clnum = 2;
  }
  ////
  ColumnRatio(column, clnum, counter);
  ColumnExport(column, clnum, interval, outputfname);
  ExtractingNonZeroValueFromColumn(outputfname);
}

void MakingColumnDataFromFile_ColumnSpecified(char inputfname[], 
					      char outputfname[],
					      double interval,
					      int target, int allNum)
{
  double value;
  double column[1000];
  double counter=0.0;
  double max=0.0;
  long clnum;
  long i,j;
  long lineNum;
  FILE* fin;

  if(FileExistCheck(inputfname) == 0){
    printf("File %s not exist!\n", inputfname);
    exit(1);
  }
  lineNum = FileLineNumberReturn(inputfname);
  fin=fopen(inputfname, "r");

  ColumnInitialize_double(column, 1000);
  for(i=0; i<lineNum; i++){
    for(j=0; j<allNum; j++){
      fscanf(fin, "%lf", &value);
      if(j==target-1){
	ColumnAdd(column, 1000, value, interval);
	if(max<value){
	  max = value;
	}
	counter=counter+1.0;
      }
    }
  }

  fclose(fin);
  clnum = (long)(max/interval)+1;
  //added on 2010/04/26
  if(clnum == 1){
    clnum = 2;
  }
  ////
  ColumnRatio(column, clnum, counter);
  ColumnExport(column, clnum, interval, outputfname);
  ExtractingNonZeroValueFromColumn(outputfname);
}

void MakingColumnDataFromFile_Positive(char inputfname[], char outputfname[],
				       double interval)
{
  double value;
  double column[1000];
  double counter=0.0;
  double max=0.0;
  long clnum;

  FILE* fin;

  fin=fopen(inputfname, "r");
  if(fin==NULL){
    printf("%s cannot be opened\n", inputfname);
    exit(EXIT_FAILURE);
  }
  ColumnInitialize_double(column, 1000);
  while(fscanf(fin, "%lf", &value)!=EOF){
    if(value>0.0){
      ColumnAdd(column, 1000, value, interval);
      if(max<value){
	max = value;
      }
      counter=counter+1.0;
    }
  }
  fclose(fin);
  clnum = (long)(max/interval)+1;
  //added on 2010/04/26
  if(clnum == 1){
    clnum = 2;
  }
  ////
  ColumnRatio(column, clnum, counter);
  ColumnExport(column, clnum, interval, outputfname);
  ExtractingNonZeroValueFromColumn(outputfname);
}

void ExtractingNonZeroValueFromColumn(char fname[])
{
  double x[10000], y[10000];
  long counter, i;
  FILE *fin, *fout;

  fin = fopen(fname, "r");
  counter = 0;
  while(fscanf(fin, "%lf", &x[counter])!=EOF){
    fscanf(fin, "%lf", &y[counter]);
    counter++;
  }
  fclose(fin);

  fout = fopen(fname, "w");
  for(i=0; i<counter; i++){
    if(y[i]>=0.000001){
      fprintf(fout, "%lf\t%lf\n", x[i], y[i]);
    }
  }
  fclose(fout);

}

/******************************/
/* THIS FUNCTION ADDS AN ENTRY
TO A COLUMN ACCORDING TO THE VALUE */
/******************************/
void ColumnAdd(double column[], long clnum, double value, double interval)
{
  long i;

  i=(long)(value/interval);
  if(i>=clnum){
    i=clnum-1;
  }
  if(i<0){
    i = 0;
  }

  column[i]=column[i]+1.0;
}

void ColumnInitialize_long(long column[], long clnum)
{
  long i;
  for(i=0;i<clnum;i++){
    column[i]=0;
  }
}

void ColumnInitialize_int(int column[], long clnum)
{
  long i;
  for(i=0;i<clnum;i++){
    column[i]=0;
  }
}

void ColumnInitialize_double(double column[], long clnum)
{
  long i;
  for(i=0;i<clnum;i++){
    column[i]=0.0;
  }
}

void ColumnRatio(double column[], long clnum, double counter)
{
  long i;
  for(i=0;i<clnum;i++){
    column[i]=column[i]/counter;
  }
}

void ColumnExport(double column[], long clnum, double interval, char fname[])
{
  long i;
  FILE* fout;

  fout=fopen(fname, "w");

  for(i=0;i<clnum;i++){
    fprintf(fout, "%lf\t%lf\n", ((double)i+0.5)*interval, column[i]);
  }

  fclose(fout);
}

void ColumnInitialize_double_increment(double column[], long clnum,
				       double initial_value, double increment)
{
  long i;
  for(i = 0; i < clnum; i++){
    column[i] = initial_value + increment * (double)i;
  }
}

void ColumnInitialize_long_increment(long column[], long clnum,
				     long initial_value, long increment)
{
  long i;
  for(i = 0; i < clnum; i++){
    column[i] = initial_value + increment * i;
  }
}

void ColumnInitialize_int_increment(int column[], long clnum,
				    int initial_value, int increment)
{
  long i;
  for(i = 0; i < clnum; i++){
    column[i] = initial_value + increment * i;
  }
}

double LeastSquaresFromFile(char fname[], int index)
{
  double *x, *y;
  long lineNum, i;
  double returnValue;
  FILE *fin;

  lineNum = FileLineNumberReturn(fname);
  x = (double *)malloc(sizeof(double)*lineNum);
  y = (double *)malloc(sizeof(double)*lineNum);
  if(x == NULL || y == NULL){
    printf("Pointer couldn't be assigned in 'LeastSquaresFromFile()'\n");
    exit(1);
  }

  fin = fopen(fname, "r");
  for(i=0; i<lineNum; i++){
    fscanf(fin, "%lf\t%lf", &x[i], &y[i]);
  }
  fclose(fin);

  returnValue = LeastSquares(x, y, lineNum, index);

  free(x);
  free(y);

  return returnValue;
}

//index = 0 returns slope, index = 1 returns interception
double LeastSquares(double x[], double y[], long clnum, int index)
{
  long i;
  double xsum = 0.0, ysum = 0.0, x2sum = 0.0, xysum = 0.0;
  double a, b;
  double retValue;

  if(clnum > 1){
    for(i=0; i<clnum; i++){
      xsum = xsum + x[i];
      ysum = ysum + y[i];
      x2sum = x2sum + (x[i]*x[i]);
      xysum = xysum +(x[i]*y[i]);
    }
    a = (clnum*xysum - xsum*ysum)/(clnum*x2sum-xsum*xsum);
    b = (x2sum*ysum-xysum*xsum)/(clnum*x2sum-xsum*xsum);
  }else{
    printf("Error in 'LeastSquares()': enough points not exist\n");
    exit(EXIT_FAILURE);
  }

  if(index == 0){
    retValue = a;
  }else{
    retValue = b;
  }

  return retValue;
}

double LeastSquaresYInterceptFixed(double x[], double y[], long clnum, double yintercept)
{
  double xysum = 0.0, xsum = 0.0, x2sum = 0.0;
  long i;

  if(clnum <= 0){
    printf("No data was found for least square fitting!\n");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<clnum; i++){
    xysum = xysum + x[i]*y[i];
    xsum = xsum + x[i];
    x2sum = x2sum + x[i]*x[i];
  }

  return (xysum-yintercept*xsum)/x2sum;
}

double LeastSquaresFromFileYInterceptFixed(char fname[], double yintercept)
{
  double *x, *y;
  long lineNum, i;
  double returnValue;
  FILE *fin;

  lineNum = FileLineNumberReturn(fname);
  x = (double *)malloc(sizeof(double)*lineNum);
  y = (double *)malloc(sizeof(double)*lineNum);
  if(x == NULL || y == NULL){
    printf("Pointer couldn't be assigned in 'LeastSquaresFromFile()'\n");
    exit(1);
  }

  fin = fopen(fname, "r");
  for(i=0; i<lineNum; i++){
    fscanf(fin, "%lf\t%lf", &x[i], &y[i]);
  }
  fclose(fin);

  returnValue = LeastSquaresYInterceptFixed(x, y, lineNum, yintercept);

  free(x);
  free(y);

  return returnValue;
}

//index = 0 returns slope, index = 1 returns interception
double LogLeastSquares(double x[], double y[], long clnum, int index)
{
  long i;
  double xvalue[10000];
  double logvalue[10000];
  long counter = 0;
  //int pos = 0;
  double retValue;

  for(i=0; i<clnum; i++){
    if(y[i]>0.0){
      xvalue[counter] = x[i];
      logvalue[counter] = log(y[i]);
      counter++;
    }
  }

  if(counter > 1){
    if(index == 0){
      retValue = LeastSquares(xvalue, logvalue, counter, index);
    }else{
      retValue = exp(LeastSquares(xvalue, logvalue, counter, index));
    }
  }else{
    //printf("Error in 'LogLeastSquares()': too few entries in array\n");
    //exit(EXIT_FAILURE);
    retValue = -10000.0; //assign unrealistic value
  }

  return retValue;
}

//THIS FUNCTION SHOWS AUTOCORRELATION FUNCTION CALCULATED
//FROM THE FILE
void AcorrFromFile(char fname[], char fname_out[])
{
  long i, j, counter;
  double time[10000], value[10000]; //maximum 10000 entries
  //double xmean[50000], ymean[50000], xymean[50000], xvar[50000], yvar[50000];
  double xmean, ymean, xymean, xvar, yvar;
  FILE *fin, *fout;

  fin = fopen(fname, "r");
  fout = fopen(fname_out, "w");
  printf("File open passed in 'AcorrFromFile'\n");

  counter = 0;
  while(fscanf(fin, "%lf", &time[counter])!=EOF){
    fscanf(fin, "%lf", &value[counter]);
    counter = counter+1;
    if(counter==10000){
      printf("Warning: Maximum 10000 points were taken in %s\n", fname);
      break;
    }
  }
  fclose(fin);
  printf("Getting data from %s done in 'AcorrFromFile'\n", fname);

  for(i=0; i<counter/2; i++){
    xmean = 0.0;
    ymean = 0.0;
    xymean = 0.0;
    xvar = 0.0;
    yvar = 0.0;
    for(j=0; j<counter-i; j++){
      xmean = xmean + value[j];
      ymean = ymean + value[j+i];
      xymean = xymean + value[j]*value[j+i];
      xvar = xvar + pow(value[j],2.0);
      yvar = yvar + pow(value[j+i],2.0);
    }
    xmean = xmean/(counter-i);
    ymean = ymean/(counter-i);
    xymean = xymean/(counter-i);
    xvar = xvar/(counter-i) - pow(xmean, 2.0);
    yvar = yvar/(counter-i) - pow(ymean, 2.0);
    fprintf(fout, "%lf\t%lf\n", time[i], Correl(xmean, ymean, pow(xvar, 0.5), pow(yvar, 0.5), xymean));
  }
  fclose(fout);

  plot_linespoints(fname_out, "p", "Time", "Autocorrelation", 0, 0, 0, 1);

}

void AcorrFromFile_Ver2(char fname[], char fname_out[], long plotIndex)
{
  long i, j, counter;
  double time[20000], value[20000]; //maximum 20000 entries
  //double xmean[50000], ymean[50000], xymean[50000], xvar[50000], yvar[50000];
  double xmean, ymean, xymean, xvar, yvar;
  FILE *fin, *fout;

  fin = fopen(fname, "r");
  fout = fopen(fname_out, "w");
  printf("File open passed in 'AcorrFromFile_Ver2'\n");

  counter = 0;
  while(fscanf(fin, "%lf", &value[counter])!=EOF){
    //fscanf(fin, "%lf", &value[counter]);
    counter = counter+1;
    if(counter==20000){
      printf("Warning: Maximum 20000 points were taken in %s\n", fname);
      break;
    }
  }
  fclose(fin);
  printf("Getting data from %s done in 'AcorrFromFile_ver2'\n", fname);

  for(i=0; i<counter/2; i++){
    xmean = 0.0;
    ymean = 0.0;
    xymean = 0.0;
    xvar = 0.0;
    yvar = 0.0;
    for(j=0; j<counter-i; j++){
      xmean = xmean + value[j];
      ymean = ymean + value[j+i];
      xymean = xymean + value[j]*value[j+i];
      xvar = xvar + pow(value[j],2.0);
      yvar = yvar + pow(value[j+i],2.0);
    }
    xmean = xmean/(counter-i);
    ymean = ymean/(counter-i);
    xymean = xymean/(counter-i);
    xvar = xvar/(counter-i) - pow(xmean, 2.0);
    yvar = yvar/(counter-i) - pow(ymean, 2.0);
    fprintf(fout, "%ld\t%lf\n", i, Correl(xmean, ymean, pow(xvar, 0.5), pow(yvar, 0.5), xymean));
  }
  fclose(fout);

  if(plotIndex == 1){
    plot_linespoints(fname_out, "p", "Time", "Autocorrelation", 0, 0, 0, 1);
  }

}

//This function returns the n-th order moment
//from the file 'fname[]'.
double MomentCalcFromFile(char fname[], int order)
{
  double value;
  double moment;
  double count;
  FILE *fin;

  fin = fopen(fname, "r");

  moment = 0.0;
  count = 0.0;
  while(fscanf(fin, "%lf", &value)!=EOF){
    moment = moment + pow(value, (double)order);
    count = count + 1.0;
  }

  if(count > 0.0){
    moment = moment/count;
  }else{
    printf("No data point found on the file '%s'!\n", fname);
    printf("The value 0.0 was returned\n");
    moment = 0.0;
  }

  fclose(fin);

  return moment;
}


double CumulantCalcFromFile(char fname[], int order)
{
  int i;
  double moment[4];
  //double cumulant[4];
  double returnValue;

  if(order > 4 || order <= 0){
    printf("Wrong 'order' was passed to 'CumulantCalcFromFile()'\n");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<order; i++){
    moment[i] = MomentCalcFromFile(fname, i+1);
  }

  switch(order){
  case 1:
    returnValue = moment[0];
    break;
  case 2:
    returnValue = moment[1] - pow(moment[0], 2.0);
    break;
  case 3:
    returnValue = (moment[2] - 3.0*moment[0]*moment[1] + 2.0*pow(moment[0],3.0))/pow(moment[1]-pow(moment[0], 2.0), 1.5);
    break;
  case 4:
    returnValue = (moment[3] - 4.0*moment[0]*moment[2] + 6.0*pow(moment[0],2.0)*moment[1]-3.0*pow(moment[0], 4.0))/pow(moment[1]-pow(moment[0],2.0), 2.0);
    break;
  default:
    break;
  }

  return returnValue;

}

double ArraySum_double(double array[], long k)
{
  long i;
  double sum = 0.0;
  for(i=0; i<k; i++){
    sum = sum + array[i];
  }
  return sum;
}

long ArraySum_long(long array[], long k)
{
  long i;
  long sum = 0;
  for(i=0; i<k; i++){
    sum = sum + array[i];
  }
  return sum;
}

int ArraySum_int(int array[], long k)
{
  long i;
  int sum = 0;
  for(i=0; i<k; i++){
    sum = sum + array[i];
  }
  return sum;
}


void GetWeibullParFile(char fname[], double *shape, double *scale)
{
  long counter;
  double value[50000];
  FILE *fin;

  fin = fopen(fname, "r");
  if(fin == NULL){
    printf("ERROR: File '%s' not found\n", fname);
    exit(EXIT_FAILURE);
  }

  counter = 0;
  while(fscanf(fin, "%lf", &value[counter])!=EOF && counter < 50000){
    counter++;
  }

  GetWeibullParameter(value, counter, shape, scale);
}

void FrechetParameterAssign(double *s, double *a,
			    double mean, double var)
{
  double dev = 20.0;
  double alpha[11];
  double checkValue[11];
  double interval = 2.0;
  int i;

  //Assigning the candidate of alpha
  alpha[0] = 3.0;
  for(i=0; i<11; i++){
    alpha[i] = alpha[0] + i*interval;
  }

  while(dev > 0.000001){
    for(i=0; i<11; i++){
      checkValue[i] = 
	tgamma(1.0 - 2.0/alpha[i])/pow(tgamma(1.0 - 1.0/alpha[i]), 2.0) - 1.0 - var/mean/mean;
      //printf("%lf\n", checkValue[i]);
    }
    i = 0;
    while(1){
      if(checkValue[i]*checkValue[i+1] > 0.0){
	i++;
      }else{
	break;
      }
      if(i==10){
	printf("Error: Root not found!\n");
	exit(1);
      }
    }
    dev = alpha[i+1] - alpha[i];
    interval = dev/10.0;
    alpha[0] = alpha[i];
    for(i=0; i<11; i++){
      alpha[i] = alpha[0] + interval*i;
    }
  }

  (*a) = (alpha[10]+alpha[0])/2.0;
  (*s) = mean/tgamma(1.0-1.0/(*a));
}


void GetWeibullParameter(double value[], long k,
			 double *shape, double *scale)
{
  double tempValue[50000];
  double cumulativeProb[50000];
  double slope, intercept;
  long i;
  char fname[256];
  FILE *fout, *gp;

  strcpy(fname, "WeibullFitting.dat");
  fout = fopen(fname, "w");

  if(k > 50000){
    printf("Too many values for sorting in 'GetWeibullParameter()'\n");
    printf("First 50000 entries will be taken\n");
    for(i=0; i<50000; i++){
      tempValue[i] = log(value[i]);
      cumulativeProb[i] = (double)(i+1)/50000.0;
      cumulativeProb[i] = log(-log(1.0-cumulativeProb[i])); //For fitting
    }
    DataSort_LowToHigh(tempValue, 50000);
    printf("Data sort done in GetWeibullParameter()\n");
    slope = LeastSquares(tempValue, cumulativeProb, 49999, 0);
    intercept = LeastSquares(tempValue, cumulativeProb, 49999, 1);
    for(i=0; i<50000; i++){
      fprintf(fout, "%lf\t%lf\n", tempValue[i], cumulativeProb[i]);
    }
  }else{
    for(i=0; i<k; i++){
      tempValue[i] = log(value[i]);
      cumulativeProb[i] = (double)(i+1)/(double)k;
      cumulativeProb[i] = log(-log(1.0-cumulativeProb[i])); //For fitting
    }
    DataSort_LowToHigh(tempValue, k);
    printf("Data sort done in GetWeibullParameter()\n");
    slope = LeastSquares(tempValue, cumulativeProb, k-1, 0);
    intercept = LeastSquares(tempValue, cumulativeProb, k-1, 1);
    for(i=0; i<k; i++){
      fprintf(fout, "%lf\t%lf\n", tempValue[i], cumulativeProb[i]);
    }
  }
  fclose(fout);

  (*shape) = slope;
  (*scale) = exp(-intercept/slope);

  printf("\n<Weibull dist fitting>\n");
  printf("Shape: %lf\n", (*shape));
  printf("Scale: %lf\n", (*scale));

  //Just for check
  gp = gpopen();
  gpStandardFormat(gp);
  gpLabels(gp, "log(x)", "log(-log(1-F[x]))");
  fprintf(gp, "f(x) = %lf*x + %lf\n", slope, intercept);
  fprintf(gp, "plot '%s' w p ps 2 ti 'Data'\n", fname);
  fprintf(gp, "replot f(x) w l lw 2 ti 'Fitting'\n");
  gpclose(gp);
}


void RandomDifferenceDataGenerator(char fname[], char fname_out[], long dataNum)
{
  long i, counter;
  double originalData[20000];
  long pos1, pos2;
  FILE *fin, *fout;

  fin = fopen(fname, "r");
  i=0;
  while(fscanf(fin, "%lf", &originalData[i])!=EOF && i<20000){
    i++;
    if(i==20000){
      printf("Maximum data number, 20000, reached in 'RandomDifferenceDataGenerator()'\n");
    }
  }
  counter = i;
  fclose(fin);

  fout = fopen(fname_out, "w");
  for(i=0; i<dataNum; i++){
    pos1 = (long)(genrand_real3()*counter);
    pos2 = (long)(genrand_real3()*counter);
    fprintf(fout, "%lf\n", fabs(originalData[pos1]-originalData[pos2]));
  }
  fclose(fout);
}



/************************
 GNUPLOT
 ***********************/
void plot_linespoints(char filename[], char plottype[], char xlabel[], 
		      char ylabel[], int xlogscale, int ylogscale, 
		      int eps_output_index, int datanum)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }

  gpStandardFormat(gp);
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);

  if(xlogscale==1){
    fprintf(gp, "set logscale x\n");
  }
  if(ylogscale==1){
    fprintf(gp, "set logscale y\n");
  }

  if(strcmp(plottype, "l") == 0){
    for(i=0;i<datanum;i++){
      if(i==0){
	fprintf(gp, "plot '%s' using 1:2 w %s lw 2 notitle\n", 
		filename, plottype);
      }else{
	fprintf(gp, "replot '%s' using 1:%d w %s lw 2 notitle\n",
		filename, i+2, plottype);
      }
    }
  }else if(strcmp(plottype, "p")==0){
    for(i=0;i<datanum;i++){
      if(i==0){
	fprintf(gp, "plot '%s' using 1:2 w %s ps 2 notitle\n", 
		filename, plottype);
      }else{
	fprintf(gp, "replot '%s' using 1:%d w %s ps 2 notitle\n",
		filename, i+2, plottype);
      }
    }
  }else{
    for(i=0;i<datanum;i++){
      if(i==0){
	fprintf(gp, "plot '%s' using 1:2 w %s lw 2 ps 2 notitle\n", 
		filename, plottype);
      }else{
	fprintf(gp, "replot '%s' using 1:%d w %s lw 2 ps 2 notitle\n",
		filename, i+2, plottype);
      }
    }
  }

  ChangeExtentionTo(fname_plt, ".plt");
  gpsave(gp, fname_plt);
  //fprintf(gp, "save '%s'\n", fname_plt);

  if(eps_output_index == 1){
    eps_output(gp, fname_plt);
  }

  fflush(gp);
  getchar();
  pclose(gp);
}

void plot_linespoints_cl(char filename[], char plottype[], char xlabel[], 
			 char ylabel[], int xlogscale, int ylogscale, 
			 int eps_output, int clpos)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }

  gpStandardFormat(gp);
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);

  if(xlogscale==1){
    fprintf(gp, "set logscale x\n");
  }
  if(ylogscale==1){
    fprintf(gp, "set logscale y\n");
  }

  fprintf(gp, "plot '%s' using 1:%d w %s lw 2 ps 2 notitle\n", 
	  filename, clpos, plottype);

  ChangeExtentionTo(fname_plt, ".plt");
  fprintf(gp, "save '%s'\n", fname_plt);

  if(eps_output==1){
    fprintf(gp, "set terminal postscript eps enhanced color\n");
    fprintf(gp, "set output '%s.eps'\n", filename);
    fprintf(gp, "replot\n");
    fprintf(gp, "set terminal x11\n");
    fprintf(gp, "replot\n");
  }

  fflush(gp);
  getchar();
  pclose(gp);
}

/*****************************/
/*  Histogram plot  */
/*****************************/
void plot_histogram(char filename[], char xlabel[], char ylabel[], 
		    float fill_level, int eps_output_index, int datanum)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }
  gpStandardFormat(gp);
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);

  for(i=0;i<datanum;i++){
    if(i==0){
      fprintf(gp, "plot '%s' using 1:2 w boxes fs solid %f notitle\n", 
	      filename, fill_level);
    }else{
      fprintf(gp, "replot '%s' using 1:%d w boxes fs solid %f notitle\n",
	      filename, i+2, fill_level);
    }
  }

  ChangeExtentionTo(fname_plt, ".plt");
  fprintf(gp, "save '%s'\n", fname_plt);

  if(eps_output_index==1){
    eps_output(gp, fname_plt);
  }

  fflush(gp);
  getchar();
  pclose(gp);
}

double returnBoxWidth(char fname[], int dataNum)
{
  double x1, x2, temp;
  int i;
  FILE *fin;

  fin = fopen(fname, "r");
  fscanf(fin, "%lf", &x1);
  for(i=0; i<dataNum; i++) fscanf(fin, "%lf", &temp);
  fscanf(fin, "%lf", &x2);
  fclose(fin);

  return x2-x1; 
}

void plot_histogram_width_specified(char filename[], char xlabel[], char ylabel[], 
				    float fill_level, int eps_output_index, int datanum,
				    double binSize, double binRatio)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }

  gpStandardFormat(gp);
  fprintf(gp, "set xtics %lf\n", binSize*5);
  fprintf(gp, "set mxtics 5\n");

  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);
  fprintf(gp, "set boxwidth %lf\n", binSize*binRatio); //Added

  for(i=0;i<datanum;i++){
    if(i==0){
      fprintf(gp, "plot '%s' using 1:2 w boxes fs solid %f notitle\n", 
	      filename, fill_level);
    }else{
      fprintf(gp, "replot '%s' using 1:%d w boxes fs solid %f notitle\n",
	      filename, i+2, fill_level);
    }
  }

  //ChangeExtentionTo(fname_plt, ".plt");
  //fprintf(gp, "save '%s'\n", fname_plt);
  gpsave(gp, fname_plt);

  if(eps_output_index==1){
    eps_output(gp, fname_plt);
  }

  fflush(gp);
  getchar();
  pclose(gp);
}

void plot_histogram_wGaussFit(char filename[], char xlabel[], char ylabel[], 
			      float fill_level, int eps_output_index, 
			      int datanum,
			      double mean, double SD, double interval)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }

  gpStandardFormat(gp);
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);

  for(i=0;i<datanum;i++){
    if(i==0){
      fprintf(gp, "plot '%s' using 1:2 w boxes fs solid %f notitle\n", 
	      filename, fill_level);
    }else{
      fprintf(gp, "replot '%s' using 1:%d w boxes fs solid %f notitle\n",
	      filename, i+2, fill_level);
    }
  }

  fprintf(gp, "m = %lf\n", mean);
  fprintf(gp, "s = %lf\n", SD);
  fprintf(gp, "f(x)=%lf/((2.0*3.14159)**0.5)/s*exp(-((x-m)**2.0)/2.0/(s**2.0))\n", interval);
  fprintf(gp, "replot f(x) w l lw 2 notitle\n");

  ChangeExtentionTo(fname_plt, ".plt");
  fprintf(gp, "save '%s'\n", fname_plt);

  if(eps_output_index == 1){
    eps_output(gp, fname_plt);
  }

  fflush(gp);
  getchar();
  pclose(gp);
}

void plot_histogram_wGammaFit(char filename[], char xlabel[], char ylabel[], 
			      float fill_level, int eps_output_index, 
			      int datanum,
			      double mean, double SD, double interval)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }

  gpStandardFormat(gp);
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);

  for(i=0;i<datanum;i++){
    if(i==0){
      fprintf(gp, "plot '%s' using 1:2 w boxes fs solid %f notitle\n", 
	      filename, fill_level);
    }else{
      fprintf(gp, "replot '%s' using 1:%d w boxes fs solid %f notitle\n",
	      filename, i+2, fill_level);
    }
  }

  fprintf(gp, "k = %lf\n", pow(mean, 2.0)/pow(SD, 2.0));
  fprintf(gp, "theta = %lf\n", pow(SD, 2.0)/mean);
  fprintf(gp, "f(x)=%lf*(x**(k-1))*exp(-x/theta)/gamma(k)/(theta**k)\n", interval);
  fprintf(gp, "replot f(x) w l lw 2 notitle\n");

  ChangeExtentionTo(fname_plt, ".plt");
  fprintf(gp, "save '%s'\n", fname_plt);

  if(eps_output_index == 1){
    eps_output(gp, fname_plt);
  }

  fflush(gp);
  getchar();
  pclose(gp);
}

void plot_histogram_wWeibullFit(char filename[], char xlabel[], char ylabel[], 
				float fill_level, int eps_output_index, 
				int datanum,
				double shape, double scale, double interval)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }

  gpStandardFormat(gp);
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);

  for(i=0;i<datanum;i++){
    if(i==0){
      fprintf(gp, "plot '%s' using 1:2 w boxes fs solid %f notitle\n", 
	      filename, fill_level);
    }else{
      fprintf(gp, "replot '%s' using 1:%d w boxes fs solid %f notitle\n",
	      filename, i+2, fill_level);
    }
  }

  fprintf(gp, "scale = %lf\n", scale);
  fprintf(gp, "shape = %lf\n", shape);
  //fprintf(gp, "f(x)=%lf/((2.0*3.14159)**0.5)/s*exp(-((x-m)**2.0)/2.0/(s**2.0))\n", interval);
  fprintf(gp, "f(x) = %lf*(shape/scale)*((x/scale)**(shape-1.0))*exp(-((x/scale)**(shape)))\n", interval);
  fprintf(gp, "replot f(x) w l lw 2 notitle\n");

  ChangeExtentionTo(fname_plt, ".plt");
  fprintf(gp, "save '%s'\n", fname_plt);

  if(eps_output_index == 1){
    eps_output(gp, fname_plt);
  }

  fflush(gp);
  getchar();
  pclose(gp);
}


FILE *gpopen(void)
{
  FILE *gp;
  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }
  return gp;
}

FILE *gpopen_nopersist(void)
{
  FILE *gp;
  gp = popen(GNUPLOT_PATH_NoPersist, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }
  return gp;
}

void gpclose(FILE *gp)
{
  fflush(gp);
  getchar();
  pclose(gp);
}

void gpplot(FILE *gp, char fname[], char argument[])
{
  fprintf(gp, "plot '%s' %s\n", fname, argument);
}

void gpcommand(FILE *gp, char argument[])
{
  fprintf(gp, "%s\n", argument);
}

void gpStandardFormat(FILE *gp)
{
  fprintf(gp, "set size ratio 0.618\n");
  fprintf(gp, "set xtics out\n");
  fprintf(gp, "set xtics nomirror\n");
  fprintf(gp, "set ytics out\n");
  fprintf(gp, "set ytics nomirror\n");
  fprintf(gp, "set border 3\n");
}

void eps_output(FILE *gp, char fname[])
{
  char fname_eps[256];

  strcpy(fname_eps, fname);
  ChangeExtentionTo(fname_eps, ".eps");

  fprintf(gp, "set term postscript eps enhanced color\n");
  fprintf(gp, "set output '%s'\n", fname_eps);
  fprintf(gp, "replot\n");
  //fprintf(gp, "set term x11\n");
  fprintf(gp, "set term qt\n"); //for qt
  fprintf(gp, "replot\n");
}

void gpsave(FILE *gp, char fname[])
{
  char fname_eps[256];

  strcpy(fname_eps, fname);
  ChangeExtentionTo(fname_eps, ".plt");

  fprintf(gp, "save '%s'\n", fname_eps);
}

void gpAxisLabels(FILE *gp, char xlabel[], char ylabel[])
{
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);
}

void ChangeExtentionTo(char fname[], char ext[])
{
  int i;
  char fname_temp[1000];

  for(i=0; i<strlen(fname); i++){
    if(fname[i] == '.'){
      fname[i] = '\0';
    }
  }

  sprintf(fname_temp, "%s%s", fname, ext);
  strcpy(fname, fname_temp);

  //printf("File name changed to '%s'\n", fname);
}

void gpLabels(FILE *gp, char xlabel[], char ylabel[])
{
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);
}

/************************
    FILE COPY
 ***********************/
void FileCopy(char fname_r[], char fname_w[])
{
  FILE *fpr;
  FILE *fpw;

  unsigned char buf[1000000];
  int i, size;

  fpr = fopen(fname_r, "rb");
  if(fname_r == NULL){
    printf("No binary file to read %s exists!\n", fname_r);
    exit(EXIT_FAILURE);
  }
  fpw = fopen(fname_w, "wb");
  if(fname_w == NULL){
    printf("Cannot open a file %s to write!\n", fname_w);
    exit(EXIT_FAILURE);
  }

  size = fread(buf, sizeof(unsigned char), 1000000, fpr);
  fwrite(buf, sizeof(unsigned char), size, fpw);

  fclose(fpr);
  fclose(fpw);

}

void FileContentShow(char fname[])
{
  FILE *fin;
  char buff[256];

  fin = fopen(fname, "r");
  if(fin == NULL){
    printf("ERROR: '%s' not exist!\n", fname);
  }else{
    while(fgets(buff, 256, fin)!=NULL){
      printf("%s", buff);
    }
  }
  fclose(fin);
}

long FileLineNumberReturn(char fname[])
{
  long counter = 0;
  FILE *fin;
  char buff[256];

  fin = fopen(fname, "r");
  if(fin == NULL){
    printf("ERROR: '%s' not exist!\n", fname);
  }else{
    while(fgets(buff, 256, fin)!=NULL){
      counter++;
    }
  }
  fclose(fin);

  return counter;
}

int FileExistCheck(char fname[])
{
  int returnValue;
  FILE *fin;
  fin = fopen(fname, "r");
  if(fin == NULL){
    returnValue = 0;
  }else{
    returnValue = 1;
  }
  fclose(fin);
  return returnValue;
}

void DataSort_LowToHigh(double value[], long k)
{
  long i, j, pos;
  double temp;

  for(i=0; i<k-1; i++){
    temp = value[i];
    pos = i;
    for(j=i+1; j<k; j++){
      if(temp > value[j]){
	temp = value[j];
	pos = j;
      }
    }
    value[pos] = value[i];
    value[i] = temp;
  }
}

void DataSort_LowToHigh_DoubleArray(double value[], double value2[], long k)
{
  long i, j, pos;
  double temp, temp2;

  for(i=0; i<k-1; i++){
    temp = value[i];
    temp2 = value2[i];
    pos = i;
    for(j=i+1; j<k; j++){
      if(temp > value[j]){
	temp = value[j];
	temp2 = value2[j];
	pos = j;
      }
    }
    value[pos] = value[i];
    value2[pos] = value2[i];
    value[i] = temp;
    value2[i] = temp2;
  }
}

void DataSort_LowToHigh_File(char fname[])
{
  double *value;
  long lineNum;
  long i;
  FILE *fin, *fout;

  lineNum = FileLineNumberReturn(fname);
  value = (double *)malloc(sizeof(double)*lineNum);
  if(value == NULL){
    printf("Pointer couldn't be assigned to value!\n");
    exit(1);
  }

  fin = fopen(fname, "r");
  for(i=0; i<lineNum; i++){
    fscanf(fin, "%lf", &value[i]);
  }
  fclose(fin);

  DataSort_LowToHigh(value, lineNum);

  fout = fopen(fname, "w");
  for(i=0; i<lineNum; i++){
    fprintf(fout, "%lf\n", value[i]);
  }
  fclose(fout);
  free(value);
}

void DataSort_HighToLow_File(char fname[])
{
  double *value;
  long lineNum;
  long i;
  FILE *fin, *fout;

  lineNum = FileLineNumberReturn(fname);
  value = (double *)malloc(sizeof(double)*lineNum);
  if(value == NULL){
    printf("Pointer couldn't be assigned to value!\n");
    exit(1);
  }

  fin = fopen(fname, "r");
  for(i=0; i<lineNum; i++){
    fscanf(fin, "%lf\n", &value[i]);
  }
  fclose(fin);

  DataSort_HighToLow(value, lineNum);

  fout = fopen(fname, "w");
  for(i=0; i<lineNum; i++){
    fprintf(fout, "%lf", value[i]);
  }
  fclose(fout);
  free(value);
}

void DataSort_LowToHigh_Long(long value[], long k)
{
  long i, j, pos;
  long temp;

  for(i=0; i<k-1; i++){
    temp = value[i];
    pos = i;
    for(j=i+1; j<k; j++){
      if(temp > value[j]){
	temp = value[j];
	pos = j;
      }
    }
    value[pos] = value[i];
    value[i] = temp;
  }
}

void DataSort_LowToHigh_Int(int value[], long k)
{
  long i, j, pos;
  int temp;

  for(i=0; i<k-1; i++){
    temp = value[i];
    pos = i;
    for(j=i+1; j<k; j++){
      if(temp > value[j]){
	temp = value[j];
	pos = j;
      }
    }
    value[pos] = value[i];
    value[i] = temp;
  }
}

void DataSort_HighToLow(double value[], long k)
{
  long i, j, pos;
  double temp;

  for(i=0; i<k-1; i++){
    temp = value[i];
    pos = i;
    for(j=i+1; j<k; j++){
      if(temp < value[j]){
	temp = value[j];
	pos = j;
      }
    }
    value[pos] = value[i];
    value[i] = temp;
  }
}

void DataSort_HighToLow_Long(long value[], long k)
{
  long i, j, pos;
  long temp;

  for(i=0; i<k-1; i++){
    temp = value[i];
    pos = i;
    for(j=i+1; j<k; j++){
      if(temp < value[j]){
	temp = value[j];
	pos = j;
      }
    }
    value[pos] = value[i];
    value[i] = temp;
  }
}

double PercentileValueReturn(double c[], long k, double percentile)
{
  long i;
  double *value;
  double pos;
  double returnValue;

  value = (double *)malloc(sizeof(double)*k);
  if(value == NULL){
    printf("'value' pointer couldn't be assigned!\n");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<k; i++){
    value[i] = c[i];
  }

  DataSort_LowToHigh(value, k);
  pos = percentile * k;

  returnValue = value[(long)pos];
  free(value);

  return returnValue;

}

double PercentileValueReturnFile(char fname[], double percentile)
{
  long lineNum;
  double valueTemp;
  double *value;
  double returnValue;
  FILE *fin;

  lineNum = FileLineNumberReturn(fname);
  value = (double *)malloc(sizeof(double)*lineNum);
  if(value == NULL){
    printf("'value' pointer couldn't be assigned\n");
    exit(EXIT_FAILURE);
  }

  fin = fopen(fname, "r");
  lineNum = 0;
  while(fscanf(fin, "%lf", &valueTemp) != EOF){
    value[lineNum] = valueTemp;
    lineNum++;
  }
  fclose(fin);

  returnValue = PercentileValueReturn(value, lineNum, percentile);
  free(value);

  return returnValue;
}

double MinValueReturnFile(char fname[])
{
  double valueTemp;
  double valueMin;
  FILE *fin;

  fin = fopen(fname, "r");

  if(fscanf(fin, "%lf", &valueMin)!=EOF){ //assigning initial value to min.
    while(fscanf(fin, "%lf", &valueTemp)!=EOF){
      if(valueTemp < valueMin){
	valueMin = valueTemp;
      }
    }
  }

  fclose(fin);
  return valueMin;
}

double MaxValueReturnFile(char fname[])
{
  double valueTemp;
  double valueMax;
  FILE *fin;

  fin = fopen(fname, "r");
  if(fscanf(fin, "%lf", &valueMax)!=EOF){ //assigning initial value to min.
    while(fscanf(fin, "%lf", &valueTemp)!=EOF){
      if(valueTemp > valueMax){
	valueMax = valueTemp;
      }
    }
  }
  fclose(fin);
  return valueMax;
}

void InputPositiveDoubleValue(double *value, char sentence[])
{
  while(1){
    printf("%s: ", sentence);
    scanf("%lf", value);
    if((*value)<=0.0){
      printf("Enter positive value!\n");
    }else{
      break;
    }
  }
}

void InputNonNegativeDoubleValue(double *value, char sentence[])
{
  while(1){
    printf("%s: ", sentence);
    scanf("%lf", value);
    if((*value)<0.0){
      printf("Enter non-negative value!\n");
    }else{
      break;
    }
  }
}

void InputPositiveLongValue(long *value, char sentence[])
{
  while(1){
    printf("%s: ", sentence);
    scanf("%ld", value);
    if((*value)<=0){
      printf("Enter positive value!\n");
    }else{
      break;
    }
  }
}

void InputNonNegativeLongValue(long *value, char sentence[])
{
  while(1){
    printf("%s: ", sentence);
    scanf("%ld", value);
    if((*value)<0){
      printf("Enter non-negative value!\n");
    }else{
      break;
    }
  }
}

void InputPositiveIntValue(int *value, char sentence[])
{
  while(1){
    printf("%s: ", sentence);
    scanf("%d", value);
    if((*value)<=0){
      printf("Enter positive value!\n");
    }else{
      break;
    }
  }
}

void InputNonNegativeIntValue(int *value, char sentence[])
{
  while(1){
    printf("%s: ", sentence);
    scanf("%d", value);
    if((*value)<0){
      printf("Enter positive value!\n");
    }else{
      break;
    }
  }
}


/************************
    STRING OPERATION
 ***********************/
void ExtractStringFromEnd(char result[], char original[], int sampleNo)
{
  int length;
  int i;

  length = strlen(original);
  if(length < sampleNo){
    printf("ERROR in 'ExtractStringFromEnd()'\n");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<sampleNo; i++){
    result[i] = original[length-sampleNo+i];
  }
  result[sampleNo] = '\0';
}

/************************
 STRING HANDLER
 ***********************/
void latterStringReturn(char name[], char namesub[], int charnum)
{
  int i;
  if(strlen(name)>charnum){
    for(i=charnum; i<=strlen(name); i++){
      namesub[i-charnum] = name[i];
    }
  }
}

int numberReturnFromString(char sentence[], int charnum)
{
  char frameString[250];
  int i;
  if(strlen(sentence)>charnum){
    for(i=charnum; i<strlen(sentence); i++){
      frameString[i-charnum] = sentence[i];
    }
  }
  return atoi(frameString);
}


