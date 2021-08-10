/***************************
 LineageAnalysis
 **************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>

#include "additionalfunctions.h"

double acorr;

void GammaIDTDistribution(void);
void ExportGammaLineageData_PreRun(double k, double theta, FILE *fout,
				   double currTime, double maxTime);
void ExportGammaLineageData(double k, double theta, FILE *fout,
			    double currTime, double maxTime,
			    int divCount);
void EffectOfRemoval(int index);
void GammaIDTDistributionWithDeath(double epsilon, double k,
				   double *popFit,
				   double *chronoFit, double *retroFit);
void ExportGammaLineageData_WithDeath_PreRun(double k, double theta, FILE *fout,
					     double currTime, double maxTime,
					     double epsilon);
void ExportGammaLineageData_WithDeath(double k, double theta, FILE *fout,
				      double currTime, double maxTime,
				      int divCount, double epsilon);
void ExportGammaTransitionData_WithDeath(double k, double theta, FILE *fout,
					 double currTime, double maxTime,
					 int divCount, double epsilon,
					 double prevIDT);
void ExportGammaTransitionData_WithDeath_PreRun(double k, double theta, FILE *fout,
						double currTime, double maxTime,
						double epsilon, double prevIDT);
void GammaTransitionWithRemoval(double epsilon, double k,
				double *popFit,
				double *chronoFit, double *retroFit);
void CalculatingCumulantGeneratingFunction(void);
void plot_delcum(void);
void plot_delcum_gzi(void);
double xPositionReturn(double *delcum, double checkValue);
void plot_legendre(void);
void plot_legendre_gzi(void);
void plot_cumulant_gzi(void);
void GziStatistics(int index);
void ExportGammaLineageData_counter(double k, double theta, long *counter,
				    double currTime, double maxTime);
void GammaIDTDistributionForSmallTreeSet(double k, long treeNum, long repeatNum, double popFit);
void GammaTransitionForSmallTreeSet(double k, long treeNum, long repeatNum, double popFit);
void plot_deviationProb(long treeNum);
void GammaTransitionWithRemoval_ForSmallTreeTest(double k, double *popFit,
						 double *chronoFit, double *retroFit);
void ExportGammaTransitionData_counter(double k, double theta, long *counter,
				       double currTime, double maxTime,
				       double prevIDT);

int main(){

  int choice;
  int loopindex = 0;

  while(loopindex == 0){
    printf("<LINEAGE ANALYSIS>\n");
    printf("1: Gamma IDT distribution\n");
    printf("2: Correlated gamma distrbution\n");
    printf("0: Exit\n");
    InputNonNegativeIntValue(&choice, "Enter number");
    switch(choice){
    case 0:
      loopindex = 1;
      break;
    case 1:
      EffectOfRemoval(0);
      break;
    case 2:
      InputPositiveDoubleValue(&acorr, "Input mother-daughter correlation");
      EffectOfRemoval(1);
      break;
    default:
      break;
    }
  }

  return 0;
}

void GziStatistics(int index) //index = 0: gamma; index = 1: gamma transition
{
  long i, j;
  FILE *fout, *gp, *fin;
  char fname_out[] = "CGF.dat";
  char fname_fcl[] = "Fcl.dat";
  char fname[] = "LineageData.dat";
  double popFit, chronoFit, retroFit;
  long treeNum = 1000; //This should be the same as the value in the pre-run functions. The value will be changed for the small tree simulation. 
  long counter;
  long divNum;
  double pcl[100], pgzi[100], prs[100];
  double cum[301], gzicum[301];
  double delcum[301], gzidelcum[301];
  double k;
  double value;
  double maxTime = 8.0;
  double popgzi;
  double sum, sumrs;

  InputPositiveDoubleValue(&k, "Enter shape parameter");  

  if(index == 0){
    GammaIDTDistributionWithDeath(0.0, k, &popFit, &chronoFit, &retroFit); //Simulates the results; Exports the data
  }else{
    GammaTransitionWithRemoval_ForSmallTreeTest(k, &popFit, &chronoFit, &retroFit);
  }
  
  printf("<h(D)>_{cl}/tau = %lf\n", chronoFit);
  printf("lambda = %lf\n", popFit);
  printf("<h(D)>_{rs}/tau = %lf\n", retroFit);

  //Calculating Pcl
  ColumnInitialize_double(pcl, 100);
  ColumnInitialize_double(pgzi, 100);
  ColumnInitialize_double(prs, 100);
  fin = fopen(fname, "r");
  while(fscanf(fin, "%ld", &divNum) != EOF){
    if(divNum < 100){
      pcl[divNum] = pcl[divNum] + pow(2.0, -divNum);
    }else{
      pcl[99] = pcl[99] + pow(2.0, -divNum);
    }
  }
  fclose(fin);
  
  for(i=0; i<100; i++){
    pcl[i] = pcl[i]/treeNum;
  }


  //Calculates cumulant generating function
  for(i=0; i<=300; i++){
    value = 0.0;
    for(j=0; j<100; j++){
      value = value + pow(2.0, ((i*0.01)-1.0)*j)*pcl[j];
    }
    cum[i] = log(value);
  }

  //Calculates the differential of cumulant generating function
  fout = fopen(fname_out, "w");
  for(i=1; i<300; i++){
    delcum[i] = (cum[i+1]-cum[i-1])/0.02;
    fprintf(fout, "%lf\t%lf\t%lf\t%lf\n", i*0.01-1.0, cum[i], delcum[i]/maxTime,
    	    (i*0.01-1.0)*delcum[i]-cum[i]);
  }
  fprintf(fout, "\n\n");
  fclose(fout);
  
  //Show s values for confirmation
  printf("Chronological position: %lf\n", xPositionReturn(delcum, chronoFit*maxTime));
  popgzi = xPositionReturn(delcum, popFit*maxTime);
  printf("Population growth rate position: %lf\n",
	 popgzi);
  printf("Retrospective position: %lf\n", xPositionReturn(delcum, retroFit*maxTime));

  //Assigning pgzi probability distribution
  for(i=0; i<100; i++){
    pgzi[i] = pow(2.0, i*popgzi)*pcl[i];
    prs[i] = pow(2.0, i)*pcl[i];
  }
  sum = ArraySum_double(pgzi, 100);
  sumrs = ArraySum_double(prs, 100);
  fout = fopen(fname_fcl, "w");
  for(i=0; i<100; i++){
    pgzi[i] = pgzi[i]/sum;
    prs[i] = prs[i]/sumrs;
    fprintf(fout, "%ld\t%lf\t%lf\t%lf\n", i, pcl[i], pgzi[i], prs[i]);
  }
  fclose(fout);

  //Calculating cum and delcum for gzi-statistics
  fout = fopen(fname_out, "a");
  for(i=0; i<=300; i++){
    value = 0.0;
    for(j=0; j<100; j++){
      value = value + pow(2.0, ((i*0.01)-1.0)*j)*pgzi[j];
    }
    gzicum[i] = log(value);
  }
  for(i=1; i<300; i++){
    gzidelcum[i] = (gzicum[i+1]-gzicum[i-1])/0.02;
    fprintf(fout, "%lf\t%lf\t%lf\t%lf\n", i*0.01-1.0, gzicum[i], gzidelcum[i]/maxTime,
    	    (i*0.01-1.0)*gzidelcum[i]-gzicum[i]);
  }
  fprintf(fout, "\n\n");
  fclose(fout);

  //Calculating cum and delcum for retrospective-statistics
  fout = fopen(fname_out, "a");
  for(i=0; i<=300; i++){
    value = 0.0;
    for(j=0; j<100; j++){
      value = value + pow(2.0, ((i*0.01)-1.0)*j)*prs[j];
    }
    gzicum[i] = log(value);
  }
  for(i=1; i<300; i++){
    gzidelcum[i] = (gzicum[i+1]-gzicum[i-1])/0.02;
    fprintf(fout, "%lf\t%lf\t%lf\t%lf\n", i*0.01-1.0, gzicum[i], gzidelcum[i]/maxTime,
    	    (i*0.01-1.0)*gzidelcum[i]-gzicum[i]);
  }
  fclose(fout);

  //Plot the distributions etc.
  plot_linespoints(fname_fcl, "lp pt 6", "Division count", "Probability", 0, 0, 1, 3);
  plot_cumulant_gzi();
  plot_delcum_gzi();
  plot_legendre_gzi();

  //Insert functions to show the deviations of population growth rate
  InputPositiveLongValue(&treeNum, "Enter number of initial cells (tentative: 10)");
  InputPositiveLongValue(&counter, "Enter repeat cycle (tentative: 10000)");
  if(index == 0){
    GammaIDTDistributionForSmallTreeSet(k, treeNum, counter, popFit);
  }else{
    GammaTransitionForSmallTreeSet(k, treeNum, counter, popFit);
  }
    
  plot_deviationProb(treeNum);  
}


void plot_deviationProb(long treeNum)
{
  char fname_out[] = "CGF.dat";
  char fname_result[] = "PopFitList.dat";
  FILE *gp;

  gp = gpopen();
  gpStandardFormat(gp);
  gpLabels(gp, "Empirical population growth rate", "1/n logP");
  fprintf(gp, "plot '%s' ind 0 using 3:(-$4) w l lw 2 ti 'gzi = 0'\n", fname_out);
  fprintf(gp, "replot '%s' ind 1 using 3:(-$4) w l lw 2 ti 'gzi = gzipop'\n", fname_out);
  fprintf(gp, "replot '%s' ind 2 using 3:(-$4) w l lw 2 ti 'gzi = 1'\n", fname_out);
  fprintf(gp, "replot '%s' using 1:2 w p pt 7 ps 2 ti 'n=%ld'\n", fname_result, treeNum);
  gpsave(gp, fname_result);
  eps_output(gp, fname_result);
  gpclose(gp);
}


void plot_cumulant_gzi(void)
{
  char fname_out[] = "CGF.dat";
  FILE *gp;

  gp = gpopen();
  gpStandardFormat(gp);
  gpLabels(gp, "s", "K(s)");
  fprintf(gp, "plot '%s' ind 0 w l lw 2 ti 'gzi = 0'\n", fname_out);
  fprintf(gp, "replot '%s' ind 1 w l lw 2 ti 'gzi = gzipop'\n", fname_out);
  fprintf(gp, "replot '%s' ind 2 w l lw 2 ti 'gzi = 1'\n", fname_out);
  gpsave(gp, fname_out);
  eps_output(gp, fname_out);
  gpclose(gp);
}


void CalculatingCumulantGeneratingFunction(void)
{
  long i, j;
  FILE *fout, *gp, *fin;
  char fname_out[] = "CGF.dat";
  char fname_fcl[] = "Fcl.dat";
  char fname[] = "LineageData.dat";
  double popFit, chronoFit, retroFit;
  long treeNum = 1000;
  long counter;
  long divNum;
  double pcl[100];
  double cum[301];
  double delcum[301];
  double k;
  double value;
  double maxTime = 8.0;

  InputPositiveDoubleValue(&k, "Enter shape parameter");  

  GammaIDTDistributionWithDeath(0.0, k, &popFit, &chronoFit, &retroFit); //Simulates the results; Exports the data
  printf("<h(D)>_{cl}/tau = %lf\n", chronoFit);
  printf("lambda = %lf\n", popFit);
  printf("<h(D)>_{rs}/tau = %lf\n", retroFit);

  //Calculating Pcl
  ColumnInitialize_double(pcl, 100);
  fin = fopen(fname, "r");
  while(fscanf(fin, "%ld", &divNum) != EOF){
    if(divNum < 100){
      pcl[divNum] = pcl[divNum] + pow(2.0, -divNum);
    }else{
      pcl[99] = pcl[99] + pow(2.0, -divNum);
    }
  }
  fclose(fin);

  fout = fopen(fname_fcl, "w");
  for(i=0; i<100; i++){
    pcl[i] = pcl[i]/treeNum;
    fprintf(fout, "%ld\t%lf\n", i, pcl[i]);
  }
  fclose(fout);

  for(i=0; i<=300; i++){
    value = 0.0;
    for(j=0; j<100; j++){
      value = value + pow(2.0, ((i*0.01)-1.0)*j)*pcl[j];
    }
    cum[i] = log(value);
  }

  fout = fopen(fname_out, "w");
  for(i=1; i<300; i++){
    delcum[i] = (cum[i+1]-cum[i-1])/0.02;
    fprintf(fout, "%lf\t%lf\t%lf\t%lf\n", i*0.01-1.0, cum[i], delcum[i]/maxTime,
	    (i*0.01-1.0)*delcum[i]-cum[i]);
  }
  fclose(fout);

  //Show s values for confirmation
  printf("Chronological position: %lf\n", xPositionReturn(delcum, chronoFit*maxTime));
  printf("Population growth rate position: %lf\n",
	 xPositionReturn(delcum, popFit*maxTime));
  printf("Retrospective position: %lf\n", xPositionReturn(delcum, retroFit*maxTime));

  plot_linespoints(fname_fcl, "p", "Division count", "Probability", 0, 0, 1, 1);
  plot_linespoints(fname_out, "l", "s", "K(s)", 0, 0, 1, 1);
  plot_delcum();
  plot_legendre();
}


double xPositionReturn(double *delcum, double checkValue)
{
  int i=2;
  double x1, x2;
  double y1, y2;

  while(delcum[i] < checkValue){
    i++;
    if(i == 299) break;
  }

  y1 = delcum[i-1];
  y2 = delcum[i];
  x1 = 0.01*(i-1)-1.0;
  x2 = 0.01*i-1.0;

  return (x2-x1)/(y2-y1)*checkValue + (x1*y2-x2*y1)/(y2-y1);  

}


void plot_delcum(void)
{
  char fname_out[] = "CGF.dat";
  char fname_export[] = "DelCGF.dat";
  FILE *gp;
  
  gp = gpopen();
  gpStandardFormat(gp);
  gpLabels(gp, "s", "dK(s)/ds/tau");
  fprintf(gp, "plot '%s' using 1:3 w l lw 2 notitle\n", fname_out);
  gpsave(gp, fname_export);
  eps_output(gp, fname_export);
  gpclose(gp);
}


void plot_delcum_gzi(void)
{
  char fname_out[] = "CGF.dat";
  char fname_export[] = "DelCGF.dat";
  FILE *gp;
  
  gp = gpopen();
  gpStandardFormat(gp);
  gpLabels(gp, "s", "dK(s)/ds/tau");
  fprintf(gp, "plot '%s' ind 0 using 1:3 w l lw 2 ti 'gzi = 0'\n", fname_out);
  fprintf(gp, "replot '%s' ind 1 using 1:3 w l lw 2 ti 'gzi = gzipop'\n", fname_out);
  fprintf(gp, "replot '%s' ind 2 using 1:3 w l lw 2 ti 'gzi = 1'\n", fname_out);
  gpsave(gp, fname_export);
  eps_output(gp, fname_export);
  gpclose(gp);
}


void plot_legendre(void)
{
  char fname_out[] = "CGF.dat";
  char fname_export[] = "Legendre.dat";
  FILE *gp;
  
  gp = gpopen();
  gpStandardFormat(gp);
  gpLabels(gp, "p", "G(p)");
  fprintf(gp, "plot '%s' using 3:4 w l lw 2 notitle\n", fname_out);
  gpsave(gp, fname_export);
  eps_output(gp, fname_export);
  gpclose(gp);
}


void plot_legendre_gzi(void)
{
  char fname_out[] = "CGF.dat";
  char fname_export[] = "Legendre.dat";
  FILE *gp;
  
  gp = gpopen();
  gpStandardFormat(gp);
  gpLabels(gp, "p", "G(p)");
  fprintf(gp, "plot '%s' ind 0 using 3:4 w l lw 2 ti 'gzi = 0'\n", fname_out);
  fprintf(gp, "replot '%s' ind 1 using 3:4 w l lw 2 ti 'gzi = gzipop'\n", fname_out);
  fprintf(gp, "replot '%s' ind 2 using 3:4 w l lw 2 ti 'gzi = 1'\n", fname_out);
  gpsave(gp, fname_export);
  eps_output(gp, fname_export);
  gpclose(gp);
}


void EffectOfRemoval(int index)
{
  int i;
  FILE *fout, *gp;
  char fname_out[] = "EffectOfRemoval_GammaIDT.dat";
  double popFit, chronoFit, retroFit;
  double epsilon;
  double noDeathPopFit;
  double noDeathRetroFit;
  double noDeathChronoFit;
  double k;

  InputPositiveDoubleValue(&k, "Enter shape parameter");  

  fout = fopen(fname_out, "w");
  for(i=0; i<=20; i++){
    epsilon = i*0.01;
    printf("Calculating epsilon = %lf\n", epsilon);
    if(index == 0){
      GammaIDTDistributionWithDeath(epsilon, k, &popFit, &chronoFit, &retroFit);
    }else{
      GammaTransitionWithRemoval(epsilon, k, &popFit, &chronoFit, &retroFit);
    }
    if(i==0){
      noDeathPopFit = popFit;
      noDeathRetroFit = retroFit;
      noDeathChronoFit = chronoFit;
    }
    fprintf(fout, "%lf\t%lf\t%lf\t%lf\t%lf\n",
    	    epsilon, popFit, chronoFit, retroFit, (popFit-noDeathPopFit)/noDeathPopFit);
  }
  fclose(fout);

  //Show plot in gnuplot
  gp = gpopen();
  gpStandardFormat(gp);
  gpLabels(gp, "epsilon", "Delta Lambda / Lambda");
  fprintf(gp, "set size ratio 1\n");
  fprintf(gp, "f(x, k) = -k*x\n");
  fprintf(gp, "plot f(x, %lf) w l lt 1 lc rgb 'red' ti 'Retrospective'\n", noDeathRetroFit/noDeathPopFit);
  fprintf(gp, "replot f(x, %lf) w l lt 1 lc rgb 'blue' ti 'Population'\n", noDeathPopFit/noDeathPopFit);
  fprintf(gp, "replot f(x, %lf) w l lt 1 lc rgb 'green' ti 'Chronological'\n", noDeathChronoFit/noDeathPopFit);
  fprintf(gp, "replot '%s' using 1:5 w p pt 7 ps 2 lc rgb 'magenta' notitle\n", fname_out);
  gpsave(gp, fname_out);
  eps_output(gp, fname_out);
  gpclose(gp);
}

void GammaIDTDistributionForSmallTreeSet(double k, long treeNum, long repeatNum, double popFit)
{
  long i, j;
  double theta;
  double maxTime = 8.0;
  char fname_out[] = "InitialLineageData.dat";
  char fname_result[] = "PopFitList.dat";
  double *initialDivisionTimeList;
  long counter, endCounter;
  long pos;
  double *popFitList;
  FILE *fout;
  int divCount;

  theta = pow(2.0, 1.0/k)-1.0;
  popFitList = (double *)malloc(sizeof(double)*repeatNum);
  if(popFitList == NULL){
    printf("popFitList pointer couldn't be assigned!\n");
    exit(1);
  }

  //Pre-run to get the information of initial division time
  printf("Pre-run is ongoing...\n");
  fout = fopen(fname_out, "w");
  for(i=0; i<1000; i++){
    ExportGammaLineageData_WithDeath_PreRun(k, theta, fout, 0.0, maxTime, 0.0);
  }
  fclose(fout);
  counter = FileLineNumberReturn(fname_out);

  //Get the information of initial division time
  initialDivisionTimeList = (double *)malloc(sizeof(double)*counter);
  if(initialDivisionTimeList == NULL){
    printf("initialDivisionTimeList pointer couldn't be assigined!\n");
    exit(1);
  }
  //Take the values to the array
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%lf", &initialDivisionTimeList[i]);
  }
  fclose(fout);

  for(i=0; i<repeatNum; i++){
    endCounter = 0;
    for(j=0; j<treeNum; j++){
      pos = (long)(genrand_real3()*counter);
      if(pos == counter) pos = counter -1;
      ExportGammaLineageData_counter(k, theta, &endCounter, initialDivisionTimeList[pos],
				     maxTime);
    }
    popFitList[i] = 1.0/maxTime*log(((double)endCounter)/((double)treeNum));
    if((i+1)%100 == 0) printf("%ld cycle finished\n", i+1);
  }

  fout = fopen(fname_result, "w");
  for(i=0; i<=1000; i++){
    counter = 0;
    for(j=0; j<repeatNum; j++){
      if(popFitList[j]>=(popFit+i*0.01)) counter = counter + 1;
    }
    if(counter > 0){
      fprintf(fout, "%lf\t%lf\n", popFit+i*0.01, log(((double)counter)/((double)repeatNum))/treeNum);
    }
  }
  fclose(fout);
  
  free(initialDivisionTimeList);
  free(popFitList);
}

void GammaTransitionForSmallTreeSet(double k, long treeNum, long repeatNum, double popFit)
{
  long i, j;
  double theta;
  double maxTime = 8.0;
  char fname_out[] = "InitialLineageData.dat";
  char fname_result[] = "PopFitList.dat";
  double *initialDivisionTimeList;
  double *initialIDTList;
  long counter, endCounter;
  long pos;
  double *popFitList;
  FILE *fout;
  int divCount;

  theta = pow(2.0, 1.0/k)-1.0;
  popFitList = (double *)malloc(sizeof(double)*repeatNum);
  if(popFitList == NULL){
    printf("popFitList pointer couldn't be assigned!\n");
    exit(1);
  }

  counter = FileLineNumberReturn(fname_out);

  //Get the information of initial division time
  initialDivisionTimeList = (double *)malloc(sizeof(double)*counter);
  initialIDTList = (double *)malloc(sizeof(double)*counter);
  if(initialDivisionTimeList == NULL){
    printf("initialDivisionTimeList or initialDivisionTimeList pointer couldn't be assigined!\n");
    exit(1);
  }
  //Take the values to the array
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%lf\t%lf", &initialDivisionTimeList[i], &initialIDTList[i]);
  }
  fclose(fout);

  for(i=0; i<repeatNum; i++){
    endCounter = 0;
    for(j=0; j<treeNum; j++){
      pos = (long)(genrand_real3()*counter);
      if(pos == counter) pos = counter -1;

      ExportGammaTransitionData_counter(k, theta, &endCounter, initialDivisionTimeList[pos],
					maxTime, initialIDTList[pos]);

    }
    popFitList[i] = 1.0/maxTime*log(((double)endCounter)/((double)treeNum));
    if((i+1)%100 == 0) printf("%ld cycle finished\n", i+1);
  }

  fout = fopen(fname_result, "w");
  for(i=0; i<=1000; i++){
    counter = 0;
    for(j=0; j<repeatNum; j++){
      if(popFitList[j]>=(popFit+i*0.01)) counter = counter + 1;
    }
    if(counter > 0){
      fprintf(fout, "%lf\t%lf\n", popFit+i*0.01, log(((double)counter)/((double)repeatNum))/treeNum);
    }
  }
  fclose(fout);
  
  free(initialDivisionTimeList);
  free(initialIDTList);
  free(popFitList);
}


void GammaIDTDistributionWithDeath(double epsilon, double k,
				   double *popFit,
				   double *chronoFit, double *retroFit)
{
  long treeNum = 1000;
  long i;
  FILE *fout;
  char fname_out[] = "LineageData.dat";
  double theta;
  double maxTime = 8.0; // duration of 8 mean generation time
  double *initialDivisionTimeList;
  long counter;
  long pos;
  int divCount;

  theta = pow(2.0, 1.0/k)-1.0;

  //Pre-run to get the information of initial division time
  fout = fopen(fname_out, "w");
  for(i=0; i<treeNum; i++){
    ExportGammaLineageData_WithDeath_PreRun(k, theta, fout, 0.0, maxTime, epsilon);
  }
  fclose(fout);
  counter = FileLineNumberReturn(fname_out);
  printf("\n%ld pre-run data exported!\n\n", counter);

  //Get the information of initial division time
  initialDivisionTimeList = (double *)malloc(sizeof(double)*counter);
  if(initialDivisionTimeList == NULL){
    printf("initialDivisionTimeList pointer couldn't be assigined!\n");
    exit(1);
  }
  //Take the values to the column
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%lf", &initialDivisionTimeList[i]);
  }
  fclose(fout);

  //Simulating the lineage
  fout = fopen(fname_out, "w");
  for(i=0; i<treeNum; i++){
    pos = (long)(genrand_real3()*counter);
    if(pos == counter) pos = counter - 1;
    ExportGammaLineageData_WithDeath(k, theta, fout, initialDivisionTimeList[pos],
				     maxTime, 0, epsilon);
  }
  fclose(fout);    
  free(initialDivisionTimeList);

  //Check the exported data
  counter = FileLineNumberReturn(fname_out);
  printf("%ld Lineage Data Exported\n", counter);
  printf("Calculated population growth rate: %lf\n",
	 1.0/maxTime*log(((double)counter)/((double)treeNum)));
  (*popFit) = 1.0/maxTime*log(((double)counter)/((double)treeNum));

  (*chronoFit) = 0.0;
  (*retroFit) = 0.0;
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%d", &divCount);
    (*chronoFit) = (*chronoFit) + divCount*log(2.0)*pow(2.0, - divCount)/treeNum;
    (*retroFit) = (*retroFit) + divCount*log(2.0)/counter;
  }
  fclose(fout);

  (*chronoFit) = (*chronoFit)/maxTime;
  (*retroFit) = (*retroFit)/maxTime;

}

void GammaTransitionWithRemoval(double epsilon, double k,
				double *popFit,
				double *chronoFit, double *retroFit)
{
  long treeNum = 1000;
  long i;
  FILE *fout;
  char fname_out[] = "LineageData.dat";
  double theta;
  double maxTime = 8.0; // duration of 8 mean generation time
  double *initialDivisionTimeList, *initialIDT;
  long counter;
  long pos;
  int divCount;

  theta = pow(2.0, 1.0/k)-1.0;

  //Pre-run to get the information of initial division time
  fout = fopen(fname_out, "w");
  for(i=0; i<treeNum; i++){
    ExportGammaTransitionData_WithDeath_PreRun(k, theta, fout, 0.0, maxTime, epsilon, gamma_rand(k, theta));
  }
  fclose(fout);
  counter = FileLineNumberReturn(fname_out);
  printf("\n%ld pre-run data exported!\n\n", counter);

  //Get the information of initial division time
  initialDivisionTimeList = (double *)malloc(sizeof(double)*counter);
  initialIDT = (double *)malloc(sizeof(double)*counter);
  if(initialDivisionTimeList == NULL || initialIDT == NULL){
    printf("initialDivisionTimeList pointer couldn't be assigined!\n");
    exit(1);
  }
  //Take the values to the column
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%lf\t%lf", &initialDivisionTimeList[i], &initialIDT[i]);
  }
  fclose(fout);

  //Simulating the lineage
  fout = fopen(fname_out, "w");
  for(i=0; i<treeNum; i++){
    pos = (long)(genrand_real3()*counter);
    if(pos == counter) pos = counter - 1;
    ExportGammaTransitionData_WithDeath(k, theta, fout, initialDivisionTimeList[pos],
					maxTime, 0, epsilon, initialIDT[pos]);
  }
  fclose(fout);    
  free(initialDivisionTimeList);
  free(initialIDT);

  //Check the exported data
  counter = FileLineNumberReturn(fname_out);
  printf("%ld Lineage Data Exported\n", counter);
  printf("Calculated population growth rate: %lf\n",
	 1.0/maxTime*log(((double)counter)/((double)treeNum)));
  (*popFit) = 1.0/maxTime*log(((double)counter)/((double)treeNum));

  (*chronoFit) = 0.0;
  (*retroFit) = 0.0;
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%d", &divCount);
    (*chronoFit) = (*chronoFit) + divCount*log(2.0)*pow(2.0, - divCount)/treeNum;
    (*retroFit) = (*retroFit) + divCount*log(2.0)/counter;
  }
  fclose(fout);

  (*chronoFit) = (*chronoFit)/maxTime;
  (*retroFit) = (*retroFit)/maxTime;

}

void GammaTransitionWithRemoval_ForSmallTreeTest(double k, double *popFit,
						 double *chronoFit, double *retroFit)
{
  long treeNum = 1000;
  long i;
  FILE *fout;
  char fname_prerun[] = "InitialLineageData.dat";
  char fname_out[] = "LineageData.dat";
  double theta;
  double maxTime = 8.0; // duration of 8 mean generation time
  double *initialDivisionTimeList, *initialIDT;
  long counter;
  long pos;
  int divCount;

  theta = pow(2.0, 1.0/k)-1.0;

  //Pre-run to get the information of initial division time
  fout = fopen(fname_prerun, "w");
  for(i=0; i<treeNum; i++){
    ExportGammaTransitionData_WithDeath_PreRun(k, theta, fout, 0.0, maxTime, 0.0, gamma_rand(k, theta));
  }
  fclose(fout);
  counter = FileLineNumberReturn(fname_prerun);
  printf("\n%ld pre-run data exported!\n\n", counter);

  //Get the information of initial division time
  initialDivisionTimeList = (double *)malloc(sizeof(double)*counter);
  initialIDT = (double *)malloc(sizeof(double)*counter);
  if(initialDivisionTimeList == NULL || initialIDT == NULL){
    printf("initialDivisionTimeList pointer couldn't be assigined!\n");
    exit(1);
  }
  //Take the values to the column
  fout = fopen(fname_prerun, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%lf\t%lf", &initialDivisionTimeList[i], &initialIDT[i]);
  }
  fclose(fout);

  //Simulating the lineage
  fout = fopen(fname_out, "w");
  for(i=0; i<treeNum; i++){
    pos = (long)(genrand_real3()*counter);
    if(pos == counter) pos = counter - 1;
    ExportGammaTransitionData_WithDeath(k, theta, fout, initialDivisionTimeList[pos],
					maxTime, 0, 0.0, initialIDT[pos]);
  }
  fclose(fout);    
  free(initialDivisionTimeList);
  free(initialIDT);

  //Check the exported data
  counter = FileLineNumberReturn(fname_out);
  printf("%ld Lineage Data Exported\n", counter);
  printf("Calculated population growth rate: %lf\n",
	 1.0/maxTime*log(((double)counter)/((double)treeNum)));
  (*popFit) = 1.0/maxTime*log(((double)counter)/((double)treeNum));

  (*chronoFit) = 0.0;
  (*retroFit) = 0.0;
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%d", &divCount);
    (*chronoFit) = (*chronoFit) + divCount*log(2.0)*pow(2.0, - divCount)/treeNum;
    (*retroFit) = (*retroFit) + divCount*log(2.0)/counter;
  }
  fclose(fout);

  (*chronoFit) = (*chronoFit)/maxTime;
  (*retroFit) = (*retroFit)/maxTime;

}


void GammaIDTDistribution(void)
{
  long treeNum = 1000;
  long i;
  FILE *fout;
  char fname_out[] = "LineageData.dat";
  double k, theta;
  double maxTime = 8.0; // duration of 8 mean generation time
  double *initialDivisionTimeList;
  long counter;
  long pos;
  int divCount;
  double chronoFitness, retroFitness;

  InputPositiveDoubleValue(&k, "Enter shape parameter");
  theta = pow(2.0, 1.0/k)-1.0;
  printf("\n<Simulation conditions>\n");
  printf("Shape parameter: %lf\n", k);
  printf("Scale parameter: %lf\n", theta);
  printf("Mean chronological fitness: %lf\n", log(2.0)/(k*theta));
  printf("Predicted population growth rate: %lf\n",
	 (pow(2.0, 1.0/k) - 1.0)/theta);
  printf("Predicted mean retrospective fitness: %lf\n",
	 pow(2.0, 1.0/k)*log(2.0)/k/theta);

  //Pre-run to get the information of initial division time
  fout = fopen(fname_out, "w");
  ExportGammaLineageData_PreRun(k, theta, fout, 0.0, maxTime);
  fclose(fout);
  counter = FileLineNumberReturn(fname_out);
  printf("\n%ld pre-run data exported!\n\n", counter);

  //Get the information of initial division time
  initialDivisionTimeList = (double *)malloc(sizeof(double)*counter);
  if(initialDivisionTimeList == NULL){
    printf("initialDivisionTimeList pointer couldn't be assigined!\n");
    exit(1);
  }
  //Take the values to the column
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%lf", &initialDivisionTimeList[i]);
  }
  fclose(fout);

  //Simulating the lineage
  fout = fopen(fname_out, "w");
  for(i=0; i<treeNum; i++){
    pos = (long)(genrand_real3()*counter);
    if(pos == counter) pos = counter - 1;
    ExportGammaLineageData(k, theta, fout, initialDivisionTimeList[pos],
			   maxTime, 0);
  }
  fclose(fout);    
  free(initialDivisionTimeList);

  //Check the exported data
  counter = FileLineNumberReturn(fname_out);
  printf("%ld Lineage Data Exported\n", counter);
  printf("Calculated population growth rate: %lf\n",
	 1.0/maxTime*log(((double)counter)/((double)treeNum)));

  chronoFitness = 0.0;
  retroFitness = 0.0;
  fout = fopen(fname_out, "r");
  for(i=0; i<counter; i++){
    fscanf(fout, "%d", &divCount);
    chronoFitness = chronoFitness + divCount*log(2.0)*pow(2.0, - divCount)/treeNum;
    retroFitness = retroFitness + divCount*log(2.0)/counter;
  }
  fclose(fout);

  printf("Chronological mean fitness: %lf\n", chronoFitness/maxTime);
  printf("Retrospective mean fitness: %lf\n", retroFitness/maxTime);

}

void ExportGammaLineageData_PreRun(double k, double theta, FILE *fout,
				   double currTime, double maxTime)
{
  double idt;

  idt = gamma_rand(k, theta);
  if(currTime + idt >= maxTime){
    fprintf(fout, "%lf\n", currTime + idt - maxTime);
  }else{
    ExportGammaLineageData_PreRun(k, theta, fout, currTime + idt, maxTime); //Sister 1
    ExportGammaLineageData_PreRun(k, theta, fout, currTime + idt, maxTime); //Sister 2
  }
  
}

void ExportGammaLineageData_WithDeath_PreRun(double k, double theta, FILE *fout,
					     double currTime, double maxTime,
					     double epsilon)
{
  double idt;

  if(genrand_real3()<=pow(2.0, -epsilon)){ //Survived
    idt = gamma_rand(k, theta);
    if(currTime + idt >= maxTime){
      fprintf(fout, "%lf\n", currTime + idt - maxTime);
    }else{
      ExportGammaLineageData_WithDeath_PreRun(k, theta, fout, currTime + idt, maxTime, epsilon); //Sister 1
      ExportGammaLineageData_WithDeath_PreRun(k, theta, fout, currTime + idt, maxTime, epsilon); //Sister 2
    }
  }
  
}


void ExportGammaTransitionData_WithDeath_PreRun(double k, double theta, FILE *fout,
						double currTime, double maxTime,
						double epsilon, double prevIDT)
{
  double idt;

  if(genrand_real3()<=pow(2.0, -epsilon)){ //Survived
    idt = NextValueGamma(k, theta, acorr, prevIDT);
    if(currTime + idt >= maxTime){
      fprintf(fout, "%lf\t%lf\n", currTime + idt - maxTime, idt);
    }else{
      ExportGammaTransitionData_WithDeath_PreRun(k, theta, fout, currTime + idt, maxTime, epsilon, idt); //Sister 1
      ExportGammaTransitionData_WithDeath_PreRun(k, theta, fout, currTime + idt, maxTime, epsilon, idt); //Sister 2
    }
  }
  
}


void ExportGammaLineageData(double k, double theta, FILE *fout,
			    double currTime, double maxTime,
			    int divCount)
{
  double idt;

  idt = gamma_rand(k, theta);
  if(currTime + idt >= maxTime){
    fprintf(fout, "%d\n", divCount);
  }else{
    ExportGammaLineageData(k, theta, fout, currTime + idt, maxTime, divCount + 1); //Sister 1
    ExportGammaLineageData(k, theta, fout, currTime + idt, maxTime, divCount + 1); //Sister 2
  }
  
}

void ExportGammaLineageData_counter(double k, double theta, long *counter,
				    double currTime, double maxTime)
{
  double idt;

  idt = gamma_rand(k, theta);
  if(currTime + idt >= maxTime){
    (*counter) = (*counter) + 1;
  }else{
    ExportGammaLineageData_counter(k, theta, counter, currTime + idt, maxTime); //Sister 1
    ExportGammaLineageData_counter(k, theta, counter, currTime + idt, maxTime); //Sister 2
  }
  
}

void ExportGammaLineageData_WithDeath(double k, double theta, FILE *fout,
				      double currTime, double maxTime,
				      int divCount, double epsilon)
{
  double idt;

  if(genrand_real3()<=pow(2.0, -epsilon)){ //survived
    idt = gamma_rand(k, theta);
    if(currTime + idt >= maxTime){
      fprintf(fout, "%d\n", divCount);
    }else{
      ExportGammaLineageData_WithDeath(k, theta, fout, currTime + idt, maxTime, divCount + 1, epsilon); //Sister 1
      ExportGammaLineageData_WithDeath(k, theta, fout, currTime + idt, maxTime, divCount + 1, epsilon); //Sister 2
    }
  }
  
}


void ExportGammaTransitionData_WithDeath(double k, double theta, FILE *fout,
					 double currTime, double maxTime,
					 int divCount, double epsilon,
					 double prevIDT)
{
  double idt;

  if(genrand_real3()<=pow(2.0, -epsilon)){ //survived
    idt = NextValueGamma(k, theta, acorr, prevIDT);
    if(currTime + idt >= maxTime){
      fprintf(fout, "%d\n", divCount);
    }else{
      ExportGammaTransitionData_WithDeath(k, theta, fout, currTime + idt, maxTime, divCount + 1, epsilon, idt); //Sister 1
      ExportGammaTransitionData_WithDeath(k, theta, fout, currTime + idt, maxTime, divCount + 1, epsilon, idt); //Sister 2
    }
  }
    
}

void ExportGammaTransitionData_counter(double k, double theta, long *counter,
				       double currTime, double maxTime,
				       double prevIDT)
{
  double idt;

  idt = NextValueGamma(k, theta, acorr, prevIDT);
  if(currTime + idt >= maxTime){
    (*counter) = (*counter) + 1;
  }else{
    ExportGammaTransitionData_counter(k, theta, counter, currTime + idt, maxTime, idt); //Sister 1
    ExportGammaTransitionData_counter(k, theta, counter, currTime + idt, maxTime, idt); //Sister 2
  }    
}


