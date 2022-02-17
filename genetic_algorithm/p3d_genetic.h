
#include<stdio.h>
#include<conio.h>
#include<fstream>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include <iostream>
#include <vector>
#include <string>

//using namespace std;

FILE* get_gp_log();

#define CONV_RATE 0;
#define CONV_FIX_ITER 1;

/* Minimize or Maximize the Objective Function */
typedef enum p3d_min_max
{
  P3D_MIN = 0,
  P3D_MAX
}P3DMinMax;

typedef enum fitness_criteria
{
  FIT_SORT_BY_VALUE,
  FIT_PROBABILITY_SORT
}P3dFitnessCriteria;

typedef struct genetic_params GeneticParams;
struct genetic_params
{
	int Population;
  int Offspring;
  int GeneLength;

  double CrossOverFactor;
	double CrossRate;
	double MutRate;
  double MutFactor;
  double ConvergenceRate;
  int iters;
  int ConvCriteria;
  P3DMinMax MinMax;

  double *GeneRange;
	double *Parents;
	double *Children;
	double *FunctionValues;

  void (*eval_obj_funcs)(GeneticParams *gene, int num);
  void (*generate_initial_population)(GeneticParams *gene);
};

/* External Functions to use to define and run the GA*/
void run_genetic_algorithm(GeneticParams *gene);
void init_genetic_params(GeneticParams *gene);
void def_genetic_params(GeneticParams *gene, int pop, int gen_len, double crs_rate, 
  double mut_rate, double crs_fact, double mut_fac, double conv_rate, int iters, 
  int conv_cri, double gene_range[], P3DMinMax minmax,
  void (*ev_obj_f)(GeneticParams *gene, int num),
  void (*gen_init_pop)(GeneticParams *gene));
void clear_genetic_params(GeneticParams *gene);

/*Internal Functions*/
void generate_parents(GeneticParams *gene);
void do_crossover(GeneticParams *gene);
void do_mutation(GeneticParams *gene);
void do_selection_simple(GeneticParams *gene);
//void evaluate_fitness(GeneticParams *gene);
void evaluate_objective_func(GeneticParams *gene, int start, int end);
int check_convergence(GeneticParams *gene);

