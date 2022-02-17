#include "p3d_genetic.h"

FILE *gp = NULL;
FILE* get_gp_log(){
  return gp;
}
/*=============================================================================*/
void init_genetic_params(GeneticParams *gene)
{
  if(!gp)
    fopen_s(&gp,"gene.txt","w");
  gene->eval_obj_funcs = NULL;
  gene->generate_initial_population = NULL;

  gene->Population = 0;
  gene->Offspring = 0;
  gene->CrossRate = 0;
  gene->CrossOverFactor = 1;
  gene->MutRate = 0;
  gene->ConvergenceRate = 0;
  gene->iters = 0;
  gene->ConvCriteria = 0;
  gene->MinMax = P3D_MIN;
  gene->MutFactor = -1;

  gene->Parents = NULL;
  gene->Children = NULL;
  gene->FunctionValues = NULL;
  gene->GeneRange = NULL;
}

/*=============================================================================*/
void clear_genetic_params(GeneticParams *gene)
{
  fclose(gp);
  gene->eval_obj_funcs = NULL;
  gene->generate_initial_population = NULL;

  gene->Population = 0;
  gene->GeneLength = 0;
  gene->CrossRate = 0;
  gene->MutRate = 0;
  gene->Offspring = 0;
  gene->ConvergenceRate = 0;
  gene->iters = 0;
  gene->ConvCriteria = 0;
  gene->CrossOverFactor = 0;
  gene->MutFactor = -1;

  free(gene->GeneRange);
  free(gene->Parents);
  free(gene->Children);
  free(gene->FunctionValues);
}

/*=============================================================================*/
void def_genetic_params(GeneticParams *gene, int pop, int gen_len, double crs_rate, 
  double mut_rate, double crs_fact, double mut_fac, double conv_rate, int iters, 
  int conv_cri, double gene_range[], P3DMinMax minmax,
  void (*ev_obj_f)(GeneticParams *gene, int num),
  void (*gen_init_pop)(GeneticParams *gene))
{
  int ii;
  gene->eval_obj_funcs = ev_obj_f;
  gene->generate_initial_population = gen_init_pop;
  gene->Population = pop;
  gene->GeneLength = gen_len;
  gene->CrossRate = crs_rate;
  gene->MutRate = mut_rate;
  gene->Offspring = (int)(gene->CrossRate*(double)gene->Population);
  gene->ConvergenceRate = conv_rate;
  gene->iters = iters;
  if(crs_fact != -1)
    gene->CrossOverFactor = crs_fact;
  if(conv_cri != -1)
    gene->ConvCriteria = conv_cri;
  if(mut_fac != -1)
    gene->MutFactor = mut_fac;
  gene->MinMax = minmax;

  gene->GeneRange = (double*)calloc(2 * gene->GeneLength, sizeof(double));
  for(ii=0; ii<gene->GeneLength*2; ii++)
  {
    gene->GeneRange[ii] = gene_range[ii];
  }
  gene->Parents = (double*)calloc((gene->Population + gene->Offspring + 2) * gene->GeneLength, sizeof(double));
  gene->Children = (double*)calloc((gene->Offspring + 2) * gene->GeneLength, sizeof(double));
  gene->FunctionValues = (double*)calloc(gene->Population + gene->Offspring + 2, sizeof(double));
}

/*=============================================================================*/
void print_gnome(GeneticParams *gene, char* s, int num)
{
  int j;
  fprintf(gp,"%s %d: [ ",s,num);
  for(j=0;j<gene->GeneLength;j++)
    fprintf(gp,"%.2lf,\t",gene->Parents[num*gene->GeneLength+j]);
  fprintf(gp,"]\tVal: %lf\n",gene->FunctionValues[num]);
}

/*=============================================================================*/
void print_gene_limits(GeneticParams *gene)
{
  int j;
  fprintf(gp,"Gene Ranges: \n");
  for(j=0;j<gene->GeneLength;j++)
  {
    fprintf(gp,"%.2lf,\t",gene->GeneRange[j*2]);
  }
  fprintf(gp,"\n");
  for(j=0;j<gene->GeneLength;j++)
  {
    fprintf(gp,"%.2lf,\t",gene->GeneRange[j*2+1]);
  }
  fprintf(gp,"\n");

}

/*=============================================================================*/
void print_iteration_info(GeneticParams *gene, int itr)
{
  int i,j;
  fprintf(gp,"\niteration: %d\n",itr);
  for(i=0;i<gene->Population+gene->Offspring;i++)
  {
    if(i<gene->Population)
      fprintf(gp,"Parent %d: [ ",i);
    else
      fprintf(gp,"Child %d: [ ",i);
    for(j=0;j<gene->GeneLength;j++)
      fprintf(gp,"%.2lf,\t",gene->Parents[i*gene->GeneLength+j]);
    fprintf(gp,"]\tVal: %lf\n",gene->FunctionValues[i]);
  }
/*
  for(i=0;i<gene->Offspring;i++)
  {
    fprintf(gp,"Child %d: [ ",i);
    for(j=0;j<gene->GeneLength;j++)
      fprintf(gp,"%.2lf,\t",gene->Children[i*gene->GeneLength+j]);
    fprintf(gp,"]\tVal: %lf\n",gene->FunctionValues[gene->Population+i]);
  }
*/
}

/*=============================================================================*/
void run_genetic_algorithm(GeneticParams *gene)
{
	int i=-1;
	int b=0;
  print_gene_limits(gene);
  if(gene->generate_initial_population == NULL)
    generate_parents(gene);
  else
    gene->generate_initial_population(gene);
  
  print_iteration_info(gene, i);
  for( i=0;i<gene->iters;i++)
  {
    do_crossover(gene);
    print_iteration_info(gene, i);
    do_mutation(gene);
    do_selection_simple(gene);
//    print_iteration_info(gene, i);
    if(check_convergence(gene))
    {
      fprintf(gp,"CONVERGENCE\n");
//      break;
    }
  }
  print_iteration_info(gene, i);

}

/*=============================================================================*/
void generate_parents(GeneticParams *gene)
{
  int ii,jj;
  double val;
  srand ( (unsigned int)time(NULL) );

  for(ii=0; ii<gene->Population ; ii++)
  {
    for(jj=0; jj<gene->GeneLength ; jj++)
    {
      val = rand()%(int)(gene->GeneRange[2*jj+1] - gene->GeneRange[2*jj]);
      val = val + gene->GeneRange[2*jj];
//      gene->ParentsV.push_back(val);
      gene->Parents[ii*gene->GeneLength+jj] = val;
    }
  }
  evaluate_objective_func(gene, 0, gene->Population);
}

/*=============================================================================*/
void do_crossover(GeneticParams *gene)
{
	int i,j,k;
	int crospoint;
	double factor;

  factor = (double)gene->GeneLength*gene->CrossOverFactor;
  factor = (int)factor;
  if(factor == 0)
    factor++;

//  fprintf(gp,"crossover:\n");
  for(i=0; i<gene->Offspring; i=i+2)
  {
    crospoint = rand()%(int)factor;
    if(crospoint == 0)
      crospoint++;
    if(crospoint == gene->GeneLength)
      crospoint--;

//    fprintf(gp,"Point: %d\n",crospoint);
/*
    print_gnome(gene, "Parent", i*gene->GeneLength);
    print_gnome(gene, "Parent", (i+1)*gene->GeneLength);
*/
		for(k=0;k<crospoint;k++)
		{
      gene->Children[i*gene->GeneLength + k] = gene->Parents[i*gene->GeneLength + k];
			gene->Children[(i+1)*gene->GeneLength + k] = gene->Parents[(i+1)*gene->GeneLength + k];
		}
		for(k=crospoint;k<gene->GeneLength;k++)
		{
			gene->Children[i*gene->GeneLength + k] = gene->Parents[(i+1)*gene->GeneLength + k];
			gene->Children[(i+1)*gene->GeneLength + k] = gene->Parents[i*gene->GeneLength + k];
		}
  }
  for(i=gene->Population, k=0; i<gene->Population + gene->Offspring ; i++,k++)
  {
    for(j=0; j<gene->GeneLength ; j++)
      gene->Parents[i*gene->GeneLength + j] = gene->Children[k*gene->GeneLength + j];
  }
  evaluate_objective_func(gene, gene->Population, gene->Population + gene->Offspring);
}

/*=============================================================================*/
static double generate_population_value(GeneticParams *gene, int gene_place)
{
  double val;
  val = rand()%(int)(gene->GeneRange[2 * gene_place + 1] - gene->GeneRange[2*gene_place]);
  val = val + gene->GeneRange[2 * gene_place];
  return val;
}

/*=============================================================================*/
void do_mutation(GeneticParams *gene)
{
	int i,j;
	int mutspot, mutcount, mutnum, num_spots;
  double factor,val;

  // to restrict the genelength which gets mutated
  // CrossOverFactor = 1 - entire gene considered
  factor = (double)gene->GeneLength * gene->CrossOverFactor;
  factor = (int)factor;
  if(factor == 0)
    factor++;
  
  mutcount = (int)(gene->MutRate * gene->Population) + 1;

  for(i=0;i<mutcount;i++)
  {
    mutnum = rand()%(gene->Population + gene->Offspring);
    if(mutnum < (int).05*gene->Population )
      mutnum = mutnum + (int).05*gene->Population ;
    mutspot = rand()%(int)factor;
    num_spots = rand()%(int)(gene->GeneLength*gene->MutFactor+1);
//    fprintf(gp,"Mutate: %d\n",mutspot);
//      print_gnome(gene, "gnome",mutnum);
    for(j=0; j<num_spots; j++)
    {
      do
      {
        val = generate_population_value(gene, mutspot);
      }while(val == gene->Parents[mutnum*gene->GeneLength + mutspot]);
      gene->Parents[mutnum*gene->GeneLength + mutspot] = generate_population_value(gene, mutspot);
      evaluate_objective_func(gene, mutnum, mutnum+1);
//      print_gnome(gene, "gnome", mutnum);
      mutspot++;
      if(mutspot == gene->GeneLength)
        mutspot = 0;
    }
  }
}

/*=============================================================================*/
void do_selection_simple(GeneticParams *gene)
{
	double z,check,check1;
	int i,j,k;

	for(i=0;i<gene->Population + gene->Offspring;i++)
	{
		for(j=i+1;j<gene->Population + gene->Offspring;j++)
		{
			z = gene->FunctionValues[i] - gene->FunctionValues[j];
      if(gene->MinMax == P3D_MAX)
        z = -1*z;
			if(z>0)
			{
				check = gene->FunctionValues[i];
				gene->FunctionValues[i] = gene->FunctionValues[j];
				gene->FunctionValues[j] = check;
        for(k=0;k<gene->GeneLength;k++)
				{
          check1 = gene->Parents[i*gene->GeneLength+k];
					gene->Parents[i*gene->GeneLength+k] = gene->Parents[j*gene->GeneLength+k];
					gene->Parents[j*gene->GeneLength+k] = check1;
				}
			}
		}
	}
}

/*=============================================================================*/
int check_convergence(GeneticParams *gene)
{
	int i,count;
	count = 0;
	for(i=0;i<gene->Population;i++)
	{
		double Diff = abs(gene->FunctionValues[i] - gene->FunctionValues[i+1]);
		
		if(Diff < 0.0001)
			count++;
	}
	int per = count/gene->Population;
	per = per*100;
//	printf("\nPer:%3d",per);//_getch();

  if(per >= gene->ConvergenceRate)
		return(1);
	else
		return(0);

}

/*=============================================================================*/
void evaluate_objective_func(GeneticParams *gene, int start, int end)
{
  int i;
  for(i=start; i<end; i++)
  {
    gene->eval_obj_funcs(gene, i);
  }
}

