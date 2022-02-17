# c_cpp_projects

## Genetic algorithm
- Simple implementation of GA in c.
- This was used for a POC, to optimize part position for 3d printing to minimize support volume  

To use:
- Include p3d_genetic.h in your project includes.
- Define the initial parameters
- Define the objective function
- Define the function to generate initial population (Optional)
- Define the lower and upper bounds of variables in gene_range
```
GeneticParams ga;
int gene_len = 2;
double gene_range[2 * gene_len] = {0, 100, 0, 100};
def define_and_run_ga() {
  init_genetic_params(&ga);
  def_genetic_params(&ga, 
    /* Population */         200,
    /*gene len    */         gene_len,      
    /* crs_rate*/            0.60, 
    /* mut_rate   */         0.20,     
    /*Cros Factor */         -1, 
    /* mut_fac    */         0.25, 
    /* conv_rate  */         0.95,    
    /* iters      */         200, 
    /* conv_cri   */         -1,      
    /* gene_range */         gene_range,
    /* Minimize/Maximize  */ P3D_MIN, 
    /* eval obj f */ &eval_obj_function,
    /* gen_init_pop */ NULL
    );
  run_genetic_algorithm(&ga);
}
```