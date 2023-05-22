#include "PDDSparseGM.hpp"
#include "Monegros_Poisson.hpp"
#include "rectangle.hpp"
//scp -r PDDSparse fbernal0@login.g100.cineca.it:/g100/home/userexternal/fbernal0/LoopStudy
int main(int argc, char *argv[]){
   std::string config("configuration.txt");
    PDDSparseGM PDDS(argc,argv,config);
    BVP bvp;
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, pfscalarN> scalarN_init;
    std::map<std::string, std::string> string_init;
    //BVP initialization 
    scalar_init["f"] = Equation_f;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    vector_init["b"] = Equation_b;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Rectangle2D, Stopping);
    bvp.BVP_init(2,scalar_init,scalarN_init,vector_init, matrix_init,string_init, Equation_RBF);
    //PDDS.Solve(bvp);
    PDDS.Solve_Subdomains(bvp);
    std::vector<double> h_vec;
    std::vector<int> N_vec(11277,2000);
    char *line_buf = NULL;
    size_t line_buf_size = 0;
    int line_count = 0;
    ssize_t line_size;
    FILE *fp = fopen("h.txt", "r");
    if (!fp)
    {
        fprintf(stderr, "Error opening h file\n");
    }
    /* Get the first line of the file. */
    line_size = getline(&line_buf, &line_buf_size, fp);
    /* Loop through until we are done with the file. */
    while (line_size > 0)
    {
        /* Increment our line count */
        line_count++;
        /* Show the line details */
        h_vec.push_back(atof(line_buf));
        /* Get the next line */
        line_size = getline(&line_buf, &line_buf_size, fp);
    }
    /* Free the allocated line buffer */
    free(line_buf);
    line_buf = NULL;
    fclose(fp);
    //PDDS.Compute_h_N(bvp,0.005,h_vec,N_vec);
    PDDS.Solve_NumVR(bvp, h_vec,N_vec);
    MPI_Finalize();
    return 0;
}