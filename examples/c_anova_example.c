// c_anova_example.c

#include <stdio.h>
#include "measurements.h"

int main() {
    // Local Variables
    const int nsamples = 5;
    const int ncategories = 4;
    double x[] = {
        8.0, 9.0, 6.0, 7.0, 3.0,
        2.0, 4.0, 3.0, 5.0, 1.0,
        3.0, 5.0, 4.0, 2.0, 3.0,
        2.0, 2.0, -1.0, 0.0, 3.0
    };
    anova_table tbl;

    // Perform the ANOVA
    c_anova(nsamples, ncategories, x, nsamples, &tbl);

    // Display the results
    printf("Between Category Results:\n");
    printf("\tDOF: %i\n\tSSQ: %f\n\tMSQ: %f\n\tF Stat: %f\n\tProbability: %f\n",
        tbl.between.dof, 
        tbl.between.sum_of_squares, 
        tbl.between.mean_of_squares,
        tbl.between.f_stat,
        tbl.between.probability
    );
    printf("Residual Results:\n");
    printf("\tDOF: %i\n\tSSQ: %f\n\tMSQ: %f\n", 
        tbl.residual.dof,
        tbl.residual.sum_of_squares,
        tbl.residual.mean_of_squares
    );
    printf("Total Results:\n");
    printf("\tDOF: %i\n\tSSQ: %f\n\tMSQ: %f\n\tMean: %f\n", 
        tbl.total.dof,
        tbl.total.sum_of_squares,
        tbl.total.mean_of_squares,
        tbl.total.mean
    );

    // End
    return 0;
}