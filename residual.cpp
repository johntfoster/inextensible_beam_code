#include <iostream>
#include <math.h>

void evaluate_dN_u1_dxi(const double xi, double* dN_u1_dxi){
    dN_u1_dxi[0] = 1. / 2. * (-1 + xi) + xi / 2.; 
    dN_u1_dxi[1] = -2. * xi; 
    dN_u1_dxi[2] = xi / 2. + (1 + xi) / 2.;
    return;
}

void evaluate_N_lambda(const double xi, double* N_lambda){
    N_lambda[0] = (1. - xi) / 2.; 
    N_lambda[1] = (1. + xi) / 2.; 
    return;
}

void evaluate_dN_u2_dxi(const double xi, double* dN_u2_dxi){
    dN_u2_dxi[0] = 3. / 4. * (-1. + xi * xi); 
    dN_u2_dxi[1] = 1. / 4. * (1. + 2. * xi - 3. * xi * xi); 
    dN_u2_dxi[2] = -(3. / 4.) * (-1 + xi * xi);
    dN_u2_dxi[3] = 1. / 4. * (1. - 2. * xi - 3. * xi * xi);
    return;
}

void evaluate_d2N_u2_dxi2(const double xi, double* d2N_u2_dxi2){
    d2N_u2_dxi2[0] = (3 * xi) / 2.; 
    d2N_u2_dxi2[1] = 1. / 2. * (1. - 3. * xi); 
    d2N_u2_dxi2[2] = -((3. * xi) / 2.);
    d2N_u2_dxi2[3] = 1. / 2. * (-1. - 3. * xi);
    return;
}


void computeElementResidual(const double youngs_mod, const double area, const double moment_of_inertia, 
                            const int int_rule, const double* int_points, const double* int_wts, 
                            const double* nodes, const double* unknowns, double* residual){
        
    double J, wtsJ, temp1, temp2, temp3, temp4;
    double dN_u1_dxi[3];
    double dN_u2_dxi[4];
    double d2N_u2_dxi2[4];
    double N_lambda[2];

    //compute Jacobian
    J = (nodes[2] - nodes[0]) / 2.;

    //Loop over integration points
    for(int i = 0; i < int_rule; ++i){

        //compute dNu1dxi
        evaluate_dN_u1_dxi(int_points[i], dN_u1_dxi);
        //compute dNu2dxi
        evaluate_dN_u2_dxi(int_points[i], dN_u2_dxi);
        //compute d2Nu2dxixi
        evaluate_d2N_u2_dxi2(int_points[i], d2N_u2_dxi2);
        //compute N_lambda
        evaluate_N_lambda(int_points[i], N_lambda);

        //multiply Jacobian into integration weights
        wtsJ = int_wts[i] * J;

        //compute dN_u1_dxi . u1
        temp1 = (dN_u1_dxi[0] * unknowns[0] + 
                 dN_u1_dxi[1] * unknowns[4] + 
                 dN_u1_dxi[2] * unknowns[5]) / J;

        //compute dN_u2_dxi . u2
        temp2 = (dN_u2_dxi[0] * unknowns[1] + 
                 dN_u2_dxi[1] * unknowns[2] + 
                 dN_u2_dxi[2] * unknowns[6] +
                 dN_u2_dxi[3] * unknowns[7]) / J;


        //compute d2N_u2_dxi2 . u2
        temp3 = (d2N_u2_dxi2[0] * unknowns[1] + 
                 d2N_u2_dxi2[1] * unknowns[2] + 
                 d2N_u2_dxi2[2] * unknowns[6] +
                 d2N_u2_dxi2[3] * unknowns[7]) / J / J;

        //compute N_lambda . lambda
        temp4 = N_lambda[0] * unknowns[3] + 
                N_lambda[1] * unknowns[8];

        //compute element residual
        residual[0] += dN_u1_dxi[0] / J * (youngs_mod * area * (temp1 - 0.5 * temp2 * temp2) 
                       - 2. * temp4 * (1. + temp1)) * wtsJ;
        residual[4] += dN_u1_dxi[1] / J * (youngs_mod * area * (temp1 - 0.5 * temp2 * temp2) 
                       - 2. * temp4 * (1. + temp1)) * wtsJ;
        residual[5] += dN_u1_dxi[2] / J * (youngs_mod * area * (temp1 - 0.5 * temp2 * temp2) 
                       - 2. * temp4 * (1. + temp1)) * wtsJ;


        residual[1] += (dN_u2_dxi[0] / J * ((youngs_mod * area * (temp1 - 0.5 * temp2 * temp2)) * temp2  
                       - 2. * temp4 * temp2) + d2N_u2_dxi2[0] / J / J * 
                       (youngs_mod * moment_of_inertia * temp3)) * wtsJ;
        residual[2] += (dN_u2_dxi[1] / J * ((youngs_mod * area * (temp1 - 0.5 * temp2 * temp2)) * temp2  
                       - 2. * temp4 * temp2) + d2N_u2_dxi2[1] / J / J * 
                       (youngs_mod * moment_of_inertia * temp3)) * wtsJ;
        residual[6] += (dN_u2_dxi[2] / J * ((youngs_mod * area * (temp1 - 0.5 * temp2 * temp2)) * temp2  
                       - 2. * temp4 * temp2) + d2N_u2_dxi2[2] / J / J * 
                       (youngs_mod * moment_of_inertia * temp3)) * wtsJ;
        residual[7] += (dN_u2_dxi[3] / J * ((youngs_mod * area * (temp1 - 0.5 * temp2 * temp2)) * temp2  
                       - 2. * temp4 * temp2) + d2N_u2_dxi2[3] / J / J * 
                       (youngs_mod * moment_of_inertia * temp3)) * wtsJ;

        residual[3] += N_lambda[0] * ( (1. + temp1) * (1. + temp1) + temp2 * temp2 - 1.);
        residual[8] += N_lambda[1] * ( (1. + temp1) * (1. + temp1) + temp2 * temp2 - 1.);

    }

    return;
}

void __computeResidual(const double youngs_mod, const double area, const double moment_of_inertia, 
                       const int int_rule, const double* int_points, const double* int_wts, 
                       const int num_elements, const double* nodes, const double* unknowns, 
                       double* &residual){


    //Loop over elements for assembly into residual
    for (int i = 0; i < num_elements; ++i, residual+=9, unknowns+=9, nodes+=3) {


        //Move the pointers back for overlap in the assembly
        if ( i != 0 ){
            residual-=4;
            unknowns-=4;
        }
        
        std::cout << "Hello world!" << std::endl;

        //Compute element residual
        computeElementResidual(youngs_mod, area, moment_of_inertia, int_rule, 
                               int_points, int_wts, nodes, unknowns, residual);

    }

    return;

}


extern "C"
{
    extern void computeResidual(const double youngs_mod, const double area, const double moment_of_inertia, 
                       const int int_rule, const double* int_points, const double* int_wts, 
                       const int num_elements, const double* nodes, const double* unknowns, 
                       double* residual){

        return __computeResidual(youngs_mod, area, moment_of_inertia, int_rule, int_points, int_wts, 
                       num_elements, nodes, unknowns, residual);
    }

}
