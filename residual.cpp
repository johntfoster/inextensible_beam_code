#include <iostream>
#include <math.h>

template<unsigned p, unsigned k> 
//p = interpolation order, k = derivative of basis function, k = 0 returns basis
//function itself
struct Bspline {
    static void evaluate(const double x, const unsigned int num_knots, const double* knot_vector, double* N) {

        if (k == 0) { // Compute the basis function

            Bspline<p - 1, 0>::evaluate(x, num_knots, knot_vector, N);

            for (int i = 0; i < num_knots - p; ++i) {

                if (fabs(knot_vector[i + p] - knot_vector[i]) < 1e-15)
                    N[i] *= 0.0;
                else
                    N[i] *= (x - knot_vector[i]) / (knot_vector[i + p] - knot_vector[i]);

                if (fabs(knot_vector[i + p + 1] - knot_vector[i + 1]) < 1e-15)
                    N[i] += 0.0;
                else
                    N[i] += (knot_vector[i + p + 1] - x) / 
                            (knot_vector[i + p + 1] - knot_vector[i + 1]) * N[i + 1];
                
            }

            return;
        
        } else { // Compute kth derivative

            Bspline<p - 1, k - 1>::evaluate(x, num_knots, knot_vector, N);

            for (int i = 0; i < num_knots - p; ++i) {

                if (fabs(knot_vector[i + p] - knot_vector[i]) < 1e-15)
                    N[i] *= 0.0;
                else
                    N[i] *= p / (knot_vector[i + p] - knot_vector[i]);

                if (fabs(knot_vector[i + p + 1] - knot_vector[i + 1]) < 1e-15)
                    N[i] += 0.0;
                else
                    N[i] -= p / (knot_vector[i + p + 1] - knot_vector[i + 1]) * N[i + 1];
            }

            return;

        }
    }
};


template<> 
struct Bspline<0, 0>{
  static void evaluate(const double x, const int num_knots, const double* knot_vector, double* N){

        for (int i = 0; i < num_knots; ++i) {

            if (x > knot_vector[i] && x < knot_vector[i + 1])
                N[i] = 1.;
            else
                N[i] = 0.;
        }
    }
};

void evaluate_N_u1_4node(const double xi, double* N_u1){
   N_u1[0] = -0.0625 + xi / 16. + (9 * xi * xi) / 16. - (9* xi * xi * xi) / 16.;
   N_u1[1] =  0.5625 - (27 * xi) / 16. - (9 * xi * xi) / 16. + (27 * xi * xi * xi) / 16.;
   N_u1[2] =  0.5625 + (27 * xi) / 16. - (9 * xi * xi) / 16. - (27 * xi * xi * xi) / 16.;
   N_u1[3] = -0.0625 - xi / 16. + (9 * xi * xi) / 16. + (9* xi * xi * xi) /16.;
   return;
}

void evaluate_dN_u1_dxi_4node(const double xi, double* dN_u1_dxi){
   dN_u1_dxi[0] = 0.0625 + (9 * xi) / 8. - (27 * xi * xi) / 16.;
   dN_u1_dxi[1] = -1.6875 - (9 * xi) / 8. + (81 * xi * xi) / 16.;
   dN_u1_dxi[2] = 1.6875 - (9 * xi) / 8. - (81 * xi * xi) / 16.;
   dN_u1_dxi[3] = -0.0625 + (9 * xi) / 8. + (27 * xi * xi) / 16.;
   return;
}

void evaluate_dN_u1_dxi(const double xi, double* dN_u1_dxi){
    dN_u1_dxi[0] = 1. / 2. * (-1 + xi) + xi / 2.; 
    dN_u1_dxi[1] = -2. * xi; 
    dN_u1_dxi[2] = xi / 2. + (1 + xi) / 2.;
    return;
}

void evaluate_dN_u2_dxi(const double xi, const double le, double* dN_u2_dxi){
    dN_u2_dxi[0] = (1 - xi) * (1 - xi) / 4. - ((1 - xi) * (2 + xi)) / 2.; 
    dN_u2_dxi[1] = (le * (1 - xi) * (1 - xi)) / 8. - (le * (1 - xi) * (1 + xi)) / 4.; 
    dN_u2_dxi[2] = ((2 - xi) * (1 + xi)) / 2. - (1 + xi) * (1 + xi) / 4.;
    dN_u2_dxi[3] = (le *(-1 + xi) * (1 + xi)) / 4. + (le * (1 + xi) * (1 + xi))/ 8.;
    return;
}

void evaluate_d2N_u2_dxi2(const double xi, const double le, double* d2N_u2_dxi2){
    d2N_u2_dxi2[0] = -1 + xi + (2. + xi) / 2.;
    d2N_u2_dxi2[1] = (le * (-4 * (1. - xi) + 2. * (1. + xi))) / 8.; 
    d2N_u2_dxi2[2] = -1 + (2 - xi) / 2. - xi;
    d2N_u2_dxi2[3] = (le * (2 * (-1 + xi) + 4 * (1 + xi))) / 8.;
    return;
}

void evaluate_N_lambda(const double xi, double* N_lambda){
    N_lambda[0] = (1. - xi) / 2.; 
    N_lambda[1] = (1. + xi) / 2.; 
    return;
}

void evaluate_N_u1(const double xi, double* N_u1){
    N_u1[0] = -xi / 2. + xi * xi / 2.;
    N_u1[1] = 1 - xi * xi;
    N_u1[2] = xi / 2. + xi * xi / 2.;
    return;
}

void evaluate_N_u2(const double xi, const double le, double* N_u2){
    N_u2[0] = 1 / 4. * (1 - xi) * (1 - xi) * (2 + xi);
    N_u2[1] = le / 8. * (1 - xi) * (1 - xi) * (1 + xi);
    N_u2[2] = 1 / 4. * (1 + xi) * (1 + xi) * (2 - xi);
    N_u2[3] = le / 8. * (1 + xi) * (1 + xi) * (xi - 1);
    return;
}

void __evaluate_constraint(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                       const double* solution, const double* nodes, 
                       double* x, double* constraint){

    double le, xi, delta_xi, temp1, temp2, J;
    double dN_u1_dxi[3];
    double dN_u2_dxi[4];

    //Evaluate at the last node
    le = nodes[num_elements] - nodes[num_elements - 2];

    J = le / 2.;

    evaluate_dN_u1_dxi(1., dN_u1_dxi);
    evaluate_dN_u2_dxi(1., le, dN_u2_dxi);

    //compute dN_u1_dxi . u1
    temp1 = (dN_u1_dxi[0] * solution[num_dof - 9] + 
             dN_u1_dxi[1] * solution[num_dof - 5] + 
             dN_u1_dxi[2] * solution[num_dof - 4]) / J;

    //compute dN_u2_dxi . u2
    temp2 = (dN_u2_dxi[0] * solution[num_dof - 8] + 
             dN_u2_dxi[1] * solution[num_dof - 7] + 
             dN_u2_dxi[2] * solution[num_dof - 3] +
             dN_u2_dxi[3] * solution[num_dof - 2]) / J;

    constraint[num_elements * num_of_eval_pts_per_element] = (1. + temp1) * (1. + temp1) + temp2 * temp2 - 1.;

    x[num_elements * num_of_eval_pts_per_element] = nodes[num_elements * 2];

    //Compute the increment of \xi for evaluation points in elements
    delta_xi = 2. / num_of_eval_pts_per_element;

    for (int i = 0; i < num_elements; ++i, solution+=5) {

        le = nodes[2 * i + 2] - nodes[2 * i];


        //Evaluate constraint at each evaluation point starting at the left
        //side of the element.
        for (int j = 0; j < num_of_eval_pts_per_element; ++j){
            
            xi = -1.0 + j * delta_xi;

            evaluate_dN_u1_dxi(xi, dN_u1_dxi);
            evaluate_dN_u2_dxi(xi, le, dN_u2_dxi);

            //compute dN_u1_dxi . u1
            temp1 = (dN_u1_dxi[0] * solution[0] + 
                     dN_u1_dxi[1] * solution[4] + 
                     dN_u1_dxi[2] * solution[5]) / J;

            //compute dN_u2_dxi . u2
            temp2 = (dN_u2_dxi[0] * solution[1] + 
                     dN_u2_dxi[1] * solution[2] + 
                     dN_u2_dxi[2] * solution[6] +
                     dN_u2_dxi[3] * solution[7]) / J;

            constraint[i * num_of_eval_pts_per_element + j] = (1. + temp1) * (1. + temp1) + temp2 * temp2 - 1.;

            x[i * num_of_eval_pts_per_element + j] = nodes[2*i] + le / 2. * (xi + 1.0);

        }

    }

    return;
}


void __evaluate_lambda(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                       const double* solution, const double* nodes, 
                       double* x, double* lambda){

    double le, xi, delta_xi;
    double N_lambda[2];

    //Evaluate at the last node
    le = nodes[num_elements] - nodes[num_elements - 2];

    evaluate_N_lambda(1.0, N_lambda);

    lambda[num_elements * num_of_eval_pts_per_element] = N_lambda[0] * solution[num_dof - 6] + 
                                                         N_lambda[1] * solution[num_dof - 1];

    x[num_elements * num_of_eval_pts_per_element] = nodes[num_elements * 2];

    //Compute the increment of \xi for evaluation points in elements
    delta_xi = 2. / num_of_eval_pts_per_element;

    for (int i = 0; i < num_elements; ++i, solution+=5) {

        le = nodes[2 * i + 2] - nodes[2 * i];

        //Evaluate lambda at each evaluation point starting at the left
        //side of the element.
        for (int j = 0; j < num_of_eval_pts_per_element; ++j){
            
            xi = -1.0 + j * delta_xi;

            evaluate_N_lambda(xi, N_lambda);

            lambda[i * num_of_eval_pts_per_element + j] = N_lambda[0] * solution[3] + 
                                                          N_lambda[1] * solution[8];

            x[i * num_of_eval_pts_per_element + j] = N_lambda[0] * nodes[2 * i] +
                                                     N_lambda[1] * nodes[2 * i + 2];
        }

    }

    return;
}


void __evaluate_u1(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                   const double* solution, const double* nodes, 
                   double* x, double* u1){

    double le, xi, delta_xi;
    double N_u1[3];

    //Evaluate at the last node
    le = nodes[num_elements] - nodes[num_elements - 2];

    evaluate_N_u1(1.0, N_u1);

    u1[num_elements * num_of_eval_pts_per_element] = (N_u1[0] * solution[num_dof - 9] + 
                                                      N_u1[1] * solution[num_dof - 5] +
                                                      N_u1[2] * solution[num_dof - 4]);

    x[num_elements * num_of_eval_pts_per_element] = nodes[num_elements * 2];

    //Compute the increment of \xi for evaluation points in elements
    delta_xi = 2. / num_of_eval_pts_per_element;

    for (int i = 0; i < num_elements; ++i, solution+=5) {

        le = nodes[2 * i + 2] - nodes[2 * i];

        //Evaluate u2 at each evaluation point starting at the left
        //side of the element.
        for (int j = 0; j < num_of_eval_pts_per_element; ++j){
            
            xi = -1.0 + j * delta_xi;

            evaluate_N_u1(xi, N_u1);

            u1[i * num_of_eval_pts_per_element + j] = N_u1[0] * solution[0] + 
                                                      N_u1[1] * solution[4] +
                                                      N_u1[2] * solution[5];

            //x[i * num_of_eval_pts_per_element + j] = nodes[i] + le / 2. * (xi + 1.0);
            x[i * num_of_eval_pts_per_element + j] = N_u1[0] * nodes[2*i] +
                                                     N_u1[1] * nodes[2*i + 1] +
                                                     N_u1[2] * nodes[2*i + 2];
        }

    }

    return;
}


void __evaluate_u2(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                   const double* solution, const double* nodes, 
                   double* x, double* u2){

    double le, xi, delta_xi;
    double N_u2[4];

    //Evaluate at the last node
    le = nodes[num_elements] - nodes[num_elements - 2];

    evaluate_N_u2(1.0, le, N_u2);

    u2[num_elements * num_of_eval_pts_per_element] = (N_u2[0] * solution[num_dof - 8] + 
                                                      N_u2[1] * solution[num_dof - 7] +
                                                      N_u2[2] * solution[num_dof - 3] +
                                                      N_u2[3] * solution[num_dof - 2]);

    x[num_elements * num_of_eval_pts_per_element] = nodes[num_elements * 2];

    //Compute the increment of \xi for evaluation points in elements
    delta_xi = 2. / num_of_eval_pts_per_element;

    for (int i = 0; i < num_elements; ++i, solution+=5) {

        le = nodes[2 * i + 2] - nodes[2 * i];

        //Evaluate u2 at each evaluation point starting at the left
        //side of the element.
        for (int j = 0; j < num_of_eval_pts_per_element; ++j){
            
            xi = -1.0 + j * delta_xi;

            evaluate_N_u2(xi, le, N_u2);

            u2[i * num_of_eval_pts_per_element + j] = N_u2[0] * solution[1] + 
                                                      N_u2[1] * solution[2] +
                                                      N_u2[2] * solution[6] +
                                                      N_u2[3] * solution[7];

            x[i * num_of_eval_pts_per_element + j] = nodes[2*i] + le / 2. * (xi + 1.0);
            
        }

    }

    return;
}



void computeElementResidual1(const double youngs_mod, const double area, const double moment_of_inertia, 
                            const int int_rule, const double* int_points, const double* int_wts, 
                            const double* nodes, const double* unknowns, double* residual){
        
    double le, J, wtsJ, temp1, temp2, temp3, temp4;
    double dN_u1_dxi[3];
    double dN_u2_dxi[4];
    double d2N_u2_dxi2[4];
    double N_lambda[2];

    //compute element length
    //compute Jacobian
    le = nodes[2] - nodes[0];
    J = le / 2.;

    //Loop over integration points
    for(int i = 0; i < int_rule; ++i){

        //compute dNu1dxi
        evaluate_dN_u1_dxi(int_points[i], dN_u1_dxi);
        //compute dNu2dxi
        evaluate_dN_u2_dxi(int_points[i], le, dN_u2_dxi);
        //compute d2Nu2dxixi
        evaluate_d2N_u2_dxi2(int_points[i], le, d2N_u2_dxi2);
        //compute N_lambda
        evaluate_N_lambda(int_points[i], N_lambda);

        //multiply Jacobian into integration weights
        wtsJ = int_wts[i] * J;

        // dofs: 0 -> u_1       - node 1
        //       1 -> u_2       - node 1
        //       2 -> \theta    - node 1
        //       3 -> \lambda   - node 1
        //
        //       4 -> u_1       - node 2
        //
        //       5 -> u_1       - node 3
        //       6 -> u_2       - node 3
        //       7 -> \theta    - node 3
        //       8 -> \lambda   - node 3 

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
        //temp4 = 0.0;

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

void computeElementResidual2(const double youngs_mod, const double area, const double moment_of_inertia, 
                            const int int_rule, const double* int_points, const double* int_wts, 
                            const double* nodes, const double* unknowns, double* residual){
        
    double le, J, wtsJ, temp1, temp2, temp3, temp4;
    double dN_u1_dxi[3];
    double dN_u2_dxi[4];
    double d2N_u2_dxi2[4];
    double N_lambda[2];

    //compute element length
    //compute Jacobian
    le = nodes[2] - nodes[0];
    J = le / 2.;

    //Loop over integration points
    for(int i = 0; i < int_rule; ++i){

        //compute dNu1dxi
        evaluate_dN_u1_dxi(int_points[i], dN_u1_dxi);
        //compute dNu2dxi
        evaluate_dN_u2_dxi(int_points[i], le, dN_u2_dxi);
        //compute d2Nu2dxixi
        evaluate_d2N_u2_dxi2(int_points[i], le, d2N_u2_dxi2);
        //compute N_lambda
        evaluate_N_lambda(int_points[i], N_lambda);

        //multiply Jacobian into integration weights
        wtsJ = int_wts[i] * J;

        // dofs: 0 -> u_1       - node 1
        //       1 -> u_2       - node 1
        //       2 -> \theta    - node 1
        //       3 -> \lambda   - node 1
        //
        //       4 -> u_1       - node 2
        //
        //       5 -> u_1       - node 3
        //       6 -> u_2       - node 3
        //       7 -> \theta    - node 3
        //       8 -> \lambda   - node 3 

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
        //temp4 = N_lambda[0] * unknowns[3] + 
                //N_lambda[1] * unknowns[8];
        temp4 = 0.0;

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

        //residual[3] += N_lambda[0] * ( (1. + temp1) * (1. + temp1) + temp2 * temp2 - 1.);
        //residual[8] += N_lambda[1] * ( (1. + temp1) * (1. + temp1) + temp2 * temp2 - 1.);
        residual[3] += 0.0;
        residual[8] += 0.0;

    }

    return;
}

void computeElementResidual3(const double youngs_mod, const double area, const double moment_of_inertia, 
                            const int int_rule, const double* int_points, const double* int_wts, 
                            const double* nodes, const double* unknowns, double* residual){
        
    double le, J, wtsJ, temp1, temp2, temp3, temp4;
    double dN_u1_dxi[3];
    double dN_u2_dxi[4];
    double d2N_u2_dxi2[4];
    double N_lambda[2];

    //compute element length
    //compute Jacobian
    le = nodes[2] - nodes[0];
    J = le / 2.;

    //Loop over integration points
    for(int i = 0; i < int_rule; ++i){

        //compute dNu1dxi
        evaluate_dN_u1_dxi(int_points[i], dN_u1_dxi);
        //compute dNu2dxi
        evaluate_dN_u2_dxi(int_points[i], le, dN_u2_dxi);
        //compute d2Nu2dxixi
        evaluate_d2N_u2_dxi2(int_points[i], le, d2N_u2_dxi2);
        //compute N_lambda
        evaluate_N_lambda(int_points[i], N_lambda);

        //multiply Jacobian into integration weights
        wtsJ = int_wts[i] * J;

        // dofs: 0 -> u_1       - node 1
        //       1 -> u_2       - node 1
        //       2 -> \theta    - node 1
        //       3 -> \lambda   - node 1
        //
        //       4 -> u_1       - node 2
        //
        //       5 -> u_1       - node 3
        //       6 -> u_2       - node 3
        //       7 -> \theta    - node 3
        //       8 -> \lambda   - node 3 

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
        //temp4 = 0.0;

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


void __computeResidual(const int beam_type_flag, const double youngs_mod, const double area, 
                       const double moment_of_inertia, const int int_rule, const double* int_points, 
                       const double* int_wts, const int num_elements, const int num_of_nodes_per_element, 
                       const int num_dofs_per_node, const int num_dof, const double* nodes, 
                       const double* unknowns, double* residual){


    //Ensure that residual is 0.0
    memset(residual, 0.0, num_dof*sizeof(double));

    //Loop over elements for assembly into residual
    for (int i = 0; i < num_elements; ++i) {

        //Compute element residual and sum into global
        if(beam_type_flag == 0){
        computeElementResidual1(youngs_mod, area, moment_of_inertia, int_rule, 
                               int_points, int_wts, &nodes[(num_of_nodes_per_element - 1) * i], 
                               &unknowns[(num_dofs_per_node) * i], 
                               &residual[(num_dofs_per_node) * i]);
        } else if (beam_type_flag == 1){
        computeElementResidual2(youngs_mod, area, moment_of_inertia, int_rule, 
                               int_points, int_wts, &nodes[(num_of_nodes_per_element - 1) * i], 
                               &unknowns[(num_dofs_per_node) * i], 
                               &residual[(num_dofs_per_node) * i]);
        } else if (beam_type_flag == 2){
        computeElementResidual3(youngs_mod, area, moment_of_inertia, int_rule, 
                               int_points, int_wts, &nodes[(num_of_nodes_per_element - 1) * i], 
                               &unknowns[(num_dofs_per_node) * i], 
                               &residual[(num_dofs_per_node) * i]);
        }

    }

    return;

}


//void __computeJacobian(const double youngs_mod, const double area, const double moment_of_inertia, 
                       //const double delta, const int int_rule, const double* int_points, const double* int_wts, 
                       //const int num_elements, const int num_dof,  const double* nodes, 
                       //const double* unknowns, double* jacobian){

  //double f[num_dof];
  //double fp[num_dof];
  //double perturbed_unknowns[num_dof];

  //double* f_ptr = f;
  //double* fp_ptr = fp;

  ////Copy unknowns to perturbed unknowns
  //memcpy(&perturbed_unknowns, &unknowns, num_dof*sizeof(double));

  ////Compute reference residual, returns values to f_ptr
  //__computeResidual2(youngs_mod, area, moment_of_inertia, int_rule, int_points, int_wts, 
                     //num_elements, num_dof, nodes, unknowns, f_ptr);



  ////Loop over all degrees of freedom (rows of jacobian)
  //for (int i = 0; i < num_dof; ++i) {

      ////Perturb a degree of freedom
      //perturbed_unknowns[i] += delta;
      
      ////Compute perturbed residual, returns values to fp_ptr
      //__computeResidual2(youngs_mod, area, moment_of_inertia, int_rule, int_points, int_wts, 
                         //num_elements, num_dof, nodes, perturbed_unknowns, fp_ptr);


      ////Compute forward difference jacobian
      //for (int j = 0; j  < num_dof; ++j){
          //jacobian[i * num_dof + j] = (fp[j] - f[j]) / delta;
      //}

      ////Reset the perturbed dof
      //perturbed_unknowns[i] = unknowns[i];

  //}

  //return;

//}

extern "C"
{
    extern void computeResidual(const int beam_type_flag, const double youngs_mod, const double area, 
                        const double moment_of_inertia, const int int_rule, const double* int_points, 
                        const double* int_wts, const int num_elements, const int num_of_nodes_per_element, 
                        const int num_dofs_per_node, const int num_dof, const double* nodes, 
                        const double* unknowns, double* residual){

        return __computeResidual(beam_type_flag, youngs_mod,  area,  moment_of_inertia, int_rule, int_points, int_wts, 
                         num_elements,  num_of_nodes_per_element, num_dofs_per_node,  num_dof, nodes, 
                         unknowns, residual);
    }

    //extern void computeJacobian(const double youngs_mod, const double area, const double moment_of_inertia, 
                       //const double delta, const int int_rule, const double* int_points, const double* int_wts, 
                       //const int num_elements, const int num_dof,  const double* nodes, 
                       //const double* unknowns, double* jacobian){

        //return __computeJacobian(youngs_mod, area, moment_of_inertia, delta, 
                                 //int_rule,  int_points,  int_wts, num_elements, 
                                 //num_dof, nodes, unknowns, jacobian);
    //}
    
    extern void evaluate_u2(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                       const double* solution, const double* nodes, 
                       double* x, double* u2){

        return __evaluate_u2(num_of_eval_pts_per_element, num_dof, num_elements, 
                       solution, nodes, x, u2);
    }

    extern void evaluate_u1(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                       const double* solution, const double* nodes, 
                       double* x, double* u1){

        return __evaluate_u1(num_of_eval_pts_per_element, num_dof, num_elements, solution, 
                       nodes, x, u1);
    }

    extern void evaluate_lambda(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                       const double* solution, const double* nodes, 
                       double* x, double* lambda){

        return __evaluate_lambda(num_of_eval_pts_per_element, num_dof, num_elements, 
                       solution, nodes, x, lambda);
    }

    extern void evaluate_constraint(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                       const double* solution, const double* nodes, 
                       double* x, double* constraint){

        return __evaluate_constraint(num_of_eval_pts_per_element, num_dof, num_elements, 
                       solution, nodes, x, constraint);
    }

    extern void evaluate_Bspline_p1(const double x, const int num_knots, 
            const double* knot_vector, double* N){

        return Bspline<1, 0>::evaluate(x, num_knots, knot_vector, N);
    }

}
