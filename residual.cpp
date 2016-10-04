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

void evaluate_N_u2(const double xi, const double le, double* N_u2){
    N_u2[0] = 1 / 4. * (1 - xi) * (1 - xi) * (2 + xi);
    N_u2[1] = le / 8. * (1 - xi) * (1 - xi) * (1 + xi);
    N_u2[2] = 1 / 4. * (1 + xi) * (1 + xi) * (2 - xi);
    N_u2[3] = le / 8. * (1 + xi) * (1 + xi) * (xi - 1);
    return;
}


void __generate_u2(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                   const double* solution, const double* nodes, 
                   double* x, double* u2){

    double le, xi, delta_xi;
    double N_u2[4];

    //Evaluate at the last node
    le = nodes[num_elements] - nodes[num_elements - 1];

    evaluate_N_u2(1.0, le, N_u2);

    u2[num_elements * num_of_eval_pts_per_element] = (N_u2[0] * solution[num_dof - 4] + 
                                                      N_u2[1] * solution[num_dof - 3] +
                                                      N_u2[2] * solution[num_dof - 2] +
                                                      N_u2[3] * solution[num_dof - 1]);

    x[num_elements * num_of_eval_pts_per_element] = nodes[num_elements];

    //Compute the increment of \xi for evaluation points in elements
    delta_xi = 2. / num_of_eval_pts_per_element;

    for (int i = 0; i < num_elements; ++i, solution+=2) {

        le = nodes[i + 1] - nodes[i];

        //Evaluate u2 at each evaluation point starting at the left
        //side of the element.
        for (int j = 0; j < num_of_eval_pts_per_element; ++j){
            
            xi = -1.0 + j * delta_xi;

            evaluate_N_u2(xi, le, N_u2);

            u2[i * num_of_eval_pts_per_element + j] = N_u2[0] * solution[0] + 
                                                      N_u2[1] * solution[1] +
                                                      N_u2[2] * solution[2] +
                                                      N_u2[3] * solution[3];

            //x[i * num_of_eval_pts_per_element + j] = N_u2[0] * nodes[i] + 
                                                     //N_u2[2] * nodes[i + 1];
            x[i * num_of_eval_pts_per_element + j] = nodes[i] + le / 2. * (xi + 1.0);
            
        }

    }

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

        
    double le, J, wtsJ, temp3;
    double d2N_u2_dxi2[4];

    //compute element length
    le = nodes[1] - nodes[0];
    //compute Jacobian
    J =  le / 2.;

    //Loop over integration points
    for(int i = 0; i < int_rule; ++i){

        //compute d2Nu2dxixi
        evaluate_d2N_u2_dxi2(int_points[i], le, d2N_u2_dxi2);

        //multiply Jacobian into integration weights
        wtsJ = int_wts[i] * J;
        
        // dofs: 0 -> u_2       - node 1
        //       1 -> \theta    - node 1
        //
        //       2 -> u_2       - node 2
        //       3 -> \theta    - node 2

        //compute d2N_u2_dxi2 . u2
        temp3 = (d2N_u2_dxi2[0] * unknowns[0] + 
                 d2N_u2_dxi2[1] * unknowns[1] + 
                 d2N_u2_dxi2[2] * unknowns[2] +
                 d2N_u2_dxi2[3] * unknowns[3]) / J / J;


        residual[0] += d2N_u2_dxi2[0] / J / J * youngs_mod * moment_of_inertia * temp3 * wtsJ;
        residual[1] += d2N_u2_dxi2[1] / J / J * youngs_mod * moment_of_inertia * temp3 * wtsJ;
        residual[2] += d2N_u2_dxi2[2] / J / J * youngs_mod * moment_of_inertia * temp3 * wtsJ;
        residual[3] += d2N_u2_dxi2[3] / J / J * youngs_mod * moment_of_inertia * temp3 * wtsJ;

    }

    return;
}


void __computeResidual(const double youngs_mod, const double area, const double moment_of_inertia, 
                        const int int_rule, const double* int_points, const double* int_wts, 
                        const int num_elements, const int num_of_nodes_per_element, 
                        const int num_dofs_per_node, const int num_dof, const double* nodes, 
                        const double* unknowns, double* residual){


    //Ensure that residual is 0.0
    memset(residual, 0.0, num_dof*sizeof(double));

    //Loop over elements for assembly into residual
    for (int i = 0; i < num_elements; ++i) {

        //Compute element residual and sum into global
        computeElementResidual2(youngs_mod, area, moment_of_inertia, int_rule, 
                               int_points, int_wts, &nodes[(num_of_nodes_per_element - 1) * i], 
                               &unknowns[(num_dofs_per_node) * i], 
                               &residual[(num_dofs_per_node) * i]);

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
    extern void computeResidual(const double youngs_mod, const double area, const double moment_of_inertia, 
                        const int int_rule, const double* int_points, const double* int_wts, 
                        const int num_elements, const int num_of_nodes_per_element, 
                        const int num_dofs_per_node, const int num_dof, const double* nodes, 
                        const double* unknowns, double* residual){

        return __computeResidual(youngs_mod,  area,  moment_of_inertia, int_rule, int_points, int_wts, 
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
    
    extern void generate_u2(int num_of_eval_pts_per_element, int num_dof, int num_elements, 
                   const double* solution, const double* nodes, 
                   double* x, double* u2){

        return __generate_u2(num_of_eval_pts_per_element, num_dof, num_elements, 
                             solution, nodes, x, u2);
    }

}
