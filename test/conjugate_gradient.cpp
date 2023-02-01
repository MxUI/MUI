#include <iostream>
#include "../linalg/conjugate_gradient.h"
#include "../linalg/ilu_preconditioner.h"
#include "../linalg/ic_preconditioner.h"
#include "../linalg/ssor_preconditioner.h"

int main() {
    std::cout << "===========================================================================" << std::endl;
    std::cout << "==================== TEST 00: 2-by-2 A_ Matrix ============================" << std::endl;
    std::cout << "====MUI Conjugate Gradient Solver w/o Preconditioner & w/o Initial Guess===" << std::endl;
    std::cout << "===========================================================================" << std::endl;

    std::pair<int, double> cgReturn0;

    int Css_size0 = 2;
    int Aas_size0 = 2;
    double cg_solve_tol0 = 1e-6;
    int cg_max_iter0 = 1000;

    mui::linalg::sparse_matrix<int,double> Css0; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas0; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_0; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref0; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff0; //< Difference between reference value and calculated value of Transformation Matrix

    Css0.resize_null(Css_size0, Css_size0);
    Aas0.resize_null(Aas_size0, 1);
    H_0.resize_null(Aas_size0, 1);
    H_ref0.resize_null(Aas_size0, 1);
    H_diff0.resize_null(Aas_size0, 1);

    H_ref0.set_value(0, 0, 0.5488);
    H_ref0.set_value(1, 0, 0.7152);

    Css0.set_value(0, 0, 2.5409);
    Css0.set_value(0, 1, -0.0113);
    Css0.set_value(1, 0, -0.0113);
    Css0.set_value(1, 1, 0.5287);

    Aas0.set_value(0, 0, 1.3864);
    Aas0.set_value(1, 0, 0.3719);

    mui::linalg::conjugate_gradient<int,double> cg0(Css0, Aas0, cg_solve_tol0, cg_max_iter0);
    cgReturn0 = cg0.solve();
    H_0 = cg0.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_0.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref0.print();

    H_diff0 = H_0 - H_ref0;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff0.print();

    std::cout << "Total CG iteration number: " << cgReturn0.first <<" with final r_norm_rel: " << cgReturn0.second << std::endl;

    std::cout << std::endl;

    std::cout << "===========================================================================" << std::endl;
    std::cout << "==================== TEST 01: 2-by-2 A_ Matrix ============================" << std::endl;
    std::cout << "=MUI Conjugate Gradient Solver with ILU Preconditioner & w/o Initial Guess=" << std::endl;
    std::cout << "===========================================================================" << std::endl;

    std::pair<int, double> cgReturn1;
    
    int Css_size1 = 2;
    int Aas_size1 = 2;
    double cg_solve_tol1 = 1e-6;
    int cg_max_iter1 = 1000;

    mui::linalg::sparse_matrix<int,double> Css1; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas1; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_1; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref1; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff1; //< Difference between reference value and calculated value of Transformation Matrix

    Css1.resize_null(Css_size1, Css_size1);
    Aas1.resize_null(Aas_size1, 1);
    H_1.resize_null(Aas_size1, 1);
    H_ref1.resize_null(Aas_size1, 1);
    H_diff1.resize_null(Aas_size1, 1);

    H_ref1.set_value(0, 0, 0.5488);
    H_ref1.set_value(1, 0, 0.7152);
    
    Css1.set_value(0, 0, 2.5409);
    Css1.set_value(0, 1, -0.0113);
    Css1.set_value(1, 0, -0.0113);
    Css1.set_value(1, 1, 0.5287);

    Aas1.set_value(0, 0, 1.3864);
    Aas1.set_value(1, 0, 0.3719);

    mui::linalg::ilu_preconditioner<int,double> M1(Css1);
    mui::linalg::conjugate_gradient<int,double> cg1(Css1, Aas1, cg_solve_tol1, cg_max_iter1, &M1);
    cgReturn1 = cg1.solve();
    H_1 = cg1.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_1.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref1.print();

    H_diff1 = H_1 - H_ref1;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff1.print();

    std::cout << "Total CG iteration number: " << cgReturn1.first <<" with final r_norm_rel: " << cgReturn1.second << std::endl;

    std::cout << std::endl;

    std::cout << "===========================================================================" << std::endl;
    std::cout << "==================== TEST 02: 3-by-3 A_ Matrix ============================" << std::endl;
    std::cout << "=MUI Conjugate Gradient Solver with SSOR Preconditioner & w/o Initial Guess" << std::endl;
    std::cout << "===========================================================================" << std::endl;

    std::pair<int, double> cgReturn2;

    int Css_size2 = 3;
    int Aas_size2 = 3;
    double cg_solve_tol2 = 1e-6;
    int cg_max_iter2 = 1000;

    mui::linalg::sparse_matrix<int,double> Css2; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas2; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_2; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref2; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff2; //< Difference between reference value and calculated value of Transformation Matrix

    Css2.resize_null(Css_size2, Css_size2);
    Aas2.resize_null(Aas_size2, 1);
    H_2.resize_null(Aas_size2, 1);
    H_ref2.resize_null(Aas_size2, 1);
    H_diff2.resize_null(Aas_size2, 1);

    H_ref2.set_value(0, 0, 0.5488);
    H_ref2.set_value(1, 0, 0.7152);
    H_ref2.set_value(2, 0, 0.6028);

    Css2.set_value(0, 0, 0.7444);
    Css2.set_value(0, 1, -0.5055);
    Css2.set_value(0, 2, -0.0851);
    Css2.set_value(1, 0, -0.5055);
    Css2.set_value(1, 1, 3.4858);
    Css2.set_value(1, 2, 0.0572);
    Css2.set_value(2, 0, -0.0851);
    Css2.set_value(2, 1, 0.0572);
    Css2.set_value(2, 2, 0.4738);

    Aas2.set_value(0, 0, -0.0043);
    Aas2.set_value(1, 0, 2.2501);
    Aas2.set_value(2, 0, 0.2798);

    mui::linalg::symmetric_successive_over_relaxation_preconditioner<int,double> M2(Css2, 1.2);
    mui::linalg::conjugate_gradient<int,double> cg2(Css2, Aas2, cg_solve_tol2, cg_max_iter2, &M2);
    cgReturn2 = cg2.solve();
    H_2 = cg2.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_2.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref2.print();

    H_diff2 = H_2 - H_ref2;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff2.print();

    std::cout << "Total CG iteration number: " << cgReturn2.first <<" with final r_norm_rel: " << cgReturn2.second << std::endl;

    std::cout << std::endl;

    std::cout << "===========================================================================" << std::endl;
    std::cout << "==================== TEST 03: 6-by-6 A_ Matrix ============================" << std::endl;
    std::cout << "=MUI Conjugate Gradient Solver with IC Preconditioner & with Initial Guess=" << std::endl;
    std::cout << "===========================================================================" << std::endl;

    std::pair<int, double> cgReturn3;

    int Css_size3 = 6;
    int Aas_size3 = 6;
    double cg_solve_tol3 = 1e-6;
    int cg_max_iter3 = 1000;

    mui::linalg::sparse_matrix<int,double> Css3; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas3; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_init3; //< Initial Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_3; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref3; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff3; //< Difference between reference value and calculated value of Transformation Matrix

    Css3.resize_null(Css_size3, Css_size3);
    Aas3.resize_null(Aas_size3, 1);
    H_init3.resize_null(Aas_size3, 1);
    H_3.resize_null(Aas_size3, 1);
    H_ref3.resize_null(Aas_size3, 1);
    H_diff3.resize_null(Aas_size3, 1);

    H_init3.set_value(0, 0, 9);
    H_init3.set_value(1, 0, 0);
    H_init3.set_value(2, 0, -2);
    H_init3.set_value(3, 0, 3);
    H_init3.set_value(4, 0, -2);
    H_init3.set_value(5, 0, 5);

    H_ref3.set_value(0, 0, 0.5488);
    H_ref3.set_value(1, 0, 0.7152);
    H_ref3.set_value(2, 0, 0.6028);
    H_ref3.set_value(3, 0, 0.5449);
    H_ref3.set_value(4, 0, 0.4237);
    H_ref3.set_value(5, 0, 0.6459);

    Css3.set_value(0, 0, 3.4430);
    Css3.set_value(0, 1, -0.3963);
    Css3.set_value(0, 2, 2.5012);
    Css3.set_value(0, 3, 0.9525);
    Css3.set_value(0, 4, 0.6084);
    Css3.set_value(0, 5, -1.2728);

    Css3.set_value(1, 0, -0.3963);
    Css3.set_value(1, 1, 0.6015);
    Css3.set_value(1, 2, -0.4108);
    Css3.set_value(1, 3, -0.1359);
    Css3.set_value(1, 4, -0.0295);
    Css3.set_value(1, 5, 0.2630);

    Css3.set_value(2, 0, 2.5012);
    Css3.set_value(2, 1, -0.4108);
    Css3.set_value(2, 2, 2.5927);
    Css3.set_value(2, 3, 0.7072);
    Css3.set_value(2, 4, 0.5587);
    Css3.set_value(2, 5, -1.0613);

    Css3.set_value(3, 0, 0.9525);
    Css3.set_value(3, 1, -0.1359);
    Css3.set_value(3, 2, 0.7072);
    Css3.set_value(3, 3, 1.1634);
    Css3.set_value(3, 4, 0.1920);
    Css3.set_value(3, 5, -0.4344);

    Css3.set_value(4, 0, 0.6084);
    Css3.set_value(4, 1, -0.0295);
    Css3.set_value(4, 2, 0.5587);
    Css3.set_value(4, 3, 0.1920);
    Css3.set_value(4, 4, 0.7636);
    Css3.set_value(4, 5, -0.3261);

    Css3.set_value(5, 0, -1.2728);
    Css3.set_value(5, 1, 0.2630);
    Css3.set_value(5, 2, -1.0613);
    Css3.set_value(5, 3, -0.4344);
    Css3.set_value(5, 4, -0.3261);
    Css3.set_value(5, 5, 1.0869);

    Aas3.set_value(0, 0, 3.0685);
    Aas3.set_value(1, 0, 0.0484);
    Aas3.set_value(2, 0, 2.5783);
    Aas3.set_value(3, 0, 1.2865);
    Aas3.set_value(4, 0, 0.8671);
    Aas3.set_value(5, 0, -0.8230);

    mui::linalg::incomplete_cholesky_preconditioner<int,double> M3(Css3);
    mui::linalg::conjugate_gradient<int,double> cg3(Css3, Aas3, cg_solve_tol3, cg_max_iter3, &M3);
    cgReturn3 = cg3.solve(H_init3);
    H_3 = cg3.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_3.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref3.print();

    H_diff3 = H_3 - H_ref3;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff3.print();

    std::cout << "Total CG iteration number: " << cgReturn3.first <<" with final r_norm_rel: " << cgReturn3.second << std::endl;

    std::cout << std::endl;

    return 0;
}