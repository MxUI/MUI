#include <iostream>
#include "../linalg/conjugate_gradient.h"
#include "../linalg/ilu_preconditioner.h"
#include "../linalg/ic_preconditioner.h"

int main() {

	std::pair<int, double> cgReturn;
    
    int Css_size = 2;
    int Aas_size = 2;
    double cg_solve_tol = 1e-6;
	int cg_max_iter = 1000;

	mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
	mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
	mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
	mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
	mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize_null(Css_size, Css_size);
    Aas.resize_null(Aas_size, 1);
    H_.resize_null(Aas_size, 1);
    H_ref.resize_null(Aas_size, 1);
    H_diff.resize_null(Aas_size, 1);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);
    
    Css.set_value(0, 0, 2.5409);
    Css.set_value(0, 1, -0.0113);
    Css.set_value(1, 0, -0.0113);
    Css.set_value(1, 1, 0.5287);

    Aas.set_value(0, 0, 1.3864);
    Aas.set_value(1, 0, 0.3719);


    mui::linalg::incomplete_cholesky_preconditioner<int,double> M(Css);
    mui::linalg::conjugate_gradient<int,double> cg(Css, Aas, cg_solve_tol, cg_max_iter, &M);
    cgReturn = cg.solve();
    H_ = cg.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << "Total CG iteration number: " << cgReturn.first <<" with final r_norm_rel: " << cgReturn.second << std::endl;

    return 0;
}
