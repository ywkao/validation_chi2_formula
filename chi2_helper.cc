#include "chi2_helper.h"

double chi2_function(double x, double y)
{
    double x0 = 85.70;
    double y0 = 174.81;
    double m00 = 305.14;
    double m01 = 282.18;
    double m10 = 282.18;
    double m11 = 572.63;

    double determinant = m00 * m11 - m01 * m10;
    double m00_inverse = m11 / determinant;
    double m01_inverse = -1. * m01 / determinant;
    double m10_inverse = -1. * m10 / determinant;
    double m11_inverse = m00 / determinant;

    double dx = x - x0;
    double dy = y - y0;

    double chi2_value = m00_inverse * dx * dx + (m01_inverse + m10_inverse) * dx * dy + m11_inverse * dy * dy;
    return chi2_value;
}

double chi2_calculator_2x2(double w_mass, double t_mass)
{
    // means & covariance matrix elements are taken from covMatrx_Era2017_M1000.json
    TVectorD vec_mean_values(2);
    
    vec_mean_values(0) = 85.70;
    vec_mean_values(1) = 174.81;

    TVectorD vec_mass(2);
    vec_mass(0) = w_mass - vec_mean_values(0);
    vec_mass(1) = t_mass - vec_mean_values(1);

    TMatrixD matrix_mass(2,1);
    matrix_mass(0,0) = w_mass - vec_mean_values(0);
    matrix_mass(1,0) = t_mass - vec_mean_values(1);

    TMatrixD matrix_massT(1,2);
    matrix_massT(0,0) = w_mass - vec_mean_values(0);
    matrix_massT(0,1) = t_mass - vec_mean_values(1);

    TMatrixD matrix(2,2);
    matrix(0,0) = 305.14;
    matrix(0,1) = 282.18;
    matrix(1,0) = 282.18;
    matrix(1,1) = 572.63;

    double chi2_value = matrix.Invert()*vec_mass*vec_mass;

    //TMatrixD C(1,1);
    //TMatrixD D(1,2);
    //D.Mult(matrix_massT, matrix.Invert());
    //C.Mult(D, matrix_mass);
    //double chi2_value_test      = C(0,0);
    //cout<<"chi2_value       ="<<chi2_value<<endl;
    //cout<<"chi2_value_test  ="<<chi2_value_test<<endl;

    return chi2_value;
};
