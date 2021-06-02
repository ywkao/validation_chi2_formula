#include <stdio.h>
#include <iostream>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1D.h>
#include "chi2_helper.h"

int main()
{
    //----------------------------------------------------------------------------------------------------
    // inverse matrix
    //----------------------------------------------------------------------------------------------------
    float m11 = 305.14;
    float m12 = 282.18;
    float m22 = 572.63;
    float determinant = m11*m22 - m12*m12;

    std::cout<< "Inverse matrix: " << std::endl;
    std::cout << "|" << m22/determinant << "  " << -1. * m12/determinant << "|" << std::endl;
    std::cout << "|" << -1. * m12/determinant << "  " << m11/determinant << "|" << std::endl;

    //----------------------------------------------------------------------------------------------------
    // validation
    //----------------------------------------------------------------------------------------------------
    TRandom3 rndm(1234);

    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    TH1D *h = new TH1D("h", "; Difference; Entries", 100, -0.1, 0.1);

    for(int i=0; i<10000; ++i)
    {
        double w_mass_random = 20*rndm.Rndm() + 80.;
        double top_mass_random = 40*rndm.Rndm() + 173.;

        double value_from_original = chi2_calculator_2x2(w_mass_random, top_mass_random);
        double value_from_function = chi2_function(w_mass_random, top_mass_random);
        double difference = (value_from_original - value_from_function) / value_from_function;
        h->Fill(difference);
    }

    h->Draw();
    c1->SaveAs("validation.png");

    return 0;
}
