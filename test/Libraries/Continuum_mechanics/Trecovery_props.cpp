/* This file is part of simcoon.
 
 simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 simcoon is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with simcoon.  If not, see <http://www.gnu.org/licenses/>.
 
 */

///@file Tconstitutive.cpp
///@brief Test for Constitutive tensors in Voigt notation
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/recovery_props.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Trecovery_props, test_check_symetries_iso)
{
    double E = 70000.;
    double nu = 0.3;
    
    int axis = 0;
    int sym_check = 0;
    std::string umat_type;
    vec props;
    
    //Test of check_symetries with L_iso
    mat Ltest = L_iso(E, nu, "Enu");
    check_symetries(Ltest, umat_type, axis, props, sym_check);
    
    EXPECT_EQ(umat_type,"ELISO");
    EXPECT_EQ(axis,0);
    EXPECT_EQ(sym_check,1);
    EXPECT_LT(fabs(E - props(0)),1.E-9);
    EXPECT_LT(fabs(nu - props(1)),1.E-9);
}

TEST(test_check_symetries_cubic, test_check_symetries_cubic_functions)
{
    double E = 1000;
    double nu = 0.2;
    double G = 500;
    
    int axis = 0;
    int sym_check = 0;
    std::string umat_type;
    vec props;
    
    //Test of check_symetries with L_cubic
    mat Ltest = L_cubic(E,nu,G,"EnuG");
    check_symetries(Ltest, umat_type, axis, props, sym_check);
    
    EXPECT_EQ(umat_type,"ELCUB");
    EXPECT_EQ(axis,0);
    EXPECT_EQ(sym_check,1);
    EXPECT_LT(fabs(E - props(0)),1.E-9);
    EXPECT_LT(fabs(nu - props(1)),1.E-9);
    EXPECT_LT(fabs(G - props(2)),1.E-9);
}

TEST(Trecovery_props, test_check_symetries_isotrans)
{
    double EL = 10000.;
    double ET = 20000.;
    double nuTL = 0.3;
    double nuTT = 0.3;
    double GLT = 8000.;
    
    int axis = 0;
    int sym_check = 0;
    std::string umat_type;
    vec props;
    
    //Test of check_symetries with L_isotrans
    mat Ltest = L_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
    check_symetries(Ltest, umat_type, axis, props, sym_check);
    
    cout << "props ELIST1 = " << props.t() << "\n";
    
    EXPECT_EQ(umat_type,"ELIST");
    EXPECT_EQ(axis,1);
    EXPECT_EQ(sym_check,1);
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);
    
    Ltest = L_isotrans(EL, ET, nuTL, nuTT, GLT, 2);
    check_symetries(Ltest, umat_type, axis, props, sym_check);
    
    cout << "props ELIST2 = " << props.t() << "\n";
    
    EXPECT_EQ(umat_type,"ELIST");
    EXPECT_EQ(axis,2);
    EXPECT_EQ(sym_check,1);
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);
    
    Ltest = L_isotrans(EL, ET, nuTL, nuTT, GLT, 3);
    check_symetries(Ltest, umat_type, axis, props, sym_check);
    
    cout << "props ELIST3 = " << props.t() << "\n";
    
    EXPECT_EQ(umat_type,"ELIST");
    EXPECT_EQ(axis,3);
    EXPECT_EQ(sym_check,1);
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);
}

TEST(Trecovery_props, test_check_symetries_ortho)
{
    double E1 = 10000.;
    double E2 = 20000.;
    double E3 = 30000.;
    double nu12 = 0.3;
    double nu13 = 0.3;
    double nu23 = 0.3;
    double G12 = 6000.;
    double G13 = 7000.;
    double G23 = 8000.;
    
    int axis = 0;
    int sym_check = 0;
    std::string umat_type;
    vec props;
    
    //Test of check_symetries with L_cubic
    mat Ltest = L_ortho(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, "EnuG");
    check_symetries(Ltest, umat_type, axis, props, sym_check);
    
    cout << "props ELORT = " << props.t() << "\n";
    
    EXPECT_EQ(umat_type,"ELORT");
    EXPECT_EQ(axis,0);
    EXPECT_EQ(sym_check,1);
    EXPECT_LT(fabs(E1 - props(0)),1.E-9);
    EXPECT_LT(fabs(E2 - props(1)),1.E-9);
    EXPECT_LT(fabs(E3 - props(2)),1.E-9);
    EXPECT_LT(fabs(nu12 - props(3)),1.E-9);
    EXPECT_LT(fabs(nu13 - props(4)),1.E-9);
    EXPECT_LT(fabs(nu23 - props(5)),1.E-9);
    EXPECT_LT(fabs(G12 - props(6)),1.E-9);
    EXPECT_LT(fabs(G13 - props(7)),1.E-9);
    EXPECT_LT(fabs(G23 - props(8)),1.E-9);
}

TEST(Trecovery_props, test_Liso_Miso_props)
{
    double E = 70000.;
    double nu = 0.3;
    
    mat Ltest = L_iso(E, nu, "Enu");
    mat Mtest = M_iso(E, nu, "Enu");

    vec props = L_iso_props(Ltest);
    
    EXPECT_LT(fabs(E - props(0)),1.E-9);
    EXPECT_LT(fabs(nu - props(1)),1.E-9);
    
    props = M_iso_props(Mtest);

    EXPECT_LT(fabs(E - props(0)),1.E-9);
    EXPECT_LT(fabs(nu - props(1)),1.E-9);
}

TEST(Trecovery_props, test_L_cubic_M_cubic_props)
{
    double C11 = 1000;
    double C12 = 400;
    double C44 = 500;
    
    double nu = 1 / (1 + C11/C12);
    double E = C11*(1 - 3*pow(nu,2) - 2*pow(nu,3))/(1-pow(nu,2));
    double G = C44;
    
    //Test of L_cubic function
    mat Ltest = L_cubic(C11,C12,C44, "Cii");
    mat Mtest = M_cubic(C11,C12,C44, "Cii");
    
    vec props = L_cubic_props(Ltest);
    
    EXPECT_LT(fabs(E - props(0)),1.E-9);
    EXPECT_LT(fabs(nu - props(1)),1.E-9);
    EXPECT_LT(fabs(G - props(2)),1.E-9);
    
    props = M_cubic_props(Mtest);
    
    EXPECT_LT(fabs(E - props(0)),1.E-9);
    EXPECT_LT(fabs(nu - props(1)),1.E-9);
    EXPECT_LT(fabs(G - props(2)),1.E-9);
    
}

TEST(Trecovery_props, test_L_isoT_M_isoT_props)
{
    double EL = 10000.;
    double ET = 20000.;
    double nuTL = 0.3;
    double nuTT = 0.3;
    double GLT = 8000.;
    
    //Test of L_cubic function
    mat Ltest = L_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
    mat Mtest = M_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
    
    vec props = L_isotrans_props(Ltest,1);
    
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);
    
    props = M_isotrans_props(Mtest,1);
    
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);

    Ltest = L_isotrans(EL, ET, nuTL, nuTT, GLT, 2);
    Mtest = M_isotrans(EL, ET, nuTL, nuTT, GLT, 2);
    
    props = L_isotrans_props(Ltest,2);
    
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);
    
    props = M_isotrans_props(Mtest,2);
    
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);
    
    Ltest = L_isotrans(EL, ET, nuTL, nuTT, GLT, 3);
    Mtest = M_isotrans(EL, ET, nuTL, nuTT, GLT, 3);
    
    props = L_isotrans_props(Ltest,3);
    
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);
    
    props = M_isotrans_props(Mtest,3);
    
    EXPECT_LT(fabs(EL - props(0)),1.E-9);
    EXPECT_LT(fabs(ET - props(1)),1.E-9);
    EXPECT_LT(fabs(nuTL - props(2)),1.E-9);
    EXPECT_LT(fabs(nuTT - props(3)),1.E-9);
    EXPECT_LT(fabs(GLT - props(4)),1.E-9);

}

