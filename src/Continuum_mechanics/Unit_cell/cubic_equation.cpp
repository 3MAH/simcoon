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

///@file cubic_equation.cpp
///@brief Class that gathers periodic equations between nodes in a cubic Representative Volume Element
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <CGAL/squared_distance_3.h> //for 3D functions
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/component.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for cubic_equation===================================

//=====Public methods for cubic_equation============================================

/*!
 \brief default constructor
 */

std::vector<std::string> cubic_equation::set_name_CD = {"CD11", "CD22", "CD33", "CD12", "CD13", "CD23"};
    
//-------------------------------------------------------------
cubic_equation::cubic_equation()
//-------------------------------------------------------------
{

}
    
/*!
 \brief Constructor
 \param cubic_mesh : a cubic mesh object
 \n\n
 */

//-------------------------------------------------------------
cubic_equation::cubic_equation(const cubic_mesh &cm, const cubic_mesh &cm_perio, const unsigned int &loading_type)
//-------------------------------------------------------------
{
    construct(cm,cm_perio,loading_type);
}

/*!
 \brief Copy constructor
 \param cm cubic_equation object to duplicate
 */
    
//------------------------------------------------------
cubic_equation::cubic_equation(const cubic_equation& eq)
//------------------------------------------------------
{
    Corner_listXmYpZm = eq.Corner_listXmYpZm;
    Corner_listXpYmZm = eq.Corner_listXpYmZm;
    Corner_listXpYpZm = eq.Corner_listXpYpZm;
    Corner_listXmYmZp = eq.Corner_listXmYmZp;
    Corner_listXmYpZp = eq.Corner_listXmYpZp;
    Corner_listXpYmZp = eq.Corner_listXpYmZp;
    Corner_listXpYpZp = eq.Corner_listXpYpZp;
    
    Edge_listXpYm = eq.Edge_listXpYm;
    Edge_listXpYp = eq.Edge_listXpYp;
    Edge_listXmYp = eq.Edge_listXmYp;
    Edge_listXpZm = eq.Edge_listXpZm;
    Edge_listXpZp = eq.Edge_listXpZp;
    Edge_listXmZp = eq.Edge_listXmZp;
    Edge_listYpZm = eq.Edge_listYpZm;
    Edge_listYpZp = eq.Edge_listYpZp;
    Edge_listYmZp = eq.Edge_listYmZp;
    
    Face_listXp = eq.Face_listXp;
    Face_listYp = eq.Face_listYp;
    Face_listZp = eq.Face_listZp;
}

/*!
 \brief Destructor
 
 Deletes cubic_equation
 */

//-------------------------------------
cubic_equation::~cubic_equation() {}
//-------------------------------------
    
//-------------------------------------------------------------
void cubic_equation::construct(const cubic_mesh &cm, const cubic_mesh &cm_perio, const unsigned int &loading_type)
//-------------------------------------------------------------
{
    std::vector<int> CD;
    std::vector<int> list_dofs;
    mat Dxyz;
    Mat<int> CD_num;
    //Mechanical
    if (loading_type == 0) {
        CD_num = {{0,3,4},{3,1,5},{4,5,2}};
        Dxyz = {{-1.*cm.Dx,-0.5*cm.Dx,-0.5*cm.Dx},{-0.5*cm.Dy,-1.*cm.Dy,-0.5*cm.Dy},{-0.5*cm.Dz,-0.5*cm.Dz,-1.*cm.Dz}};
        CD = {1010011, 1020022, 1030033, 1040012, 1050013, 1060023};
        list_dofs = {1,2,3};
    }
    //Thermomechanical
    else if (loading_type == 1) {
        CD_num = {{0,3,4,6},{3,1,5,7},{4,5,2,8}};
        Dxyz = {{-1.*cm.Dx,-0.5*cm.Dx,-0.5*cm.Dx},{-0.5*cm.Dy,-1.*cm.Dy,-0.5*cm.Dy},{-0.5*cm.Dz,-0.5*cm.Dz,-1.*cm.Dz},{-1.0*cm.Dx,-1.0*cm.Dy,-1.0*cm.Dz}};
        CD = {1010011, 1020022, 1030033, 1040012, 1050013, 1060023, 1010111, 1020211, 1030311};
        list_dofs = {1,2,3,11};
    }
    //Thermal
    else {
        CD_num.resize(1,3);
        CD_num(0,1) = 6;
        CD_num(0,1) = 7;
        CD_num(0,2) = 8;
        Dxyz.resize(3,1);
        CD_num(0,0) = -1.0*cm.Dx;
        CD_num(1,0) = -1.0*cm.Dy;
        CD_num(2,0) = -1.0*cm.Dz;
        CD = {1010111, 1020211, 1030311};
        list_dofs = {11};
    }
    
    //construct the constrain drivers
    Point null(0.,0.,0.);
    for (auto n:CD) {
        Node node;
        node.number = n;
        node.coords = null;
        Constrain_Driver_123.push_back(node);
    }
    
    unsigned int size_list_dofs = list_dofs.size();
    Corner_listXmYpZm.resize(size_list_dofs);
    Corner_listXpYmZm.resize(size_list_dofs);
    Corner_listXpYpZm.resize(size_list_dofs);
    Corner_listXmYmZp.resize(size_list_dofs);
    Corner_listXmYpZp.resize(size_list_dofs);
    Corner_listXpYmZp.resize(size_list_dofs);
    Corner_listXpYpZp.resize(size_list_dofs);
    
    Edge_listXpYm.resize(size_list_dofs);
    Edge_listXpYp.resize(size_list_dofs);
    Edge_listXmYp.resize(size_list_dofs);
    Edge_listXpZm.resize(size_list_dofs);
    Edge_listXpZp.resize(size_list_dofs);
    Edge_listXmZp.resize(size_list_dofs);
    Edge_listYpZm.resize(size_list_dofs);
    Edge_listYpZp.resize(size_list_dofs);
    Edge_listYmZp.resize(size_list_dofs);
    Face_listXp.resize(size_list_dofs);
    Face_listYp.resize(size_list_dofs);
    Face_listZp.resize(size_list_dofs);
    
    unsigned int p1=0;
    unsigned int p2=0;
    unsigned int p3=0;
    component temp;
    
    for (unsigned int i=0; i<list_dofs.size(); i++) {
        //Corner_listXpYmZm
        p1 = 0;
        temp.node = *cm.Corner_listXpYmZm;
        temp.dof = list_dofs[i];
        temp.coef = 1.0;
        Corner_listXpYmZm[i].components.push_back(temp);
        
        temp.node = *cm.Corner_listXmYmZm;
        temp.dof = list_dofs[i];
        temp.coef = -1.0;
        Corner_listXpYmZm[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p1,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p1,i);
        Corner_listXpYmZm[i].components.push_back(temp);
        
        //Corner_listXmYpZm
        p1 = 1;
        temp.node = *cm.Corner_listXmYpZm;
        temp.dof = list_dofs[i];
        temp.coef = 1.0;
        Corner_listXmYpZm[i].components.push_back(temp);
        
        temp.node = *cm.Corner_listXmYmZm;
        temp.dof = list_dofs[i];
        temp.coef = -1.0;
        Corner_listXmYpZm[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p1,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p1,i);
        Corner_listXmYpZm[i].components.push_back(temp);
        
        //Corner_listXmYmZp
        p1 = 2;
        temp.node = *cm.Corner_listXmYmZp;
        temp.dof = list_dofs[i];
        temp.coef = 1.0;
        Corner_listXmYmZp[i].components.push_back(temp);
        
        temp.node = *cm.Corner_listXmYmZm;
        temp.dof = list_dofs[i];
        temp.coef = -1.0;
        Corner_listXmYmZp[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p1,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p1,i);
        Corner_listXmYmZp[i].components.push_back(temp);
    
        //Corner_listXpYpZm
        p1 = 0;
        p2 = 1;
        temp.node = *cm.Corner_listXpYpZm;
        temp.dof = list_dofs[i];
        temp.coef = 1.0;
        Corner_listXpYpZm[i].components.push_back(temp);
        
        temp.node = *cm.Corner_listXmYmZm;
        temp.dof = list_dofs[i];
        temp.coef = -1.0;
        Corner_listXpYpZm[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p1,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p1,i);
        Corner_listXpYpZm[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p2,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p2,i);
        Corner_listXpYpZm[i].components.push_back(temp);
        
        //Corner_listXpYmZp
        p1 = 0;
        p2 = 2;
        temp.node = *cm.Corner_listXpYmZp;
        temp.dof = list_dofs[i];
        temp.coef = 1.0;
        Corner_listXpYmZp[i].components.push_back(temp);
        
        temp.node = *cm.Corner_listXmYmZm;
        temp.dof = list_dofs[i];
        temp.coef = -1.0;
        Corner_listXpYmZp[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p1,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p1,i);
        Corner_listXpYmZp[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p2,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p2,i);
        Corner_listXpYmZp[i].components.push_back(temp);
    
        //Corner_listXmYpZp
        p1 = 1;
        p2 = 2;
        temp.node = *cm.Corner_listXmYpZp;
        temp.dof = list_dofs[i];
        temp.coef = 1.0;
        Corner_listXmYpZp[i].components.push_back(temp);
        
        temp.node = *cm.Corner_listXmYmZm;
        temp.dof = list_dofs[i];
        temp.coef = -1.0;
        Corner_listXmYpZp[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p1,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p1,i);
        Corner_listXmYpZp[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p2,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p2,i);
        Corner_listXmYpZp[i].components.push_back(temp);
    
        //Corner_listXpYpZp
        p1 = 0;
        p2 = 1;
        p3 = 2;
        temp.node = *cm.Corner_listXpYpZp;
        temp.dof = list_dofs[i];
        temp.coef = 1.0;
        Corner_listXpYpZp[i].components.push_back(temp);
        
        temp.node = *cm.Corner_listXmYmZm;
        temp.dof = list_dofs[i];
        temp.coef = -1.0;
        Corner_listXpYpZp[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p1,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p1,i);
        Corner_listXpYpZp[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p2,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p2,i);
        Corner_listXpYpZp[i].components.push_back(temp);
        
        temp.node = Constrain_Driver_123[CD_num(p3,i)];
        temp.dof = 1;
        temp.coef = Dxyz(p3,i);
        Corner_listXpYpZp[i].components.push_back(temp);
        
        //Edge_listXpYm
        Edge_listXpYm[i].resize(cm_perio.Edge_listXpYm->size());
        for (unsigned int j=0; j<cm_perio.Edge_listXmYm->size(); j++) {
            p1 = 0;
            temp.node = (*cm_perio.Edge_listXpYm)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listXpYm[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listXmYm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listXpYm[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listXpYm[i][j].components.push_back(temp);
        }

        //Edge_listXpYp
        Edge_listXpYp[i].resize(cm_perio.Edge_listXpYp->size());
        for (unsigned int j=0; j<cm_perio.Edge_listXmYm->size(); j++) {
            p1 = 0;
            p2 = 1;
            temp.node = (*cm_perio.Edge_listXpYp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listXpYp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listXmYm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listXpYp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listXpYp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p2,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p2,i);
            Edge_listXpYp[i][j].components.push_back(temp);
        }

        //Edge_listXmYp
        Edge_listXmYp[i].resize(cm_perio.Edge_listXmYp->size());
        for (unsigned int j=0; j<cm_perio.Edge_listXmYm->size(); j++) {
            p1 = 1;
            temp.node = (*cm_perio.Edge_listXmYp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listXmYp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listXmYm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listXmYp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listXmYp[i][j].components.push_back(temp);
        }

        //Edge_listXpZm
        Edge_listXpZm[i].resize(cm_perio.Edge_listXpZm->size());
        for (unsigned int j=0; j<cm_perio.Edge_listXmZm->size(); j++) {
            p1 = 0;
            temp.node = (*cm_perio.Edge_listXpZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listXpZm[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listXmZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listXpZm[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listXpZm[i][j].components.push_back(temp);
        }
        
        //Edge_listXpZp
        Edge_listXpZp[i].resize(cm_perio.Edge_listXpZp->size());
        for (unsigned int j=0; j<cm_perio.Edge_listXmZm->size(); j++) {
            p1 = 0;
            p2 = 2;
            temp.node = (*cm_perio.Edge_listXpZp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listXpZp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listXmZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listXpZp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listXpZp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p2,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p2,i);
            Edge_listXpZp[i][j].components.push_back(temp);
        }
        
        //Edge_listXmZp
        Edge_listXmZp[i].resize(cm_perio.Edge_listXmZp->size());
        for (unsigned int j=0; j<cm_perio.Edge_listXmZm->size(); j++) {
            p1 = 2;
            temp.node = (*cm_perio.Edge_listXmZp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listXmZp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listXmZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listXmZp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listXmZp[i][j].components.push_back(temp);
        }

        //Edge_listYpZm
        Edge_listYpZm[i].resize(cm_perio.Edge_listYpZm->size());
        for (unsigned int j=0; j<cm_perio.Edge_listYmZm->size(); j++) {
            p1 = 1;
            temp.node = (*cm_perio.Edge_listYpZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listYpZm[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listYmZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listYpZm[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listYpZm[i][j].components.push_back(temp);
        }
        
        //Edge_listYpZp
        Edge_listYpZp[i].resize(cm_perio.Edge_listYpZp->size());
        for (unsigned int j=0; j<cm_perio.Edge_listYmZm->size(); j++) {
            p1 = 1;
            p2 = 2;
            temp.node = (*cm_perio.Edge_listYpZp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listYpZp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listYmZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listYpZp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listYpZp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p2,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p2,i);
            Edge_listYpZp[i][j].components.push_back(temp);
        }
        
        //Edge_listYmZp
        Edge_listYmZp[i].resize(cm_perio.Edge_listYmZp->size());
        for (unsigned int j=0; j<cm_perio.Edge_listYmZm->size(); j++) {
            p1 = 2;
            temp.node = (*cm_perio.Edge_listYmZp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Edge_listYmZp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Edge_listYmZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Edge_listYmZp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Edge_listYmZp[i][j].components.push_back(temp);
        }
        
        //Face_listXp
        Face_listXp[i].resize(cm_perio.Face_listXp->size());
        for (unsigned int j=0; j<cm_perio.Face_listXm->size(); j++) {
            p1 = 0;
            temp.node = (*cm_perio.Face_listXp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Face_listXp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Face_listXm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Face_listXp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Face_listXp[i][j].components.push_back(temp);
        }
        
        //Face_listYp
        Face_listYp[i].resize(cm_perio.Face_listYp->size());
        for (unsigned int j=0; j<cm_perio.Face_listYm->size(); j++) {
            p1 = 1;
            temp.node = (*cm_perio.Face_listYp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Face_listYp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Face_listYm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Face_listYp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Face_listYp[i][j].components.push_back(temp);
        }
        
        //Face_listZp
        Face_listZp[i].resize(cm_perio.Face_listZp->size());
        for (unsigned int j=0; j<cm_perio.Face_listZm->size(); j++) {
            p1 = 2;
            temp.node = (*cm_perio.Face_listZp)[j];
            temp.dof = list_dofs[i];
            temp.coef = 1.0;
            Face_listZp[i][j].components.push_back(temp);
            
            temp.node = (*cm_perio.Face_listZm)[j];
            temp.dof = list_dofs[i];
            temp.coef = -1.0;
            Face_listZp[i][j].components.push_back(temp);
            
            temp.node = Constrain_Driver_123[CD_num(p1,i)];
            temp.dof = 1;
            temp.coef = Dxyz(p1,i);
            Face_listZp[i][j].components.push_back(temp);
        }
    }
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const cubic_equation& c_eq)
//--------------------------------------------------------------------------
{
    s << "Display info on the cubic equations:\n";
    
    for (int i=0; i<3; i++) {
        s << "DOF : " << i+1 << "\n";

        for (auto n:c_eq.Face_listXp[i]) {
            s << "Equations between periodic face Xp and Xm :";
            s << n;
        }
        for (auto n:c_eq.Face_listYp[i]) {
            s << "Equations between periodic face Yp and Ym :";
            s << n;
        }
        for (auto n:c_eq.Face_listZp[i]) {
            s << "Equations between periodic face Zp and Zm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listXpYm[i]) {
            s << "Equations between periodic edge XpYm and XmYm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listXpYp[i]) {
            s << "Equations between periodic edge XpYp and XmYm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listXmYp[i]) {
            s << "Equations between periodic edge XmYp and XmYm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listXpZm[i]) {
            s << "Equations between periodic edge XpZm and XmZm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listXpZp[i]) {
            s << "Equations between periodic edge XpZp and XmZm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listXmZp[i]) {
            s << "Equations between periodic edge XmZp and XmZm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listYpZm[i]) {
            s << "Equations between periodic edge YpZm and YpZm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listYpZp[i]) {
            s << "Equations between periodic edge YpZp and YmZm :";
            s << n;
        }
        for (auto n:c_eq.Edge_listYmZp[i]) {
            s << "Equations between periodic edge YmZp and YmZm :";
            s << n;
        }
        
        s << "Equations between periodic corner XpYmZm and XmYmZm :";
        s << c_eq.Corner_listXpYmZm[i];
        s << "Equations between periodic corner XpYpZm and XmYmZm :";
        s << c_eq.Corner_listXpYpZm[i];
        s << "Equations between periodic corner XmYpZm and XmYmZm :";
        s << c_eq.Corner_listXmYpZm[i];
        s << "Equations between periodic corner XmYmZp and XmYmZm :";
        s << c_eq.Corner_listXmYmZp[i];
        s << "Equations between periodic corner XpYmZp and XmYmZm :";
        s << c_eq.Corner_listXpYmZp[i];
        s << "Equations between periodic corner XpYpZp and XmYmZm :";
        s << c_eq.Corner_listXpYpZp[i];
        s << "Equations between periodic corner XmYpZp and XmYmZm :";
        s << c_eq.Corner_listXmYpZp[i];
    }

    return s;
}
    
} //namespace simcoon
