#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <armadillo>

namespace bn = boost::python::numpy;
namespace bp = boost::python;

namespace simpy {

	//======================================
	class Umat_fedoo
		//======================================
	{
	private:
		int id_umat; //id of the umat
		arma::mat list_props; //mat where each col is prop for one pg
		arma::mat list_statev; //mat where each col is prop for one pg		
		arma::mat list_etot, list_Detot, list_kirchoff, list_PKII, list_cauchy;
		arma::cube list_DR;
		arma::mat list_Wm; //Energy [Wm, Wm_r, Wm_ir, Wm_d]
		arma::cube list_Lt, list_L, list_Lt_start; //tangeant matrix and elastic matrix
		arma::mat list_statev_start, list_Wm_start, list_PKII_start; 
		arma::cube listF0, listF1;

	protected:

	public:
		int nprops, nstatev; //number of properties and state variable
		int nb_points; //number of material points (Gauss points in FE analysis) in wich the umat is applied
		int corate; //corate (0 = Jaumann, 1 = Green Naghdi, 2 = Log strain with cumulated strain, 3 = Log strain with recomputed strain)
		int ndi, nshr, ncomp; //number of normal strain component, and shear strain compotenent (default = 3) and total number of strain component
		double Time,DTime;
		bool nlgeom = false;

		Umat_fedoo(const std::string&, const bn::ndarray&, const int&, const int&, const int&); 	//default constructor

		void Initialize(const double&, bn::ndarray&, const bool&);
		void set_start();
		void to_start();
		void compute_Detot(const double&, const bn::ndarray&);
		bp::tuple Run(const double&);
		bn::ndarray Get_props();
		bn::ndarray Get_PKII();
		bn::ndarray Get_Cauchy();
		bn::ndarray Get_Kirchhoff();
		bn::ndarray Get_etot();
		bn::ndarray Get_Detot();
		bn::ndarray Get_statev();
		bn::ndarray Get_L();
		bn::ndarray Get_Lt();
		bn::ndarray Get_DR();
		bn::ndarray Get_Wm();
		bn::ndarray Get_F0();
		bn::ndarray Get_F1();
																															//phase_multi(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&); //Constructor with parameters
		//phase_multi(const phase_multi&);	//Copy constructor
		//virtual ~phase_multi();

		//void to_start();
		//void set_start();

		//virtual phase_multi& operator = (const phase_multi&);

		//friend std::ostream& operator << (std::ostream&, const phase_multi&);
	};

}