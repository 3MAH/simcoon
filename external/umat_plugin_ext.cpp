#include <iostream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>

using namespace std;
using namespace arma;

#if defined(_WIN32) || defined(_WIN64)
    #define LIB_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) || defined(__clang__)
    #if __GNUC__ >= 4
        #define LIB_EXPORT __attribute__((visibility("default")))
    #else
        #define LIB_EXPORT
    #endif
#else
    #define LIB_EXPORT
#endif

class LIB_EXPORT umat_plugin_ext : public umat_plugin_ext_api {
public:

    std::string name() const {
        return "umext";
    }

    void umat_external_M(const std::string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt) {

        UNUSED(umat_name);
        UNUSED(Etot);
        UNUSED(DR);
        UNUSED(nprops);
        UNUSED(nstatev);
        UNUSED(statev);
        UNUSED(Time);
        UNUSED(DTime);
        UNUSED(nshr);
        UNUSED(tnew_dt);

        double T_init = statev(0);

        //From the props to the material properties
        double E = props(0);
        double nu = props(1);
        double alpha = props(2);

        //Elastic stiffness tensor
        L = simcoon::L_iso(E, nu, "Enu");

        ///@brief Initialization
        if(start)
        {
            T_init = T;
            sigma = zeros(6);

            Wm = 0.;
            Wm_r = 0.;
            Wm_ir = 0.;
            Wm_d = 0.;
        }

        vec sigma_start = sigma;

        //Compute the elastic strain and the related stress
        vec Eel = Etot + DEtot - alpha*(T+DT-T_init);
        sigma = simcoon::el_pred(L, Eel, ndi);

        Lt = L;

        //Computation of the mechanical and thermal work quantities
        Wm += 0.5*sum((sigma_start+sigma)%DEtot);
        Wm_r += 0.5*sum((sigma_start+sigma)%DEtot);
        Wm_ir += 0.;
        Wm_d += 0.;

        statev(0) = T_init;
    }
    
};

extern "C" LIB_EXPORT umat_plugin_ext_api* create_api()
{
    return new umat_plugin_ext();
}

extern "C" LIB_EXPORT void destroy_api(umat_plugin_ext_api* p)
{
    delete p;
}


