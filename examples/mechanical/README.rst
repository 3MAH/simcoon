Mechanical Constitutive Laws
-----------------------------------------------

Below are examples illustrating Simcoon's mechanical constitutive laws library.

**Elastic Models:**

- **ELISO** - Isotropic elasticity
- **ELIST** - Transversely isotropic elasticity
- **ELORT** - Orthotropic elasticity

**Plasticity Models:**

- **EPICP** - Plasticity with isotropic hardening (power-law)
- **EPKCP** - Plasticity with combined isotropic and kinematic hardening
- **EPCHA** - Plasticity with Chaboche hardening (cyclic plasticity); also runs the dedicated kernel kept under the name ``EPCHA_closed`` next to the default modular engine
- **EPJCK** - Johnson-Cook plasticity (strain-rate sensitivity and thermal softening)

**Viscoelastic Models:**

- **ZENER** - Poynting-Thomson (Zener) model
- **ZENER_N** - Generalized Zener model (N Kelvin branches)
- **PRONY_N** - Prony series (Generalized Maxwell)

**Modular Models (composable ``MODUL`` UMAT):**

- **MODUL** - Elasticity block + von Mises plasticity with Voce hardening
- **MODUL_finite** - The same composition under finite strain (Hencky hyperelasto-plasticity)
- **MODUL_hyper_visco** - Yeoh hyperelastic block + Prony branches (finite-strain viscoelasticity)
- **MODUL_closest_point** - Closest-point projection (``integrator="closest_point"``) vs the cutting-plane loop under a non-proportional figure-eight: accuracy, Newton convergence, cost

**Shape Memory Alloys:**

- **SMA_T** - Superelastic model (transformation only)
- **SMA_TR** - Transformation + martensite reorientation on a non-proportional tension-shear path