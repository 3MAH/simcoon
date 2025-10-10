The Recovery Props Library
========================

The recovery props library provides a set of function to check and evaluate the properties of stiffness and compliance tensors.

.. default-domain:: cpp

.. function:: void check_symetries(mat L, string umat_type, int axis, vec props, int maj_sym)

    Check the symmetries of a stiffness matrix L, and fill the vector of material properties. 
    Depending on the symmetry found, the string umat_type, the axis of symmetry (if applicable) the vector of material properties, and the  major symmetry maj_sym (L_ij = L_ji ?).
    If the major symmetry condition is not fulfilled, the check of symmetries if performed on the symmetric part of L
    
	.. csv-table:: Material Symmetries considered
   		:header: "Symmetry", "umat_type", "axis", "size of props"
   		:widths: 60, 20, 30,20

		"Fully anisotropic", "ELANI", "0", "0"
   		"Monoclinic", "ELMON", "1,2 or 3", "0"   
   		"Orthotropic", "ELORT", "0", "9"      
   		"Cubic", "ELCUB", "0", "3"  
   		"Transversely isotropic", "ELITR", "1,2 or 3", "5"  
   		"Isotropic", "ELISO", "0", "2"        
    
    .. code-block:: cpp

        check_symetries(L, umat_name, axis, props, maj_sym);