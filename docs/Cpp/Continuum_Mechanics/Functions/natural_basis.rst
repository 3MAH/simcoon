Natural Basis
====================

.. default-domain:: cpp

.. cpp:class:: natural_basis

    A class that provides the basis vectors (covariant and contravariant) for the natural frame of reference, with the associated tools to consider its evolution in a differential variety

    .. code-block:: cpp

        //======================================
        class natural_basis
        //======================================
        {
            private:

            protected:

            public :
            
                std::vector<arma::vec> g_i; // Covariant Vectors
                std::vector<arma::vec> g0i; // Contravariant Vectors
            
                arma::mat g_ij; // Covariant components of the metric tensor
                arma::mat g0ij; // Contravariant components of the metric tensor
                
                natural_basis();     //default constructor
                natural_basis(const std::vector<arma::vec> &); //Constructor with parameters
                natural_basis(const natural_basis &);    //Copy constructor
                virtual ~natural_basis();
            
                virtual void update(const std::vector<arma::vec> &); //update with a new set of covariant vectors
                virtual void from_F(const arma::mat &F); //update using the transformation gradient
            
                virtual natural_basis& operator = (const natural_basis&);
                
                friend std::ostream& operator << (std::ostream&, const natural_basis&);
        };

The natural basis class provides the objects that allows to work within a material system coordinatea, i.e.

.. cpp:member:: std::vector<arma::vec> g_i

The three *covariant* vectors
