=begin
= ALF: a collection of routines for computing associated Legendre polynomials (ALFs)
ALF is an extension library for GSL to compute associated Legendre polynomialsdeveloped by Patrick Alken.
Ruby/GSL includes interfaces to it if ALF is installed found by extconf.

The class and method descriptions below are based on references from the document of ALF (alf-1.0/doc/alf.texi) by P.Alken.

== Module structure
  GSL::ALF (module)
       GSL::ALF::Workspace (Class)

== Creating ALF workspace
--- GSL::ALF::Workspace.alloc(lmax)
--- GSL::ALF.alloc(lmax)
Creates a workspace for computing associated Legendre polynomials (ALFs). The maximum ALF degree is specified by lmax. The size of this workspace is O(lmax).

== Methods
--- GSL::ALF::Workspace#params(csphase, cnorm, norm)
    Sets various parameters for the subsequent computation of ALFs. If ((|csphase|)) is set to a non-zero value, the Condon-Shortley phase of (-1)^m will be applied to the associated Legendre functions. The Condon-Shortley phase is included by default. If ((|cnorm|)) is set to zero, the real normalization of the associated Legendre functions will be used. The default is to use complex normalization. The norm parameter defines the type of normalization which will be used. The possible values are as follows.
    * ALF::NORM_SCHMIDT: Schmidt semi-normalized associated Legendre polynomials S_l^m(x). (default)
    * ALF::NORM_SPHARM: Associated Legendre polynomials Y_l^m(x) suitable for the calculation of spherical harmonics.
    * ALF::NORM_ORTHO: Fully orthonormalized associated Legendre polynomials N_l^m(x).
    * ALF::NORM_NONE:: Unnormalized associated Legendre polynomials P_l^m(x).
--- GSL::ALF::Workspace#Plm_array(x)
--- GSL::ALF::Workspace#Plm_array(lmax, x)
--- GSL::ALF::Workspace#Plm_array(x, result)
--- GSL::ALF::Workspace#Plm_array(lmax, x, result)
--- GSL::ALF::Workspace#Plm_array(x, result, deriv)
--- GSL::ALF::Workspace#Plm_array(lmax, x, result, deriv)
    Compute all associated Legendre polynomials P_l^m(x) and optionally their first derivatives dP_l^m(x)/dx for 0 <= l <= lmax, 0 <= m <= l. The value of ((|lmax|)) cannot exceed the previously specified lmax parameter to (({ALF.alloc})), but may be less. If ((|lmax|)) is not given, the parameter to (({ALF.alloc()})) is used. The results are stored in ((|result|)), an instance of (({GSL::Vector})). Note that this vector must have enough length to store all the values for the polynomial P_l^m(x), and the length required can be known using (({ALF::array_size(lmax)})). If a vector is not given, a new vector is created and returned. 
     The indices of ((|result|)) (and ((|deriv|)) corresponding to the associated Legendre function of degree ((|l|)) and order ((|m|)) can be obtained by calling (({ALF::array_index(l, m)})).
--- GSL::ALF::Workspace#Plm_deriv_array(x)
--- GSL::ALF::Workspace#Plm_deriv_array(lmax, x)
--- GSL::ALF::Workspace#Plm_deriv_array(x, result, deriv)
--- GSL::ALF::Workspace#Plm_deriv_array(lmax, x, result, deriv)
    Compute all associated Legendre polynomials P_l^m(x) and their first derivatives dP_l^m(x)/dx for 0 <= l <= lmax, 0 <= m <= l. 

--- GSL::ALF::array_size(lmax)
    Returns the size of arrays needed for the array versions of P_l^m(x).
--- GSL::ALF::array_index(l, m)
    Returns the array index of results of (({Plm_array()})) and (({Plm_deriv_array()})) corresponding to P_l^m(x) and dP_l^m(x)/dx respectively. The index is given by l(l + 1)/2 + m.

=end
