========================================================================
    Fortran Console Application : "Adams_Fixed_Step" Project Overview
========================================================================

The Intel Fortran Console Application Wizard has created this 
"Adams_Fixed_Step" project for you as a starting point.

This file contains a summary of what you will find in each of the files 
that make up your project.

Adams_Fixed_Step.vfproj
    This is the main project file for Fortran projects generated using an 
    Application Wizard.  It contains information about the version of 
    Intel Fortran that generated the file, and information about the 
    platforms, configurations, and project features selected with the 
    Application Wizard.

Adams_Fixed_Step.f90
    This is the main source file for the Fortran Console application. 
    It contains the program entry point.

accl_potl_module.f90
    This is a module to compute the accelaration and gravity potential of
    one satellite for one given epoch. The accelaration is caused by gravity,
    relativistic effect and n-body perturbation.

    Identifier:
        type(sate_char):
            This is a derived structure illustrating the dynamic characteristics
            including position, velocity, accelaration and potential.
        gl:
            This is a derived structure too for assigning the structure of the 
            leading satellite.
        gt:
            This is a derived structure too for assigning the structure of the 
            tracking satellite.
coor_trans(directory)
    This is a directory which contains subroutines to transform coordinates
    to coordinates.

ddeabm_module.f90:
    This is a module which contains the adams-bashford-moulton integrator 13rd.

/////////////////////////////////////////////////////////////////////////////
Other notes:
    The ddeabm module comes from https://github.com/sausage02/Adams-Bashford.
    The original ddeabm FORTRAN 77 code comes from http://www.netlib.org/slatec/src/.

    If you have any questions, please contact chelsea_blue@foxmail.com.
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
Reference:
    1. L. F. Shampine, M. K. Gordon, "Solving ordinary differential equations
        with ODE, STEP, and INTRP", Report SLA-73-1060, Sandia Laboratories, 1973.
    2. L. F. Shampine, M. K. Gordon, "Computer solution of ordinary differential
        equations, the initial value problem", W. H. Freeman and Company, 1975.
    3. L. F. Shampine, H. A. Watts, "DEPAC - Design of a user oriented package of 
        ode solvers", Report SAND79-2374, Sandia Laboratories, 1979.
    4. H. A. Watts, "A smoother interpolant for DE/STEP, INTRP and DEABM: II", 
        Report SAND84-0293, Sandia Laboratories, 1984.
    5. R. P. Brent, "An algorithm with guaranteed convergence for finding a zero 
        of a function", The Computer Journal, Vol 14, No. 4., 1971.
    6. R. P. Brent, "Algorithms for minimization without derivatives", 
        Prentice-Hall, Inc., 1973.
/////////////////////////////////////////////////////////////////////////////

