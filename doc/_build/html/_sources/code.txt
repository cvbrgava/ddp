Auto generated documentation
=============================
numericals module
-----------------
This module contains all the functions dealing with integration generation of the state space matrices

.. automodule:: numericals
   :members:

currenteqs module
------------------
This module contains current equations for both PMOS and NMOS transistors. They are based on LEVEL1 MOS current equation models. These have been verified before adding to this file.

.. automodule:: currenteqs
   :members:

config module
-------------------
This is the customizable module. The file contains state space equations and details regarding the SPICE files used in simulation.

.. automodule:: config
   :members:

getparameters module
---------------------
This is the interface between project and the SPICE files. From the spice waveform files and netlist this module extracts all the parameters needed for the project to function.

.. automodule:: getparameters
   :members:

Parsers module 
----------------
This module can be used to parse netlist and equations.

.. automodule:: parsers
   :members:

MOR by Kalman Decomposition and TBR module 
--------------------------------------------

This module was developed mainly to provide a platform for Kalman Decomposition and Truncated Balanced Reduction. In the process of which the functions developed here span a wide range of control system analysis procedures.

.. automodule:: kalman_tbr
   :members:


MOR by Moment Matching module
-----------------------------

This module contains all the functions to be used for MOR by Moment matching

.. automodule:: arnoldi
   :members:

