4 input NAND Gate
=======================

The 4 input NAND Gate is a 4 dimensional system, with a relatively redundate state space equations.

The state space equations for this system are in the config file `Resources`_

For implementing current equations, we used **level1 MOS** equations. The net list of the circuit is provide in `Resources`_ section

The observed ouput of the sytem is provided here.

.. figure:: ../reports/images/4inpNAND.png

TPWL Integration
----------------

Piece-wise linearization and integration gives the following results. 

.. figure:: ../reports/images/4inpNANDModel3.png

The validity of the model is checked by giving input different from the test input used for generation of the linearization points. In the following plots, the solid line represents the SPICE output and the broken line represents the ouput from TPWL model.

.. figure:: ../notebook/Data/4inpNAND/4inpNAND/0_5urisen001.png

.. figure:: ../notebook/Data/4inpNAND/4inpNAND/2urisen001.png


TPWL-MOR Integration 
----------------------

Once the model TPWL model is verified, we go ahead with the MOR using the Moment matching method with orthogonal projection. The following is the response obtained with MOR of the system. 

The linearization for all the plots has been done using the step function with a 1us rise time. The response for the MOR model has been calculated with inputs of different rise time. The following compares 0.25 us rise time step response and 2 us rise time step response for MOR model.

.. figure:: ../notebook/Data/4inpNAND/4inpNAND/0_25urisen001_2_2.png

The validity of the MOR system is verified using an input considerably different from the reference input used to build the linearization points and hence the complete MOR system. Here the response of the MOR system built using a 1us rise time step input, with a sinusoidal input has shown.

.. figure:: ../notebook/Data/4inpNAND/4inpNAND/sin1000k_2_0n001_2_2.png

Resources
-----------

   :download:`4 input NAND gate config file <../notebook/config/config_4inpNAND.py>`

   :download:`Net list  <../notebook/Data/4inpNAND.net>`

