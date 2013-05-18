NAND implementation of XOR gate
==========================================

This is a 8 dimensional system.

The state space equations for this system are in the config file `Resources`_

For implementing current equations, we used **level1 MOS** equations. The net list of the circuit is provide in `Resources`_ section


TPWL Integration
----------------

Piece-wise linearization and integration gives the following results. 

.. figure:: ../reports/images/XOR_voltage_xorout.png


TPWL-MOR Integration 
----------------------

Once the model TPWL model is verified, we go ahead with the MOR using the Moment matching method with orthogonal projection. The following is the response obtained with MOR of the system. 

The error obtained from the MOR system, is heavily dependant on the order of the reduced system. The following plots depict the difference observed.

.. figure:: ../reports/images/XOR_voltagexorout_4_2.png

.. figure:: ../reports/images/XOR_voltagexorout_7_2.png




Resources
-----------

   :download:`NAND implementation of XOR gate config file <../notebook/config/config_xor.py>`

   :download:`Net list  <../notebook/Data/XOR.net>`
