Invertor cascaded with 4 input NAND Gate
==========================================

This system has a 5 dimensional state space and is similar to the one in earlier section.

The state space equations for this system are in the config file `Resources`_

For implementing current equations, we used **level1 MOS** equations. The net list of the circuit is provide in `Resources`_ section


TPWL Integration
----------------

Piece-wise linearization and integration gives the following results. 

.. figure:: ../reports/images/inv_4inpNAND_TPWL1.png


TPWL-MOR Integration 
----------------------

Once the model TPWL model is verified, we go ahead with the MOR using the Moment matching method with orthogonal projection. The following is the response obtained with MOR of the system. 

The error obtained from the MOR system, is heavily dependant on the order of the reduced system. The following plots depict the difference observed.

.. figure:: ../notebook/Data/4inpNAND/inv_4inpNAND_voltagenandout_2_3.png

.. figure:: ../notebook/Data/4inpNAND/inv_4inpNAND_voltagenandout_3_2.png

.. figure:: ../notebook/Data/4inpNAND/inv_4inpNAND_voltagenandout_4_2.png



Resources
-----------

   :download:`Invertor cascaded with 4 input NAND Gate config file <../notebook/config/config_inv_4inpNAND.py>`

   :download:`Net list  <../notebook/Data/inv_4inpNAND.net>`
