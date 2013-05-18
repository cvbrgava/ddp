4 input NAND Gate cascaded with Invertor
==========================================

This is the first step towards working on cascaded systems. There are two ways of working on cascaded systems.
	
	- Build the complete state space and reduce
	- Break the cascaded system into component systems and reduce each system. 

We have now followed the first procedure, but will change the code to include the second functionality.

The state space equations for this system are in the config file `Resources`_

For implementing current equations, we used **level1 MOS** equations. The net list of the circuit is provide in `Resources`_ section


TPWL Integration
----------------

Piece-wise linearization and integration gives the following results. 

.. figure:: ../reports/images/4inpNAND_inv1.png


TPWL-MOR Integration 
----------------------

Once the model TPWL model is verified, we go ahead with the MOR using the Moment matching method with orthogonal projection. The following is the response obtained with MOR of the system. 

The error obtained from the MOR system, is heavily dependant on the order of the reduced system. The following plots depict the difference observed.

.. figure:: ../notebook/Data/4inpNAND/4inpNAND_inv_voltageinvout_2_3.png

.. figure:: ../notebook/Data/4inpNAND/4inpNAND_inv_voltageinvout_3_2.png

Resources
-----------

   :download:`4 input NAND Gate cascaded with Invertor config file <../notebook/config/config_4inpNAND_inv.py>`

   :download:`Net list  <../notebook/Data/4inpNAND_inv.net>`

