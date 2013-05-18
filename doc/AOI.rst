And-Or-Invert Gate
==================

This system has a 3 dimensional state space.

The state space equations for this system are in the config file `Resources`_

For implementing current equations, we used **level1 MOS** equations. The net list of the circuit is provide in `Resources`_ section


TPWL Integration
----------------

Piece-wise linearization and integration gives the following results. 

.. figure:: ../reports/images/AOIgaten002.png


TPWL-MOR Integration 
----------------------

Once the model TPWL model is verified, we go ahead with the MOR using the Moment matching method with orthogonal projection. The following is the response obtained with MOR of the system. 

The error obtained from the MOR system, is heavily dependant on the order of the reduced system. The following plots depict the difference observed.

.. figure:: ../reports/images/AOI_voltagen002_1_2.png

.. figure:: ../reports/images/AOI_voltagen002_2_2.png

.. figure:: ../reports/images/AOI_voltagen002_3_4.png



Resources
-----------

   :download:`Invertor cascaded with 4 input NAND Gate config file <../notebook/config/config_AOI.py>`

   :download:`Net list  <../notebook/Data/AOI.net>`
