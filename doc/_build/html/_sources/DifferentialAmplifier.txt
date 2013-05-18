Differential Amplifier
=======================

The Differential Amplifier with common mode feedback is a fairly complicated 13 dimensional system. 

The state space equations for this system are in the config file `Resources`_

For implementing current equations, we used **level1 MOS** equations. The net list of the circuit is provide in `Resources`_ section

TPWL Integration
----------------

Piece-wise linearization and integration gives the following results. 

.. figure:: ../reports/images/Diffamp_PWL/n001.png
.. figure:: ../reports/images/Diffamp_PWL/n002.png
.. figure:: ../reports/images/Diffamp_PWL/n005.png
.. figure:: ../reports/images/Diffamp_PWL/cfmb2.png
.. figure:: ../reports/images/Diffamp_PWL/vom.png
.. figure:: ../reports/images/Diffamp_PWL/vop.png
.. figure:: ../reports/images/Diffamp_PWL/n009.png
.. figure:: ../reports/images/Diffamp_PWL/n007.png
.. figure:: ../reports/images/Diffamp_PWL/cmfb1.png
.. figure:: ../reports/images/Diffamp_PWL/out1neg.png
.. figure:: ../reports/images/Diffamp_PWL/out1pos.png
.. figure:: ../reports/images/Diffamp_PWL/ig9.png
.. figure:: ../reports/images/Diffamp_PWL/ig7.png


Resources
-----------

   :download:`Differential Amplifier config file <../notebook/lib/config.py>`

   :download:`Net list  <../notebook/Data/diffamp.net>`

