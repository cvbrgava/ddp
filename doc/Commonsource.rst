Common Source Amplifier
=======================

A common source is a very simple system with 2 state variables and hence, a 2D phase space. The results obtained from this system will be easy to visualize in such a state space.

.. figure:: ../reports/images/commonsource.png
   :align: center
   :width: 300px
   :height: 300px

   **Common Source Amplifier**

The state space equations for this system are as follows.

   :math:`\dot{v} & =\left(\left(\frac{V_{dd}-v}{R_{3}}\right)-I_{d}\left(M1\right)-I_{C}\right)*\frac{1}{C_{para}}`

   :math:`\dot{I_{C}} & =\left(\frac{\left(\left(\frac{V_{dd}-v}{R_{3}}\right)-I_{d}\left(M1\right)-I_{C}\right)*\frac{1}{C_{para}}-I_{c}}{R_{4}}\right)`

For implementing current equations, we used **level1 MOS** equations which can be found in `Resources`_

Non-Linear Integration with state-based boundaries
---------------------------------------------------

Non-linear models with transistions made using proximity to points where the models are built. The following are the results obtained. 

.. note:: This is the result from an older version of the program.

.. figure:: ../reports/images/cmnsrc_nonl_intg_state_bound.png
 


TPWL Integration
----------------

Piece-wise linearization and integration gives the following results. 

.. figure:: ../reports/images/cmnsrcModel3.png

The validity of the TPWL models is verified by the response to inputs other than the test input. For all the figures here, red dots represent the linearization points and the solid green line is the SPICE response, while the one in broken lines is the response from TPWL model. The linearization points are for the test input with 1 us rise time.

The following is the response to a 2 us rise time input.

.. figure:: ../notebook/Data/cmnsrc/2urisedrain.png

The response to a 3 us rise time input.

.. figure:: ../notebook/Data/cmnsrc/3urisedrain.png

The response to a 0.25 us rise time input.

.. figure:: ../notebook/Data/cmnsrc/0_25urisedrain.png

The response to a sinusoidal input.

.. figure:: ../notebook/Data/cmnsrc/sin1000kdrain.png

Resources
-----------

   :download:`common source config file <../notebook/config/config_cmnsrc.py>`

   :download:`level1 MOS library <../lib/level1MOS.lib>`


