Motivation for Order Reduction
===============================

This is a project which merges Dynamical systems with Electrical Engineering. We here realize that naive choice of state variables is not computationally efficient. We can come up with a linear combination of state variables which make describing the system much more efficient than it is right now.

.. figure:: ../reports/images/osciltr.png
   :align: left

   Voltage Controlled Differential Oscillator

.. figure:: ../reports/images/2Dphasespace.png
   :align: center
   :width: 350px
   :height: 350px

   State space of differential oscillator

This system has a 3D state space but is actually a **2D annular ring in 3D space**. So we would not need all the three dimensions to describe this system. Rather a plane roughly passing through the annular ring will capture the details of the circuit.

There are some circuits like this one, which inherently have a solution spanning a lower dimensional space in the actual state space. It is evident that dimensionality of such circuits can be reduced. 

.. figure:: ../reports/images/commonsource.png
   :align: left
   :width: 300px
   :height: 300px

   Commonsource amplifier

.. figure:: ../reports/images/1Dphasespace.png
   :align: center
   :width: 300px
   :height: 300px

   State space of common source amplifier

.. :math:`\alpha > \beta`

There are few others which don;t inherently behave in the similar fashion, instead they show this feature under specific circumstances. The *common source amplifier*, which is a **2-Dimensional** circuit, produces an almost **1-Dimensional** response to a step input. 
This implies that we can again exploit this functionality and reduce the system for this particular input. Hence the name **Trajectory based Order reduction.** 

.. seealso:: `Non-Linear Dynamics and Chaos <http://libgen.org/book/index.php?md5=6931A869174964929BA8726903086767>`_ by Steven Strogatz on understanding systems from state space perspective.
   
