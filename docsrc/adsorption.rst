:orphan:

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <div class="col-md-10">
        <div style="text-align: justify; text-justify: inter-word;">

Adsorption Isotherms
====================

.. code-block:: python

    import hgana as hga

Initialize System
-----------------

.. code-block:: python

    ads = hga.Adsorption([10, 10, 10])

Add molecules

.. code-block:: python

    ads.add_mol([1, 10, 20], is_move=False)
    ads.add_mol([2, 3, 4, 5]+[x for x in range(1, 100+1, 20)])

Set interactions

.. code-block:: python

    ads.set_interaction(0, 1, -15)

Run Simulation
--------------

.. code-block:: python

    results = ads.run(298, 100000, 10000, "output/ads.obj", binding=[{"host": 0, "guest": 1}], pb_f=[1000, 50], n_print=1000, is_parallel=True)

``Running equilibration phase ...``

``Running production phase ...``

``10000/10000  - acc/rej=2.321, p_b(0,1)=0.9652+-0.0397``


Plot Results
------------

.. code-block:: python

    ads.plot("output/ads.obj", 1, 0, (0, 1))



.. figure::  /pics/adsorption.svg
  :align: center
  :width: 70%
  :name: fig1


.. raw:: html

        </div>
      </div>
    </div>
  </div>
