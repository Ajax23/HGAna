:orphan:

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <div class="col-md-10">
        <div style="text-align: justify; text-justify: inter-word;">

Binding Affinity
================

.. code-block:: python

  import pandas as pd
  import moldyn as md

.. code-block:: python

  md.cd.affinity.hist("COLVAR", [1, 2], ["Centers of Mass", "Oxygenes"])

.. figure::  /pics/affinity.svg
  :align: center
  :width: 50%
  :name: fig1


.. code-block:: python

  md.cd.affinity.sample("COLVAR", "count.obj", [1, 0.5], [2, 0.46], is_force=True)


.. code-block:: python

  md.cd.affinity.number("count.obj", 298.15, 31.3707e-27)

.. raw:: html

  <div class="nboutput nblast">
    <div class="output_area rendered_html">
      <table border="1" class="dataframe">
        <thead>
          <tr style="text-align: right;">
            <th></th>
            <th>kJ/mol</th>
            <th>kcal/mol</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th>dG</th>
            <td>-14.028719</td>
            <td>-3.351099</td>
          </tr>
          <tr>
            <th>dG_O1</th>
            <td>-14.022662</td>
            <td>-3.349652</td>
          </tr>
          <tr>
            <th>dG_O2</th>
            <td>0.882953</td>
            <td>0.210915</td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>


.. code-block:: python

  tables = [md.cd.affinity.time("data/cyclodextrin/count.obj", cutoff, 298.15, 31.3707e-27) for cutoff in [100*x for x in range(11)]]

  table = pd.concat(tables)

  display(table)

.. raw:: html

  <div class="nboutput nblast">
    <div class="output_area rendered_html">
      <table border="1" class="dataframe">
        <thead>
          <tr style="text-align: right;">
            <th></th>
            <th>Cutoff (ps)</th>
            <th>dG (kJ/mol)</th>
            <th>dG (kcal/mol)</th>
            <th>k_on (dm^3/mol/s)</th>
            <th>k_off (1/s)</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th>0</th>
            <td>0</td>
            <td>-14.035023</td>
            <td>-3.352605</td>
            <td>1.073377e+10</td>
            <td>3.729713e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>100</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>200</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>300</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>400</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>500</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>600</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>700</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>800</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>900</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
          <tr>
            <th>0</th>
            <td>1000</td>
            <td>-8.353826</td>
            <td>-1.995514</td>
            <td>3.106633e+08</td>
            <td>1.068000e+07</td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>


.. raw:: html

        </div>
      </div>
    </div>
  </div>
