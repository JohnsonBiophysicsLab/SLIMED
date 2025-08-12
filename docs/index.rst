.. Continuum membrane model documentation master file, created by 
   M. Ying on Aug. 4, 2022.

Documentation for NERDSS - Continuum Membrane & Dynamics
=============================================

Continuum membrane & dynamics is tool in NERDSS to extend the scope of the model to non-linear membranes established with triangular mesh and optimized using and energy function. 

1. Installation
------------

Continuum membrane & dynamics is currently an independent model but intended to be integrated in NERDSS in the future. Please refer to `User guide`_ for NERDSS. Additional dependencies for continuum membrane model are listed below. Note that continuum membrane & dynamics model no longer needs `Armadillo`_ linear algebra package to run.

#. C++ 11
#. `GSL`_ (GNU scientific library)
#. OpenMP (only required in need of parallelization)

2. Compile Continuum Membrane
--------------------------
To compile and run continuum membrane lowest energy conformation search with OpenMP parallelization. See section 5 for model details.

.. code-block:: console

  make omp
  ./bin/continuum_membrane

* if parallelization is not needed:

   .. code-block:: console
   
      make serial
      ./bin/continuum_membrane
      
To compile and run membrane Brownian dynamics simulation. See section 6 for model details.

   .. code-block:: console
   
      make dyna
      ./bin/membrane_dynamics

3. Input Parameters
----------------------

The input file is ``continuum_membrane/input.params``. Parameters are broken down to geometric parameters, physical properties, insertion mode, and advanced parameters.

+--------------------------------------------------+----------------------------------------------------------------------------------+
| Parameter                                        | Description                                                                      |
+------------+-------------------------------------+----------------------------------------------------------------------------------+
| Geometric  | ``lMeshSide``                       | Target side length of the triangular mesh (nm).                                  |
| Parameters |                                     | This only servers as a reference scale.                                          |
|            |                                     | The mesh side length set up by the algorithm may vary.                           |
|            +--------------+----------------------+----------------------------------------------------------------------------------+
|            | Sphere model | ``isSphere``         | Set ``true`` to enable sphere mode.                                              |
|            |              +----------------------+----------------------------------------------------------------------------------+
|            |              | ``rSphere``          | Target radius of sphere (nm).                                                    |
|            |              |                      | This is the radius of spherical frame to set up the triangular mesh.             |
|            |              |                      | The radius of the resulting membrane represented by the triangular mesh may vary |
+------------+--------------+----------------------+----------------------------------------------------------------------------------+
| Physical   | ``c0Insertion``                     | Curvature of the membrane at the insertion area.                                 |
| Properties +-------------------------------------+----------------------------------------------------------------------------------+
|            | ``c0Membrane``                      | Spontaneous curvature of the membrane.                                           |
|            +-------------------------------------+----------------------------------------------------------------------------------+
|            | ``kcMembraneBending``               | Membrane bending constant in the energy function (pN*nm).                        |
|            +-------------------------------------+----------------------------------------------------------------------------------+
|            | ``usMembraneStretching``            | Membrane streching modulus in the energy function (pN/nm).                       |
|            +-------------------------------------+----------------------------------------------------------------------------------+
|            | ``uvVolumeConstraint``              | Volume constraint coefficient in the energy function (pN/nm^2).                  |
+------------+-------------------------------------+----------------------------------------------------------------------------------+
| Insertion  | ``isInsertionIncluded``             | Set ``true`` to include insertion.                                               |
| Mode       +-------------------------------------+----------------------------------------------------------------------------------+
|            | ``sigma``                           | 2*sigma (nm) is the length scale of decaying insertion curvature,                |
|            |                                     | or in other words expansion of non-spontaneous curvature due to insertion.       |
+------------+--------------+----------------------+----------------------------------------------------------------------------------+
| Advacned   | Optimization | ``numMaxIterations`` | Number of maximum iterations allowed.                                            |
| Parameters |              +----------------------+----------------------------------------------------------------------------------+
|            |              | ``criterionForce``   | Force criteria to determine if adequate optimization is accomplished (pN).       |
|            +--------------+----------------------+----------------------------------------------------------------------------------+
|            | Algorithm    | ``gaussQuadratureN`` | Default Gauss Quadrature used in integral approximation.                         |
+------------+--------------+----------------------+----------------------------------------------------------------------------------+

4. Triangular Mesh Setup
-----------------------
The first step for continuum membrane is to set up the triangular mesh model to approximate the geometry of the membrane. A brief framework is generated by dividing the geometric framework given by the geometric parameters (such as ``rSphere`` in sphere mode) into large triangular cells. Next, Loop's  subdivision method (`F. Cirak et al., 2000`_) is applied to further divide the brief framework into smaller cells to better approximate the given geometry.

5. Energy Function and Lowest Energy Search
-----------------------

The goal for the lowest energy search model is to minimize the membrane energy evaluated by the energy function, which is the sum of membrane bending energy, area constraint energy (or elastic area change energy), and volume constraint energy:

.. math::

   E = E_B + E_S + E_V = \int_S \frac{1}{2}\kappa (2H-C_0)^2 dS + \frac{1}{2} \mu_S \frac{(S-S_0)^2}{S_0} + \frac{1}{2} \mu_V \frac{(V-V_0)^2}{V_0}

where:

- :math:`\kappa` : membrane bending constant ``kcMembraneBending``
- :math:`H` : mean membrane culvature
- :math:`C_0` : spontaneous curvature of the membrane ``c0Membrane``
- :math:`\mu_S` : membrane streching modulus ``usMembraneStretching``
- :math:`S` : global membrane area
- :math:`S_0` : target membrane area
- :math:`\mu_V` : volume constraint coefficient ``uvVolumeConstraint``
- :math:`V` : global volume
- :math:`V_0` : target volume
 
6. Membrane Brownian Dynamics
-----------------------

Membrane Brownian Dynamics model runs a step-wise simulation of the moving membrane surface with the following equation:

.. math::

   \Delta X = -\frac{D\Delta t}{k_b T} \nabla E + \sqrt{2D\Delta t} (N(0,1))

where:

- :math:`\Delta X`: displacement of point on limit surface
- :math:`D`: diffusion constant of the membrane
- :math:`\Delta t`: time step
- :math:`N(0,1)`: standard normal distribution

Note that the displacement of membrane according to the equation above is performed on the limit surface, not the control mesh.
In this case, a conversion matrix helps to convert between triangular mesh and limit surface, as currently the points on the limit surface
represented by the mesh point are chosen to represent the surface.

.. math::

   M_{s} = C M_{m}

.. math::

   M_{m} = C^{-1} M_{s}

7. Boundary Conditions
-----------------------

Three types of boundary conditions are provided currently in both models. Note that "ghost vertices" are defined as points on the boundary of the triangular mesh that only serve to provide reference when calculating limit surface on the boundary, as calculating position of a point on the limit surface require the coordinates of 12 neighboring vertices (if regular). However, the "ghost vertices" themselves do not correspond to real points on the surface.

- Fixed: 2 rings of ghost vertices are fixed in space
- Periodic: 3 rings of ghost vertices that mimics the movement of the vertices on the opposite side of the membrane.
- Free: 2 rings of ghost vertices are generated after movement by forming parallelogram extend from the real points on the control mesh

8. Cite Continuum Membrane
-----------------------

If you use or modify continuum membrane model, in addition to citing NERDSS, please be kind and cite us:

1. Continuum Membrane Implementation
Fu, Y., Yogurtcu, O.N., Kothari, R., Thorkelsdottir, G., Sodt, A.J. & Johnson, M.E. (2019) An implicit lipid model for efficient reaction-diffusion simulations of protein binding to surfaces of arbitrary topology. *J Chem Phys.* 151 (12), 124115. doi:`10.1063/1.5120516`_

2. Membrane energies and insertion
Fu, Y., Zeno, W., Stachowiak, J. & Johnson, M.E. A continuum membrane model predicts curvature sensing by helix insertion. Submitted (2021) Available on `bioRxiv`_

.. _`User guide`: https://github.com/mjohn218/NERDSS/blob/master/NERDSS_USER_GUIDE.pdf
.. _`Armadillo`: http://arma.sourceforge.net/
.. _`GSL`: https://www.gnu.org/software/gsl/
.. _`10.1063/1.5120516`: https://pubmed.ncbi.nlm.nih.gov/31575182/
.. _`bioRxiv`: https://www.biorxiv.org/content/10.1101/2021.04.22.440963v1.full
.. _`F. Cirak et al., 2000`: http://multires.caltech.edu/pubs/thinshell.pdf
