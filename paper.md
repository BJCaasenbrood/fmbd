---
title: '\texttt{fmbd}: a compact formulation of flexible multibody dynamics with rigid and reduced elastic bodies'
author:
  - B. J. Caasenbrood
  - Claude (Opus 5.5, Anthropic)
date: 15 November 2026
abstract: |
  We present the formulation behind fmbd, a compact MATLAB toolbox for the dynamics of multibody systems composed of rigid and linearly elastic bodies. The bodies are defined by surface geometries and meshed with tetrahedral finite elements (\texttt{GMsh}). Elastic bodies are described under a moving body frame and reduced by Craig–Bampton component mode synthesis; all mass-related terms are assembled from inertia invariants computed once, so that the cost per time step is mesh-size independent. Bodies are connected by joints and transmissions, written as holonomic constraints with rotation constraints based on axis dot products. The resulting differential-algebraic equations (DEA) are integrated with a Lie group generalised-α method, solved for static equilibrium with a constraint-split Newton method. Furthermore, models can be linearised to a modal state-space model for frequency analysis and control design. Every equation in this paper is cited by number in the source code.
bibliography: paper.bib
csl: acm.csl
link-citations: true
# PDF output (pandoc + xelatex), ACM-like single column; see paper-style.tex
documentclass: article
papersize: letter
fontsize: 12pt
geometry: textwidth=5.5in, textheight=8.6in, centering
mainfont: Linux Libertine O
indent: true
linestretch: 1.05
numbersections: false
---

## Introduction

Many mechanisms combine stiff parts that move far with slender parts that deform: robot arms with compliant links, flexure stages, grippers. Their dynamics is described by the equations of
flexible multibody systems [@shabana2013]. For such systems, a floating frame of reference combined with a reduced set of deformation modes offers a good balance between accuracy and cost, as long as the deformation relative to the frame stays small.

\texttt{fmbd} implements this approach in a deliberately compact form. The aim is a code that can be read alongside its theory: each step of the formulation corresponds to a few lines of code, and the code cites the equations of this paper by number, e.g. `% Eq. (4.6)`, with the same symbols. The toolbox builds a system from STL parts, meshed with Gmsh, and provides time simulation, static equilibrium, linearisation, modal analysis and state-space models.

All analyses share one set of equations of motion, a differential-algebraic system (4.6),

$$\mathbf{M}(\mathbf{q})\,\dot{\boldsymbol{\nu}} + \mathbf{G}(\mathbf{q})^\top\boldsymbol{\lambda} = \mathbf{F}(\mathbf{q}, \boldsymbol{\nu}, t), \mathbf{g}(\mathbf{q}) = 0,$$

in which $\mathbf{M}$ is the configuration-dependent mass matrix, $\mathbf{F}$ the vector of applied and velocity-dependent forces, $\mathbf{g}$ the constraints with Jacobian $\mathbf{G}$, and $\boldsymbol{\lambda}$ the Lagrange multipliers. The paper is organised accordingly. Sections 1–6 derive the terms of this
system: the coordinates (Section 1), the contributions of rigid and elastic bodies to $\mathbf{M}$ and $\mathbf{F}$ (Sections 2–4), the constraints (Section 5) and the applied forces (Section 6). Sections 7–9 treat its solution: time integration (Section 7), the initial configuration and static equilibrium (Section 8), and linearisation (Section 9). Equations are numbered by section.

The main ingredients are established: Craig–Bampton reduction [@craig1968], the floating frame formulation [@shabana2013] and generalised-α time integration [@chung1993] in its Lie group form
[@bruls2012]. The contributions of this paper lie in how they are combined and made compact:

- all inertia terms of an elastic unit, including the quadratic velocity forces, written in
  closed form in a small set of invariants (Sections 3.3 and 4);
- rotation constraints of joints as dot products of frame axes, with their Jacobians
  (Section 5.2);
- an initial configuration that satisfies prescribed joint positions exactly and keeps the
  remaining joint coordinates as close as possible to zero (Section 8.1);
- a Newton method for static equilibrium that splits each step by the constraints, so that free
  mechanisms with a singular stiffness are handled (Section 8.2);
- a linearisation that uses the same force evaluation as the simulation, followed by a modal
  state-space model with static correction (Section 9).

**Notation.** Vectors and matrices are bold: lower case for vectors ($\mathbf{v}$,
$\boldsymbol{\omega}$), upper case for matrices ($\mathbf{M}$, $\boldsymbol{\Phi}$). Scalars are plain
($m$, $h$, $\theta$), and so are single components and matrix entries ($d_i$, $M_{s,pq}$).
Subscripts are labels and are not bold. A few vectors keep the capital letter they have in the
code: the forces $\mathbf{F}$, the first moment $\mathbf{S}$, the joint offset $\mathbf{D}$, the node
positions $\mathbf{X}$ and shapes $\mathbf{U}$, and increments $\boldsymbol{\Delta}$.

## 1 Coordinates

The equations of motion are written in the velocities $\boldsymbol{\nu}$, while the configuration $\mathbf{q}$ contains rotations, so we first need both and a way to update one from the other.

**Units.** A system is built from **units**. A rigid body is a unit, and so is a group of elastic bodies tied together. Each unit has a frame with position $\mathbf{r}_i$ and rotation $\mathbf{R}_i$ in the world: the columns of $\mathbf{R}_i$ are the axes of the frame, so $\mathbf{R}_i$ turns a vector from unit axes into world axes. An elastic unit also has modal coordinates $\boldsymbol{\eta}_i$, the amplitudes of its deformation shapes (Section 3). Per unit,

$$\mathbf{q}_i = (\mathbf{r}_i, \mathbf{R}_i, \boldsymbol{\eta}_i), \qquad \boldsymbol{\nu}_i = \begin{bmatrix} \mathbf{v}_i \\ \boldsymbol{\omega}_i \\ \dot{\boldsymbol{\eta}}_i \end{bmatrix}, \tag{1.1}$$

where $\mathbf{v}_i = \dot {\mathbf{r}}_i$ is in world axes and the angular velocity $\boldsymbol{\omega}_i$ is in the axes of the
unit itself. We write $\hat {\mathbf{a}}$ for the skew matrix with $\hat {\mathbf{a}}\, \mathbf{b} = \mathbf{a} \times \mathbf{b}$. A vector $\mathbf{y}$ fixed
in the unit is $\mathbf{R}\mathbf{y}$ in the world and changes as $\mathbf{R}(\boldsymbol{\omega}\times \mathbf{y})$, for any $\mathbf{y}$, so

$$\dot {\mathbf{R}}_i = \mathbf{R}_i\,\hat{\boldsymbol{\omega}}_i.$$

Using $\boldsymbol{\omega}$ in unit axes has a practical reason: the inertia of a body is constant in its own
axes. A rigid unit has its frame at its centre of mass, an elastic unit at its first coupling
(Section 3.1). The system vectors $\mathbf{q}$ and $\boldsymbol{\nu}$ stack those of all units.

**Increments.** $\mathbf{q}$ contains rotation matrices and $\boldsymbol{\nu}$ does not, so $\mathbf{q}$ cannot be updated by
adding. A step $\boldsymbol{\Delta} = [\boldsymbol{\Delta}_v; \boldsymbol{\Delta}_\omega; \boldsymbol{\Delta}_\eta]$, with the same layout as $\boldsymbol{\nu}$, is
applied as

$$\mathbf{q} \oplus \boldsymbol{\Delta} = \left(\mathbf{r} + \boldsymbol{\Delta}_v,\; \mathbf{R} \exp(\hat{\boldsymbol{\Delta}}_\omega),\; \boldsymbol{\eta} + \boldsymbol{\Delta}_\eta\right). \tag{1.2}$$

Positions and modal coordinates are added. The frame is turned by the rotation vector
$\boldsymbol{\Delta}_\omega$: about the axis $\boldsymbol{\Delta}_\omega/|\boldsymbol{\Delta}_\omega|$, in unit axes, over the angle
$|\boldsymbol{\Delta}_\omega|$. Moving with a constant $\boldsymbol{\nu}$ for a time $h$ gives $\mathbf{q}\oplus h\boldsymbol{\nu}$. Every solver
in Sections 7–9 changes $\mathbf{q}$ this way.

The same idea defines derivatives with respect to $\mathbf{q}$. For any function $f(\mathbf{q})$,

$$\frac{\partial f}{\partial \mathbf{q}} := \left.\frac{\partial f(\mathbf{q}\oplus\boldsymbol{\Delta})}{\partial\boldsymbol{\Delta}}\right|_{\boldsymbol{\Delta} = 0},
\qquad\longrightarrow\qquad \dot f = \frac{\partial f}{\partial \mathbf{q}}\,\boldsymbol{\nu}.$$

**Rotation vectors.** A rotation vector $\boldsymbol{\phi}$ stands for a turn about the axis $\boldsymbol{\phi}/\theta$ over
the angle $\theta = |\boldsymbol{\phi}|$. The exponential map (Rodrigues' formula) gives its rotation matrix,
and $\mathbf{T}$ relates the rate of $\boldsymbol{\phi}$ to the angular velocity:

$$\exp(\hat{\boldsymbol{\phi}}) = \mathbf{I} + \frac{\sin\theta}{\theta}\hat{\boldsymbol{\phi}} + \frac{1 - \cos\theta}{\theta^2}\hat{\boldsymbol{\phi}}^2,
\qquad
\mathbf{T}(\boldsymbol{\phi}) = \mathbf{I} - \frac{1 - \cos\theta}{\theta^2}\hat{\boldsymbol{\phi}} + \frac{\theta - \sin\theta}{\theta^3}\hat{\boldsymbol{\phi}}^2. \tag{1.3}$$

If $\mathbf{R} = \exp(\hat{\boldsymbol{\phi}})$, then $\boldsymbol{\omega} = \mathbf{T}(\boldsymbol{\phi})\,\dot{\boldsymbol{\phi}}$. For small $\boldsymbol{\phi}$, $\mathbf{T} \approx \mathbf{I}$ and the two
are the same. They also agree when the axis is fixed, because then $\hat{\boldsymbol{\phi}}\,\dot{\boldsymbol{\phi}} = 0$. They
differ when the axis changes direction during a large rotation. Going the other way,
$\dot{\boldsymbol{\phi}} = \mathbf{T}^{-1}(\boldsymbol{\phi})\,\boldsymbol{\omega}$ with

$$\mathbf{T}^{-1}(\boldsymbol{\phi}) = \mathbf{I} + \tfrac12\hat{\boldsymbol{\phi}} + \left(\frac{1}{\theta^2} - \frac{1 + \cos\theta}{2\theta\sin\theta}\right)\hat{\boldsymbol{\phi}}^2. \tag{1.4}$$

$\log$ is the inverse of $\exp$: it returns the rotation vector of a rotation matrix. We use these
for the rotation of couplings (3.3) and joints (5.2).

**Reference configuration.** This is the configuration the bodies are drawn in: all $\mathbf{R}_i = \mathbf{I}$, all
$\boldsymbol{\eta}_i = 0$ and all joint coordinates are zero.

## 2 Bodies

The mass and stiffness in (4.6) start at the bodies. An elastic body gives finite element matrices,
which Section 3 reduces; a rigid body adds its part of $\mathbf{M}$ and $\mathbf{F}$ directly.

### 2.1 Finite elements
Elastic bodies are meshed with Tet4 or Tet10 elements (Gmsh node order), with linear,
small-strain elasticity. Inside an element the displacement is interpolated from the node
displacements $\mathbf{u}_e$ with the shape functions, $\mathbf{u}(\mathbf{x}) = \mathbf{N}(\mathbf{x})\,\mathbf{u}_e$. The strains are derivatives of
the displacement, so they are also linear in $\mathbf{u}_e$: $\boldsymbol{\varepsilon} = \mathbf{B}\,\mathbf{u}_e$, in Voigt order
$\boldsymbol{\varepsilon} = [\varepsilon_{xx}\ \varepsilon_{yy}\ \varepsilon_{zz}\ \gamma_{yz}\ \gamma_{xz}\
\gamma_{xy}]^\top$. The stresses are $\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$, with $\mathbf{D}$ from the material.

The element matrices follow from the energies. The strain energy is
$\tfrac12\int\boldsymbol{\varepsilon}^\top\boldsymbol{\sigma}\,d\Omega = \tfrac12 \mathbf{u}_e^\top \mathbf{K}_e \mathbf{u}_e$ and the kinetic energy is
$\tfrac12\int\rho\,\dot {\mathbf{u}}^\top\dot {\mathbf{u}}\,d\Omega = \tfrac12\dot {\mathbf{u}}_e^\top \mathbf{M}_e\dot {\mathbf{u}}_e$, with

$$\mathbf{K}_e = \int_{\Omega_e} \mathbf{B}^\top \mathbf{D} \mathbf{B} \, d\Omega, \qquad \mathbf{M}_e = \int_{\Omega_e} \rho \mathbf{N}^\top \mathbf{N} \, d\Omega. \tag{2.1}$$

The integrals use collapsed (Duffy) Gauss rules, exact for straight-sided elements. Adding the
elements gives $\mathbf{K}$ and $\mathbf{M}$ of the body. Each body can have Rayleigh damping $\alpha \mathbf{M} + \beta \mathbf{K}$.

The three directions of a node use the same shape functions, so the mass matrix repeats one
number per node pair: the block for nodes $p$ and $q$ is $M_{s,pq}\,\mathbf{I}$, with
$M_{s,pq} = \int\rho N_p N_q\,d\Omega$. This **scalar mass matrix** $\mathbf{M}_s$ has one row per node. It
turns integrals over the mass into sums over the nodes, $\int\rho\,\mathbf{y}\cdot \mathbf{z}\,d\Omega =
\sum_{pq} M_{s,pq}\,\mathbf{y}_p\cdot \mathbf{z}_q$, which is how all mass properties below are computed.

### 2.2 Rigid bodies
A rigid body is described by its mass $m$, centre of mass $\mathbf{c}$ and inertia $\mathbf{J}$ about $\mathbf{c}$. From a
mesh with node positions $\mathbf{X}_p$ they are

$$m = \sum_{pq} M_{s,pq}, \qquad \mathbf{c} = \frac1m\sum_{pq} M_{s,pq}\,\mathbf{X}_q, \qquad
\mathbf{J} = \operatorname{tr}(\mathbf{Q})\,\mathbf{I} - \mathbf{Q},$$

with $\mathbf{Q} = \sum_{pq} M_{s,pq}\,(\mathbf{X}_p - \mathbf{c})(\mathbf{X}_q - \mathbf{c})^\top.$ They can also come from the closed surface of the body, or be given directly.

With the frame at $\mathbf{c}$, translation and rotation separate into Newton's and Euler's equations,
$m\dot {\mathbf{v}} = \mathbf{f}$ and $\mathbf{J}\dot{\boldsymbol{\omega}} + \boldsymbol{\omega}\times \mathbf{J}\boldsymbol{\omega} = \boldsymbol{\tau}$. $\mathbf{J}$ is constant because $\boldsymbol{\omega}$ is in
body axes. Moving the gyroscopic term to the right gives the constant mass matrix and the inertial
force of a rigid unit:

$$\mathbf{M}_i = \begin{bmatrix} m\mathbf{I} & 0 \\ 0 & \mathbf{J} \end{bmatrix}, \qquad
\mathbf{F}_i^{\mathrm{in}} = \begin{bmatrix} 0 \\ -\boldsymbol{\omega} \times \mathbf{J}\boldsymbol{\omega} \end{bmatrix}. \tag{2.2}$$

## 3 Elastic units

A meshed body has far too many coordinates to put into equation of motion (4.6) directly. An elastic unit therefore uses a moving frame [@shabana2013]: the frame takes the large motion and a few modes $\boldsymbol{\eta}$ take the deformation, which stays small relative to the frame. Below we reduce the mesh to those modes, and compute once the quantities from which Section 4 builds $\mathbf{M}(\mathbf{q})$ and $\mathbf{F}$ at every step.

### 3.1 Reduction
Let $\mathbf{u}$ be the node displacements of all meshes in the unit, relative to the unit frame, with $\mathbf{K}$
and $\mathbf{M}$ the matrices (2.1) of its bodies side by side. The reduction takes two steps.

**Step 1: dependent nodes.** Some nodes cannot move freely:

- **couplings** (RBE2): a node moves rigidly with its coupling point, which has a displacement
  $\mathbf{u}_c$ and a small rotation $\boldsymbol{\theta}_c$: $\mathbf{u} = \mathbf{u}_c + \boldsymbol{\theta}_c \times (\mathbf{x} - \mathbf{x}_c)$;
- **ties**: a slave node follows the master face it lies on, $\mathbf{u}_s = \sum_j N_j(\boldsymbol{\xi})\, \mathbf{u}_{m,j}$, with
  the face shape functions $N_j$ at the projected point $\boldsymbol{\xi}$;
- **node supports**: $\mathbf{u} = 0$.

All three are linear, so they can be written as $\mathbf{u} = \mathbf{T} \mathbf{p}$, with $\mathbf{p}$ the coordinates that remain:
the free nodes and six per coupling. (This $\mathbf{T}$ is not the one of (1.3).) Putting $\mathbf{u} = \mathbf{T}\mathbf{p}$ into
the energies $\tfrac12 \mathbf{u}^\top \mathbf{K} \mathbf{u}$ and $\tfrac12\dot {\mathbf{u}}^\top \mathbf{M}\dot {\mathbf{u}}$ gives the matrices in $\mathbf{p}$:
$\mathbf{K}_r = \mathbf{T}^\top \mathbf{K} \mathbf{T}$ and $\mathbf{M}_r = \mathbf{T}^\top \mathbf{M} \mathbf{T}$.

An elastic interface is a spring instead of a tie, so it removes no coordinates. It acts on the
gap between a slave node and the master face, $\mathbf{e} = \mathbf{u}_s - \sum_j N_j \mathbf{u}_{m,j} = \mathbf{G}_t \mathbf{u}$. The spring
energy $\tfrac12 \mathbf{e}^\top \mathbf{D}_s \mathbf{e}$ adds $\mathbf{K} \leftarrow \mathbf{K} + \mathbf{G}_t^\top \mathbf{D}_s \mathbf{G}_t$, with
$\mathbf{D}_s = A\,[k_n \mathbf{n}\mathbf{n}^\top + k_t(\mathbf{I} - \mathbf{n}\mathbf{n}^\top)]$ for a slave node of area $A$ and normal $\mathbf{n}$: stiffness
$k_n$ along the normal and $k_t$ across it, per unit area.

**Step 2: modes.** The unit frame must be attached to the material somewhere. We clamp the first
coupling and put the frame there: this coupling does not deform, and its motion is the motion of
the frame. (A unit with node supports is clamped at those instead, and its frame is fixed to the
ground.) The other couplings are **interfaces**; their six coordinates each are kept, because
joints attach to them. Split $\mathbf{p}$ into the interface coordinates $\mathbf{p}_b$ and the interior
coordinates $\mathbf{p}_i$. The interior is described with two kinds of shapes [@craig1968]:

- **constraint modes**: move one interface coordinate by one and let the interior settle
  statically. From $\mathbf{K}_{ii}\mathbf{p}_i + \mathbf{K}_{ib}\mathbf{p}_b = 0$ this gives $\mathbf{p}_i = -\mathbf{K}_{ii}^{-1}\mathbf{K}_{ib}\,\mathbf{p}_b$;
- **fixed-interface modes**: the vibration modes $\boldsymbol{\Phi}$ with all interfaces held,
  $\mathbf{K}_{ii}\boldsymbol{\Phi} = \mathbf{M}_{ii}\boldsymbol{\Phi}\,\boldsymbol{\Omega}^2$. The first `Modes` of them are kept.

Together,

$$\mathbf{p} \approx \mathbf{V}\boldsymbol{\eta}, \qquad \mathbf{V} = \begin{bmatrix} \mathbf{I} & 0 \\ -\mathbf{K}_{ii}^{-1}\mathbf{K}_{ib} & \boldsymbol{\Phi} \end{bmatrix},
\qquad (\mathbf{K}_{ii} - \omega^2 \mathbf{M}_{ii})\boldsymbol{\Phi} = 0. \tag{3.1}$$

So $\boldsymbol{\eta}$ holds the six coordinates of each interface, followed by the amplitudes of the kept
modes. The constraint modes make the static response to interface motion exact; the
fixed-interface modes add the dynamics up to their highest frequency.

A load $\mathbf{f}$ on a node set (not on a coupling) would deform the unit in a way these shapes may miss.
Its static shape $\mathbf{K}_{ii}^{-1}\mathbf{f}_i$ is therefore added to $\mathbf{V}$, after removing the part that $\mathbf{V}$
already contains (it is made $\mathbf{M}$-orthogonal to $\mathbf{V}$). The static response to the load is then
exact.

The reduced matrices are $\mathbf{K}_c = \mathbf{V}^\top \mathbf{K}_r \mathbf{V}$, $\mathbf{M}_c = \mathbf{V}^\top \mathbf{M}_r \mathbf{V}$ and $\mathbf{D}_c = \mathbf{V}^\top \mathbf{D}_r \mathbf{V}$, and
the shapes at the nodes are $\boldsymbol{\Phi}_u = \mathbf{T}\mathbf{V}$, so $\mathbf{u} = \boldsymbol{\Phi}_u\boldsymbol{\eta}$.

### 3.2 Kinematics
A node sits at $\mathbf{X}_0$ in the reference configuration, measured from the frame origin. With the
deformation, its position in the unit frame is $\mathbf{u} = \mathbf{X}_0 + \boldsymbol{\Phi}_u\boldsymbol{\eta}$ (from here on $\mathbf{u}$ is this
position, not the displacement). The frame places it in the world, and differentiating with
$\dot {\mathbf{R}} = \mathbf{R}\hat{\boldsymbol{\omega}}$ gives its velocity:

$$\mathbf{x} = \mathbf{r} + \mathbf{R}\,(\mathbf{X}_0 + \boldsymbol{\Phi}_u \boldsymbol{\eta}), \qquad
\dot {\mathbf{x}} = \mathbf{v} + \mathbf{R}\left(\boldsymbol{\omega} \times \mathbf{u} + \boldsymbol{\Phi}_u\dot{\boldsymbol{\eta}}\right). \tag{3.2}$$

The three terms are the velocity of the frame, the velocity from the frame turning, and the
deformation velocity.

A coupling moves the same way. Let $\boldsymbol{\Psi} = [\boldsymbol{\Psi}_t; \boldsymbol{\Psi}_r]$ be its six rows of $\mathbf{V}$, so that
$\boldsymbol{\Psi}_t\boldsymbol{\eta}$ is its displacement and $\boldsymbol{\Psi}_r\boldsymbol{\eta}$ its small rotation in the unit frame. With its
reference position $\mathbf{x}_{c0}$,

$$\mathbf{r}_c = \mathbf{r} + \mathbf{R}\,(\mathbf{x}_{c0} + \boldsymbol{\Psi}_t\boldsymbol{\eta}), \qquad \mathbf{R}_c = \mathbf{R}\exp(\widehat{\boldsymbol{\Psi}_r\boldsymbol{\eta}}). \tag{3.3}$$

For the first coupling $\mathbf{x}_{c0} = 0$ and $\boldsymbol{\Psi} = 0$: it is the frame.

### 3.3 Inertia invariants
The mass matrix and the inertial forces are sums over all nodes, and they change with $\boldsymbol{\eta}$.
Looping over the mesh at every time step would be slow. But these sums are quadratic in the node
positions $\mathbf{u}$, and $\mathbf{u}$ is linear in $\boldsymbol{\eta}$, so they can be built from a few small matrices that
are computed once.

Take $\mathbf{U} = [\mathbf{X}_0, \boldsymbol{\Phi}_u]$ with columns $\mathbf{U}_a$, $a = 0, \dots, n$: column 0 is the reference shape
and column $j$ is mode $j$. Each column lists a 3-vector $\mathbf{U}_{a,p}$ per node $p$. With
$\boldsymbol{\xi} = [1; \boldsymbol{\eta}]$, the node positions are $\mathbf{u}_p = \sum_a \xi_a \mathbf{U}_{a,p}$. The invariants are

$$\mathbf{Q}_{ab} = \mathbf{U}_a^\top \mathbf{M}_s \mathbf{U}_b \in \mathbb{R}^{3\times3}, \qquad \mathbf{s}_a = \mathbf{U}_a^\top \mathbf{M}_s \mathbf{1} \in \mathbb{R}^3, \tag{3.4}$$

that is, $\mathbf{Q}_{ab} = \sum_{pq} M_{s,pq}\,\mathbf{U}_{a,p}\,\mathbf{U}_{b,q}^\top$ and $\mathbf{s}_a = \sum_{pq} M_{s,pq}\,\mathbf{U}_{a,p}$.
$\mathbf{Q}_{ab}$ is the mass-weighted product of shapes $a$ and $b$, and $\mathbf{s}_a$ the mass-weighted sum of
shape $a$. Any sum over the nodes then becomes a small sum over shapes, for example

$$\sum_{pq} M_{s,pq}\,\mathbf{u}_p\,\mathbf{u}_q^\top = \sum_{ab}\xi_a\xi_b\,\mathbf{Q}_{ab}.$$

Some terms need cross products instead of outer products. Those follow from the same matrices:
if $\mathbf{Q} = \sum_{pq} M_{s,pq}\, \mathbf{y}_p \mathbf{z}_q^\top$, then $\sum_{pq} M_{s,pq}\, \mathbf{y}_p \times \mathbf{z}_q = \operatorname{cr}(\mathbf{Q})$
with

$$\operatorname{cr}(\mathbf{Q}) = [Q_{23} - Q_{32};\ Q_{31} - Q_{13};\ Q_{12} - Q_{21}], \tag{3.5}$$

because the first component of $\mathbf{y}\times \mathbf{z}$ is $y_2z_3 - y_3z_2$, and so on.

Two combinations are used often: the first moment $\mathbf{S}$, which is the mass times the position of
the centre of mass in the frame, and the inertia $\mathbf{J}$ about the frame origin:

$$\mathbf{S}(\boldsymbol{\eta}) = \sum_a \xi_a \mathbf{s}_a, \qquad \mathbf{Q}_{uu}(\boldsymbol{\eta}) = \sum_{ab}\xi_a\xi_b \mathbf{Q}_{ab}, \qquad
\mathbf{J}(\boldsymbol{\eta}) = \operatorname{tr}(\mathbf{Q}_{uu})\, \mathbf{I} - \mathbf{Q}_{uu}. \tag{3.6}$$

Both depend on $\boldsymbol{\eta}$: the centre of mass and the inertia shift as the unit deforms.

## 4 Equations of motion of a unit

Here the pieces come together: the kinematics (3.2) and the invariants (3.4) give the mass matrix
and the forces of an elastic unit, and stacking all units gives the equations of motion (4.6).

### 4.1 Mass matrix
By (3.2), the velocity of node $p$ is linear in the unit velocities $\boldsymbol{\nu} = [\mathbf{v}; \boldsymbol{\omega}; \dot{\boldsymbol{\eta}}]$:

$$\dot {\mathbf{x}}_p = \mathbf{L}_p\,\boldsymbol{\nu}, \qquad \mathbf{L}_p = [\,\mathbf{I},\ -\mathbf{R}\hat {\mathbf{u}}_p,\ \mathbf{R}\boldsymbol{\Phi}_p\,].$$

Here $\boldsymbol{\Phi}_p$ are the three rows of $\boldsymbol{\Phi}_u$ that belong to node $p$, and the middle block is
$\mathbf{R}(\boldsymbol{\omega} \times \mathbf{u}_p) = -\mathbf{R}\hat {\mathbf{u}}_p\,\boldsymbol{\omega}$. The kinetic energy is a sum over node pairs, weighted
with the scalar mass matrix:

$$E_k = \tfrac12\sum_{pq} M_{s,pq}\,\dot {\mathbf{x}}_p\cdot\dot {\mathbf{x}}_q = \tfrac12\,\boldsymbol{\nu}^\top \mathbf{M}\,\boldsymbol{\nu},
\qquad \mathbf{M} = \sum_{pq} M_{s,pq}\, \mathbf{L}_p^\top \mathbf{L}_q.$$

Multiplying out $\mathbf{L}_p^\top \mathbf{L}_q$ block by block, every sum over the nodes is one of the invariants of
Section 3.3. The result is the symmetric mass matrix [@shabana2013]

$$\mathbf{M}(\mathbf{q}) = \begin{bmatrix}
m\mathbf{I} & -\mathbf{R}\hat {\mathbf{S}} & \mathbf{R} \mathbf{S}_1 \\
\cdot & \mathbf{J} & \mathbf{M}_{\omega\eta} \\
\cdot & \cdot & \mathbf{M}_c
\end{bmatrix}, \qquad
\mathbf{M}_{\omega\eta,j} = \sum_a \xi_a \operatorname{cr}(\mathbf{Q}_{a,j}), \tag{4.1}$$

with $\mathbf{S}_1 = [\mathbf{s}_1 \dots \mathbf{s}_n]$ and $\mathbf{Q}_{a,j}$ the invariant of $\mathbf{U}_a$ and mode $j$. The diagonal holds
the mass $m$, the inertia $\mathbf{J}$ about the frame origin (3.6) and the modal mass $\mathbf{M}_c$. The
off-diagonal blocks couple the motions: $\mathbf{S}$ is the mass times the position of the centre of mass
in the frame, so $-\mathbf{R}\hat {\mathbf{S}}$ is zero only when the frame sits at the centre of mass; $\mathbf{R} \mathbf{S}_1$ and
$\mathbf{M}_{\omega\eta}$ are the linear and angular momentum of each mode.

### 4.2 Inertial forces
Differentiating (3.2) once more gives the acceleration of node $p$. It has a part with the
accelerations $\dot{\boldsymbol{\nu}}$ and a part that only depends on the velocities:

$$\ddot {\mathbf{x}}_p = \mathbf{L}_p\,\dot{\boldsymbol{\nu}} + \mathbf{R}\,\mathbf{c}_p, \qquad
\mathbf{c}_p = \underbrace{\boldsymbol{\omega}\times(\boldsymbol{\omega}\times \mathbf{u}_p)}_{\text{centripetal}}
    + \underbrace{2\,\boldsymbol{\omega}\times \mathbf{z}_p}_{\text{Coriolis}}, \qquad \mathbf{z}_p = \boldsymbol{\Phi}_p\dot{\boldsymbol{\eta}}.$$

The centripetal term is the acceleration of a point that goes around with the frame. The Coriolis
term appears because the node also moves inside the rotating frame, with the deformation velocity
$\mathbf{z}_p$.

The first part gives $\mathbf{M}\dot{\boldsymbol{\nu}}$. The second part is projected onto the coordinates in the same
way, with $\mathbf{L}_p^\top$, and gives the velocity forces

$$\mathbf{h} = \sum_{pq} M_{s,pq}\, \mathbf{L}_p^\top \mathbf{R}\,\mathbf{c}_q = \begin{bmatrix} \mathbf{h}_v \\ \mathbf{h}_\omega \\ \mathbf{h}_\eta \end{bmatrix},$$

which enter the force vector (4.5) as $-\mathbf{h}$. The three blocks of $\mathbf{L}_p$ give three sums:

$$\mathbf{h}_v = \mathbf{R}\sum_{pq} M_{s,pq}\,\mathbf{c}_q, \qquad
\mathbf{h}_\omega = \sum_{pq} M_{s,pq}\,\mathbf{u}_p\times \mathbf{c}_q, \qquad
h_{\eta,j} = \sum_{pq} M_{s,pq}\,\mathbf{U}_{j,p}\cdot \mathbf{c}_q,$$

where $\mathbf{U}_{j,p}$ is the displacement of node $p$ in mode $j$. Each sum is worked out below, so that
only the invariants remain and no loop over the nodes is needed.

**Force.** Summing $\mathbf{c}_q$ over the nodes replaces $\mathbf{u}_q$ by the first moment $\mathbf{S}$ and $\mathbf{z}_q$ by
$\mathbf{S}_1\dot{\boldsymbol{\eta}}$:

$$\mathbf{h}_v = \mathbf{R}\left[\boldsymbol{\omega}\times(\boldsymbol{\omega}\times \mathbf{S}) + 2\,\boldsymbol{\omega}\times(\mathbf{S}_1\dot{\boldsymbol{\eta}})\right]. \tag{4.2}$$

The first term is the centrifugal force of the whole mass, placed at the centre of mass. It is
zero when the frame is at the centre of mass.

**Moment.** Both terms of $\mathbf{c}_q$ have the form $\boldsymbol{\omega}\times \mathbf{y}_q$, so we need sums of
$\mathbf{u}_p\times(\boldsymbol{\omega}\times \mathbf{y}_q)$. The identity
$\mathbf{u} \times (\boldsymbol{\omega} \times \mathbf{y}) = \boldsymbol{\omega}\,(\mathbf{u} \cdot \mathbf{y}) - \mathbf{y}\,(\mathbf{u} \cdot \boldsymbol{\omega})$ turns these into

$$\sum_{pq} M_{s,pq}\,\mathbf{u}_p\times(\boldsymbol{\omega}\times \mathbf{y}_q) = \boldsymbol{\omega}\operatorname{tr}(\mathbf{W}) - \mathbf{W}^\top\boldsymbol{\omega},
\qquad \mathbf{W} = \sum_{pq} M_{s,pq}\,\mathbf{u}_p\,\mathbf{y}_q^\top,$$

because $\sum M_{s,pq}\,\mathbf{u}_p\cdot \mathbf{y}_q = \operatorname{tr}(\mathbf{W})$ and
$\sum M_{s,pq}\,\mathbf{y}_q\,(\mathbf{u}_p\cdot\boldsymbol{\omega}) = \mathbf{W}^\top\boldsymbol{\omega}$. For the centripetal term $\mathbf{y}_q = \boldsymbol{\omega}\times \mathbf{u}_q$,
which gives $\mathbf{W}_1 = -\mathbf{Q}_{uu}\hat{\boldsymbol{\omega}}$. For the Coriolis term $\mathbf{y}_q = \mathbf{z}_q$, which gives
$\mathbf{W}_2 = \mathbf{Q}_{uz} = \sum_a\sum_k \xi_a\dot\eta_k \mathbf{Q}_{a,k}$. Together,

$$\mathbf{h}_\omega = \boldsymbol{\omega}\operatorname{tr}(\mathbf{W}_1) - \mathbf{W}_1^\top\boldsymbol{\omega} + 2\left[\boldsymbol{\omega}\operatorname{tr}(\mathbf{W}_2) - \mathbf{W}_2^\top\boldsymbol{\omega}\right],
\qquad \mathbf{W}_1 = -\mathbf{Q}_{uu}\hat{\boldsymbol{\omega}}, \quad \mathbf{W}_2 = \mathbf{Q}_{uz}. \tag{4.3}$$

The $\mathbf{W}_1$ part is the familiar gyroscopic moment: $\mathbf{Q}_{uu}$ is symmetric, so $\operatorname{tr}(\mathbf{W}_1) = 0$
and $-\mathbf{W}_1^\top\boldsymbol{\omega} = -\boldsymbol{\omega}\times \mathbf{Q}_{uu}\,\boldsymbol{\omega} = \boldsymbol{\omega}\times \mathbf{J}\boldsymbol{\omega}$, as for a rigid body (2.2),
but with the inertia $\mathbf{J}(\boldsymbol{\eta})$ of the deformed unit. The $\mathbf{W}_2$ part is the moment that arises
because the unit deforms while it turns.

**Modal forces.** For mode $j$ we need $\mathbf{U}_{j,p}\cdot \mathbf{c}_q$. For the centripetal term, use
$\boldsymbol{\omega}\times(\boldsymbol{\omega}\times \mathbf{u}) = \boldsymbol{\omega}\,(\boldsymbol{\omega}\cdot \mathbf{u}) - |\boldsymbol{\omega}|^2 \mathbf{u}$. For the Coriolis term, use
$\mathbf{U}\cdot(\boldsymbol{\omega}\times \mathbf{z}) = -\boldsymbol{\omega}\cdot(\mathbf{U}\times \mathbf{z})$ and (3.5) for the sum of cross products:

$$h_{\eta,j} = \boldsymbol{\omega}^\top \mathbf{Q}_{j,u}\,\boldsymbol{\omega} - |\boldsymbol{\omega}|^2\operatorname{tr}(\mathbf{Q}_{j,u})
- 2\,\boldsymbol{\omega}^\top\operatorname{cr}(\mathbf{Q}_{j,z}), \tag{4.4}$$

with $\mathbf{Q}_{j,u} = \sum_b \xi_b \mathbf{Q}_{j,b}$ and $\mathbf{Q}_{j,z} = \sum_k \dot\eta_k \mathbf{Q}_{j,k}$. The first two terms
are the centrifugal load on mode $j$: a spinning unit is pulled outward and deforms. The last term
is the Coriolis load of the other modes' velocities.

### 4.3 Gravity, elastic and damping forces
The other forces on the nodes are projected with $\mathbf{L}_p^\top$ in the same way. Gravity gives every
node the acceleration $\mathbf{g}$ (world axes). With $\mathbf{g}_b = \mathbf{R}^\top \mathbf{g}$ in the unit frame, the three blocks
of $\mathbf{L}_p$ give

$$\sum_{pq} M_{s,pq}\,\mathbf{L}_p^\top \mathbf{g} = \begin{bmatrix} m\,\mathbf{g} \\ \mathbf{S}\times \mathbf{g}_b \\ \mathbf{S}_1^\top \mathbf{g}_b \end{bmatrix}:$$

the weight, the moment of the weight about the frame origin (zero when the frame is at the centre
of mass), and the load of the weight on each mode. The elastic and damping forces only act on the
modal coordinates, $-\mathbf{K}_c\boldsymbol{\eta} - \mathbf{D}_c\dot{\boldsymbol{\eta}}$. They are internal forces, so they give no net force or
moment on the frame. Collecting everything, with the velocity forces of 4.2,

$$\mathbf{F} = \begin{bmatrix} m\mathbf{g} - \mathbf{h}_v \\ \mathbf{S}\times \mathbf{g}_b - \mathbf{h}_\omega \\ \mathbf{S}_1^\top \mathbf{g}_b - \mathbf{h}_\eta - \mathbf{K}_c\boldsymbol{\eta} - \mathbf{D}_c\dot{\boldsymbol{\eta}} \end{bmatrix}. \tag{4.5}$$

For a rigid unit the same reads $\mathbf{F} = [m\mathbf{g};\ -\boldsymbol{\omega}\times \mathbf{J}\boldsymbol{\omega}]$, by (2.2). Section 6 adds the
forces of springs, dampers and inputs.

A unit on its own would obey $\mathbf{M}\dot{\boldsymbol{\nu}} = \mathbf{F}$. In a system, the joints add the constraints
$\mathbf{g}(\mathbf{q}) = 0$ of Section 5, and with them constraint forces $-\mathbf{G}^\top\boldsymbol{\lambda}$: one unknown multiplier
$\lambda_k$ per constraint, just large enough to keep it satisfied. Stacking all units,

$$\mathbf{M}(\mathbf{q})\,\dot{\boldsymbol{\nu}} + \mathbf{G}(\mathbf{q})^\top\boldsymbol{\lambda} = \mathbf{F}(\mathbf{q}, \boldsymbol{\nu}, t), \qquad \mathbf{g}(\mathbf{q}) = 0. \tag{4.6}$$

$\mathbf{M}$ is block diagonal, one block per unit; the units only meet in the constraints. The unknowns
are $\dot{\boldsymbol{\nu}}$ and $\boldsymbol{\lambda}$. $\mathbf{G} = \partial \mathbf{g} / \partial \mathbf{q}$ is the derivative in the sense of
Section 1, so $\dot {\mathbf{g}} = \mathbf{G}\boldsymbol{\nu}$. (Gravity and the constraints are both called $\mathbf{g}$, as in the code;
the constraints always carry their argument.)

## 5 Constraints

Joints connect the units and take away motion. A constraint is a scalar function $g_k(\mathbf{q})$ that
must stay zero, for example "the distance between these two points along $x$". Its time
derivative is linear in the velocities,

$$\dot {\mathbf{g}} = \mathbf{G}\,\boldsymbol{\nu},$$

so row $k$ of $\mathbf{G}$ tells which combination of velocities constraint $k$ forbids. The force that
enforces it acts along the same row: $\mathbf{G}^\top\boldsymbol{\lambda}$ in (4.6), where $\lambda_k$ is the reaction
force or moment in the locked direction. So for every constraint we need two things, its value and
its row of $\mathbf{G}$. We get them in three steps: the motion of the points that joints attach to (5.1),
the relative motion in a joint (5.2), and from that the locked directions.

### 5.1 Anchors
Joints, transmissions and signals attach to **anchors**: a unit frame or a coupling (3.3). An
anchor has a position $\mathbf{r}_a$ and a rotation $\mathbf{R}_a$. Its velocity and angular velocity (in its own
axes) are linear in $\boldsymbol{\nu}$, which defines the anchor Jacobians:

$$\dot {\mathbf{r}}_a = \mathbf{J}_v\,\boldsymbol{\nu}, \qquad \boldsymbol{\omega}_a = \mathbf{J}_\omega\,\boldsymbol{\nu}, \qquad \mathbf{R}_a^\top\dot {\mathbf{R}}_a = \hat{\boldsymbol{\omega}}_a.$$

For a unit frame this is simple: $\dot {\mathbf{r}}_a = \mathbf{v}$ and $\boldsymbol{\omega}_a = \boldsymbol{\omega}$. A coupling also moves with
the deformation. Write $\mathbf{u}_t = \mathbf{x}_{c0} + \boldsymbol{\Psi}_t\boldsymbol{\eta}$ for its position in the unit frame and
$\boldsymbol{\psi} = \boldsymbol{\Psi}_r\boldsymbol{\eta}$ for its rotation vector, so $\mathbf{r}_a = \mathbf{r} + \mathbf{R} \mathbf{u}_t$ and $\mathbf{R}_a = \mathbf{R}\exp(\hat{\boldsymbol{\psi}})$.
Differentiating both,

$$\dot {\mathbf{r}}_a = \mathbf{v} + \mathbf{R}(\boldsymbol{\omega}\times \mathbf{u}_t) + \mathbf{R}\boldsymbol{\Psi}_t\dot{\boldsymbol{\eta}}, \qquad
\boldsymbol{\omega}_a = \exp(\hat{\boldsymbol{\psi}})^\top\boldsymbol{\omega} + \mathbf{T}(\boldsymbol{\psi})\boldsymbol{\Psi}_r\dot{\boldsymbol{\eta}}.$$

The angular velocity has two parts: that of the unit, turned into the axes of the coupling, and
that of the coupling relative to the unit, by (1.3). Reading off the columns for
$(\mathbf{v}, \boldsymbol{\omega}, \dot{\boldsymbol{\eta}})$:

$$\mathbf{J}_v = [\,\mathbf{I},\ -\mathbf{R}\hat {\mathbf{u}}_t,\ \mathbf{R}\boldsymbol{\Psi}_t\,], \qquad
\mathbf{J}_\omega = [\,0,\ \exp(\hat{\boldsymbol{\psi}})^\top,\ \mathbf{T}(\boldsymbol{\psi})\boldsymbol{\Psi}_r\,]. \tag{5.1}$$

### 5.2 Joints
A joint sits between two anchors, A and B. On each it has a joint frame: a point ($\mathbf{a}$ on A, $\mathbf{b}$ on
B) and an orientation ($\mathbf{Q}_A$, $\mathbf{Q}_B$), both fixed in the anchor. In the world, the two joint
frames are at

$$\mathbf{p}_A = \mathbf{r}_A + \mathbf{R}_A\,\mathbf{a}, \quad \mathbf{R}_A \mathbf{Q}_A \qquad\text{and}\qquad \mathbf{p}_B = \mathbf{r}_B + \mathbf{R}_B\,\mathbf{b}, \quad \mathbf{R}_B \mathbf{Q}_B,$$

and they coincide in the reference configuration.

**Joint coordinates.** The joint is described by where frame B is and how it is turned, as seen
from frame A:

$$\mathbf{d} = \mathbf{Q}_A^\top \mathbf{R}_A^\top (\mathbf{p}_B - \mathbf{p}_A), \quad \mathbf{p}_X = \mathbf{r}_X + \mathbf{R}_X \mathbf{x}, \qquad
\boldsymbol{\phi} = \log \mathbf{R}_{\mathrm{rel}}, \quad \mathbf{R}_{\mathrm{rel}} = \mathbf{Q}_A^\top \mathbf{R}_A^\top \mathbf{R}_B \mathbf{Q}_B. \tag{5.2}$$

$\mathbf{d}$ is the offset between the two points, in the axes of frame A. $\mathbf{R}_{\mathrm{rel}}$ is the
rotation of frame B relative to frame A, and $\boldsymbol{\phi}$ is its rotation vector: the axis times the
angle. Both are zero in the reference configuration. For a hinge about axis 3, for example,
$\boldsymbol{\phi} = [0; 0; \theta]$ with $\theta$ the hinge angle.

**Their rates.** We also need $\dot {\mathbf{d}} = \mathbf{J}_d\,\boldsymbol{\nu}$ and $\dot{\boldsymbol{\phi}} = \mathbf{J}_\phi\,\boldsymbol{\nu}$. A joint point moves
with its anchor, $\dot {\mathbf{p}}_X = \dot {\mathbf{r}}_X + \mathbf{R}_X(\boldsymbol{\omega}_X\times \mathbf{x})$, so with (5.1) it has the Jacobian
$\mathbf{J}_{v,X} - \mathbf{R}_X\hat {\mathbf{x}} \mathbf{J}_{\omega,X}$. For $\mathbf{d}$, both the offset $\mathbf{D} = \mathbf{p}_B - \mathbf{p}_A$ and the axes of A
change:

$$\dot {\mathbf{d}} = \mathbf{Q}_A^\top\big[\underbrace{\mathbf{R}_A^\top(\dot {\mathbf{p}}_B - \dot {\mathbf{p}}_A)}_{\text{points move}}
\;\underbrace{-\;\boldsymbol{\omega}_A\times(\mathbf{R}_A^\top \mathbf{D})}_{\text{axes of A turn}}\big].$$

For $\boldsymbol{\phi}$, the relative angular velocity of B with respect to A, in the axes of frame B, is
$\boldsymbol{\omega}_{\mathrm{rel}} = \mathbf{Q}_B^\top\boldsymbol{\omega}_B - \mathbf{R}_{\mathrm{rel}}^\top \mathbf{Q}_A^\top\boldsymbol{\omega}_A$: the angular velocity
of B minus that of A, with the latter first turned into the axes of B. By (1.3) and (1.4),
$\dot{\boldsymbol{\phi}} = \mathbf{T}^{-1}(\boldsymbol{\phi})\,\boldsymbol{\omega}_{\mathrm{rel}}$. Using $-\boldsymbol{\omega}\times \mathbf{y} = \hat {\mathbf{y}}\,\boldsymbol{\omega}$, the two rates
are

$$\dot {\mathbf{d}} = \mathbf{Q}_A^\top\left[\widehat{\mathbf{R}_A^\top \mathbf{D}}\,\boldsymbol{\omega}_A + \mathbf{R}_A^\top(\dot {\mathbf{p}}_B - \dot {\mathbf{p}}_A)\right], \qquad
\dot{\boldsymbol{\phi}} = \mathbf{T}^{-1}(\boldsymbol{\phi})\left[\mathbf{Q}_B^\top\boldsymbol{\omega}_B - \mathbf{R}_{\mathrm{rel}}^\top \mathbf{Q}_A^\top\boldsymbol{\omega}_A\right]. \tag{5.3}$$

Replacing $\boldsymbol{\omega}_X$ and $\dot {\mathbf{p}}_X$ by their Jacobians times $\boldsymbol{\nu}$ gives $\mathbf{J}_d$ and $\mathbf{J}_\phi$.

**Free and locked directions.** A joint type says which of the six components of $[\mathbf{d}; \boldsymbol{\phi}]$ are
free. The free ones are the joint coordinates $\mathbf{s}$ that springs, dampers, inputs and outputs act on
(Section 6), with the matching rows of $\mathbf{J}_d$ and $\mathbf{J}_\phi$ as their Jacobian $\mathbf{J}_s$. The locked
ones become constraints.

**Translation constraints.** A locked translation is simply $g = d_i = 0$, and its row of $\mathbf{G}$ is
row $i$ of $\mathbf{J}_d$.

**Rotation constraints.** The same could be done with $\phi_i = 0$, but the components of $\boldsymbol{\phi}$
mix when the joint turns far about a free axis, and $\log$ is singular at 180°. Instead we use the
axes of the two joint frames. Let $\mathbf{e}_i^A = \mathbf{R}_A \mathbf{Q}_A \mathbf{e}_i$ be axis $i$ of frame A in the world, and
$\mathbf{e}_j^B = \mathbf{R}_B \mathbf{Q}_B \mathbf{e}_j$ axis $j$ of frame B. A rotation constraint keeps two such axes
perpendicular:

$$g = \mathbf{e}_i^A\cdot \mathbf{e}_j^B = 0.$$

Table 5.1 lists which pairs $(i, j)$ are used. Take a hinge about axis 3. Axis 3 of B must stay
parallel to axis 3 of A, so it must stay perpendicular to axes 1 and 2 of A: the pairs $(1,3)$ and
$(2,3)$. This holds at any hinge angle, because turning about axis 3 does not move axis 3. With
the three translations, the hinge has five constraints and one free coordinate.

**Table 5.1** Axis pairs $(i, j)$ of the rotation constraints.

| free rotations | pairs | meaning |
|---|---|---|
| none | $(2,3), (3,1), (1,2)$ | the frames stay aligned |
| one, about $k$ | $(i, k)$ for $i \neq k$ | axis $k$ of B stays on axis $k$ of A |
| two, $k$ locked | $(i, j)$ with $(i, j, k)$ cyclic | axis $i$ of A stays normal to axis $j$ of B |
| all | none | — |

For the row of $\mathbf{G}$, write $\mathbf{u} = \mathbf{R}_A \mathbf{u}_a$ and $\mathbf{w} = \mathbf{R}_B \mathbf{w}_b$ for the two axes, with $\mathbf{u}_a = \mathbf{Q}_A \mathbf{e}_i$
and $\mathbf{w}_b = \mathbf{Q}_B \mathbf{e}_j$ fixed in the anchors. An axis fixed in an anchor changes as
$\dot {\mathbf{u}} = \mathbf{R}_A(\boldsymbol{\omega}_A\times \mathbf{u}_a)$, so

$$\dot g = \dot {\mathbf{u}}\cdot \mathbf{w} + \mathbf{u}\cdot\dot {\mathbf{w}}
= (\boldsymbol{\omega}_A\times \mathbf{u}_a)\cdot(\mathbf{R}_A^\top \mathbf{w}) + (\mathbf{R}_B^\top \mathbf{u})\cdot(\boldsymbol{\omega}_B\times \mathbf{w}_b)
= -(\mathbf{R}_A^\top \mathbf{w}\times \mathbf{u}_a)\cdot\boldsymbol{\omega}_A - (\mathbf{R}_B^\top \mathbf{u}\times \mathbf{w}_b)\cdot\boldsymbol{\omega}_B,$$

where the last step uses $(\boldsymbol{\omega}\times \mathbf{y})\cdot \mathbf{x} = -(\mathbf{x}\times \mathbf{y})\cdot\boldsymbol{\omega}$. With
$\boldsymbol{\omega}_X = \mathbf{J}_{\omega,X}\,\boldsymbol{\nu}$, the row of $\mathbf{G}$ is

$$\frac{\partial g}{\partial \boldsymbol{\nu}} = -(\mathbf{R}_B^\top \mathbf{u} \times \mathbf{w}_b)^\top \mathbf{J}_{\omega,B} - (\mathbf{R}_A^\top \mathbf{w} \times \mathbf{u}_a)^\top \mathbf{J}_{\omega,A}. \tag{5.4}$$

### 5.3 Transmissions and supports
A transmission couples several coordinates $s_i$ through its stretch

$$\delta = \sum_i c_i s_i, \qquad \dot\delta = \sum_i c_i \mathbf{J}_{s_i}\,\boldsymbol{\nu}.$$

A belt between a motor and a link with ratio $N$, for example, has
$\delta = \theta_{\mathrm{motor}}/N - \theta_{\mathrm{link}}$. The $s_i$ are joint coordinates
(5.2) or anchor coordinates: the displacement and the rotation vector $\log \mathbf{R}_a$ of a reference
point, measured from the reference configuration. A rigid transmission is the constraint
$\delta = 0$, with $\sum_i c_i \mathbf{J}_{s_i}$ as its row of $\mathbf{G}$. An elastic one is a spring on $\delta$
(Section 6) and adds no constraint.

A support works the same way. Supporting a reference point is a joint between it and the ground.
Supporting a joint coordinate is the constraint $s - s_0 = 0$, with row $\mathbf{J}_s$.

All rows together form $\mathbf{g}$ and $\mathbf{G}$. Some may say the same thing twice, for instance when two
hinges share an axis. Such redundant rows are found once, with a pivoted QR of $\mathbf{G}^\top$, and
dropped, so that $\mathbf{G}$ has full rank.

## 6 Applied forces

The rest of $\mathbf{F}$ in (4.6) comes from springs, dampers and inputs. They all act on some coordinate
$s$: a joint coordinate (5.2), a stretch (5.3) or the position of a point. Its rate is
$\dot s = \mathbf{J}_s\boldsymbol{\nu}$. A force $f$ on $s$ has the power $f\,\dot s = (\mathbf{J}_s^\top f)\cdot\boldsymbol{\nu}$, so in the
equations of motion it appears as

$$\mathbf{J}_s^\top f.$$

This one rule gives all applied forces. A spring with stiffness $k$ and a damper with coefficient
$c$ pull back with $f = -(k\,s + c\,\dot s)$:

$$\mathbf{F} \leftarrow \mathbf{F} - \mathbf{J}_s^\top(k\,s + c\,\dot s), \qquad \dot s = \mathbf{J}_s\boldsymbol{\nu}. \tag{6.1}$$

An input $u$ is a force or torque on a coordinate with Jacobian row $\mathbf{J}_u$, and adds $\mathbf{J}_u^\top u$.
A moment $\boldsymbol{\mu}$ in world axes on an anchor has the power $\boldsymbol{\mu}\cdot(\mathbf{R}_a\boldsymbol{\omega}_a)$, with
$\boldsymbol{\omega}_a = \mathbf{J}_\omega\boldsymbol{\nu}$ from (5.1), and adds $\mathbf{J}_\omega^\top \mathbf{R}_a^\top \boldsymbol{\mu}$.

## 7 Time integration

With all parts of (4.6) in place, `simulate` steps it through time while keeping $\mathbf{g}(\mathbf{q}) = 0$. We
use the generalised-α method [@chung1993], in the form for rotations of [@bruls2012].

**The method.** It is easiest to read as Newmark's method with one addition. Newmark advances the
configuration and velocities over a step $h$ with a weighted mean of the old and new
accelerations; the weights are $\beta$ and $\gamma$. Generalised-α does not use the true
accelerations $\dot{\boldsymbol{\nu}}$ for this, but auxiliary accelerations $\mathbf{a}$, a slightly shifted average of
them with the weights $\alpha_m$ and $\alpha_f$. The shift damps vibrations that are too fast for
the step size, which a stiff model always has, while slow motion is hardly affected.

One number sets how much damping: $\rho_\infty$, the factor by which the fastest vibrations
shrink per step (1 is no damping, 0 removes them in one step; option `Rho`). The four parameters
follow from it, chosen so that the method stays second-order accurate:

$$\alpha_m = \frac{2\rho_\infty - 1}{\rho_\infty + 1}, \quad
\alpha_f = \frac{\rho_\infty}{\rho_\infty + 1}, \quad
\gamma = \tfrac12 + \alpha_f - \alpha_m, \quad
\beta = \tfrac14\left(\gamma + \tfrac12\right)^2. \tag{7.1}$$

The auxiliary accelerations are defined by

$$(1 - \alpha_m)\,\mathbf{a}_{n+1} + \alpha_m \mathbf{a}_n = (1 - \alpha_f)\,\dot{\boldsymbol{\nu}}_{n+1} + \alpha_f\dot{\boldsymbol{\nu}}_n, \tag{7.2}$$

and the Newmark update, with the increment (1.2) for the configuration, is

$$\mathbf{q}_{n+1} = \mathbf{q}_n \oplus h\left[\boldsymbol{\nu}_n + h\left(\tfrac12 - \beta\right)\mathbf{a}_n + h\beta\, \mathbf{a}_{n+1}\right], \qquad
\boldsymbol{\nu}_{n+1} = \boldsymbol{\nu}_n + h(1 - \gamma)\mathbf{a}_n + h\gamma\, \mathbf{a}_{n+1}. \tag{7.3}$$

**One step.** The unknowns of a step are the new accelerations $\dot{\boldsymbol{\nu}}_{n+1}$ and multipliers
$\boldsymbol{\lambda}_{n+1}$. Given a guess for $\dot{\boldsymbol{\nu}}_{n+1}$, (7.2) gives $\mathbf{a}_{n+1}$ and (7.3) gives $\mathbf{q}_{n+1}$
and $\boldsymbol{\nu}_{n+1}$. These must satisfy (4.6) at the new time:

$$\mathbf{r} = \mathbf{M}\dot{\boldsymbol{\nu}}_{n+1} + \mathbf{G}^\top\boldsymbol{\lambda}_{n+1} - \mathbf{F} = 0, \qquad \mathbf{g}(\mathbf{q}_{n+1}) = 0.$$

This is solved by Newton iteration, which needs the derivatives of the two residuals with respect
to the unknowns. From (7.2) and (7.3), a change $\Delta\dot{\boldsymbol{\nu}}$ of the new accelerations changes

$$\mathbf{a}_{n+1}\ \text{by}\ c_a\,\Delta\dot{\boldsymbol{\nu}}, \qquad \boldsymbol{\nu}_{n+1}\ \text{by}\ h\gamma c_a\,\Delta\dot{\boldsymbol{\nu}}, \qquad
\mathbf{q}_{n+1}\ \text{by}\ h^2\beta c_a\,\Delta\dot{\boldsymbol{\nu}}, \qquad c_a = \frac{1-\alpha_f}{1-\alpha_m}.$$

With the tangents $\mathbf{K}_t = -\partial \mathbf{F}/\partial \mathbf{q}$ and $\mathbf{C}_t = -\partial \mathbf{F}/\partial\boldsymbol{\nu}$, the force
residual changes by $(\mathbf{M} + h\gamma c_a \mathbf{C}_t + h^2\beta c_a \mathbf{K}_t)\,\Delta\dot{\boldsymbol{\nu}} + \mathbf{G}^\top\Delta\boldsymbol{\lambda}$,
and the constraints by $h^2\beta c_a\,\mathbf{G}\,\Delta\dot{\boldsymbol{\nu}}$. Dividing the constraint equation by
$h^2\beta c_a$ makes the matrix symmetric and keeps it well scaled for small steps. Each iteration
then solves

$$\mathbf{S}\begin{bmatrix}\Delta\dot{\boldsymbol{\nu}} \\ \Delta\boldsymbol{\lambda}\end{bmatrix} = -\begin{bmatrix} \mathbf{r} \\ \mathbf{g}/(h^2\beta c_a)\end{bmatrix},
\qquad \mathbf{S} = \begin{bmatrix} \mathbf{M} + h\gamma c_a \mathbf{C}_t + h^2\beta c_a \mathbf{K}_t & \mathbf{G}^\top \\ \mathbf{G} & 0 \end{bmatrix}. \tag{7.4}$$

$\mathbf{K}_t$ holds the stiffness of the springs (6.1), $\mathbf{J}_s^\top k \mathbf{J}_s$, and the modal stiffness $\mathbf{K}_c$.
$\mathbf{C}_t$ holds the dampers, $\mathbf{J}_s^\top c \mathbf{J}_s$, and the derivatives of the velocity forces of 4.2. Both
are computed once per step and reused in its iterations. That is allowed: the tangents only
affect how fast the iteration converges, not the solution it converges to.

**Start.** The first step needs $\dot{\boldsymbol{\nu}}_0$ and $\boldsymbol{\lambda}_0$. They follow from (4.6) together with
the constraints differentiated twice, $\ddot {\mathbf{g}} = \mathbf{G}\dot{\boldsymbol{\nu}} + \dot {\mathbf{G}}\boldsymbol{\nu} = 0$:

$$\begin{bmatrix} \mathbf{M} & \mathbf{G}^\top \\ \mathbf{G} & 0 \end{bmatrix}\begin{bmatrix}\dot{\boldsymbol{\nu}}_0 \\ \boldsymbol{\lambda}_0\end{bmatrix} = \begin{bmatrix} \mathbf{F} \\ -\dot {\mathbf{G}}\boldsymbol{\nu}\end{bmatrix},
\qquad \mathbf{a}_0 = \dot{\boldsymbol{\nu}}_0.$$

**Velocity projection.** The step enforces the constraints on the positions, $\mathbf{g}(\mathbf{q}_{n+1}) = 0$. The
velocities should then satisfy $\mathbf{G}\boldsymbol{\nu} = 0$, but they do so only up to the integration error. After
each step we therefore correct them with the smallest change $\Delta\boldsymbol{\nu}$, measured in kinetic
energy $\tfrac12\Delta\boldsymbol{\nu}^\top \mathbf{M}\Delta\boldsymbol{\nu}$, that gives $\mathbf{G}(\boldsymbol{\nu} + \Delta\boldsymbol{\nu}) = 0$:

$$\begin{bmatrix} \mathbf{M} & \mathbf{G}^\top \\ \mathbf{G} & 0 \end{bmatrix}\begin{bmatrix}\Delta\boldsymbol{\nu} \\ \boldsymbol{\mu}\end{bmatrix} = \begin{bmatrix} 0 \\ -\mathbf{G}\boldsymbol{\nu}\end{bmatrix}. \tag{7.5}$$

## 8 Initial configuration and equilibrium

Every analysis starts from a configuration that satisfies $\mathbf{g}(\mathbf{q}) = 0$; finding one is the first
task. The second is a state where the system rests, (4.6) with $\boldsymbol{\nu} = \dot{\boldsymbol{\nu}} = 0$, which is the
point to linearise about in Section 9.

### 8.1 Initial configuration
In the reference configuration the joints fit, but the user may ask for other joint positions
(option `Position`), and then the bodies have to be moved until everything fits again. The
equations to solve are

- **hard**: the constraints $\mathbf{g}(\mathbf{q}) = 0$ and the given joint positions $s_i = s_i^*$. Together we
  call them $\mathbf{g}_h(\mathbf{q}) = 0$, with Jacobian $\mathbf{G}_h$;
- **soft**: all other joint coordinates $\mathbf{s}$ should stay as close to zero as the hard equations
  allow.

The soft part matters because the hard equations usually leave freedom. In an arm with only the
elbow angle given, nothing says what the shoulder should do; the soft equations keep it where it
was.

Newton's method linearises the hard equations, $\mathbf{g}_h + \mathbf{G}_h\boldsymbol{\Delta} = 0$. There are fewer equations
than coordinates, so there are many solutions. The smallest one is

$$\boldsymbol{\Delta}_h = -\mathbf{G}_h^+ \mathbf{g}_h,$$

with $^+$ the pseudo-inverse. All other solutions are $\boldsymbol{\Delta} = \boldsymbol{\Delta}_h + \mathbf{N}_h \mathbf{z}$, where the columns
of $\mathbf{N}_h$ span the null space of $\mathbf{G}_h$: the directions that do not change the hard equations. We
use this freedom $\mathbf{z}$ for the soft equations. After the step, the soft coordinates are
$\mathbf{s} + \mathbf{J}_s(\boldsymbol{\Delta}_h + \mathbf{N}_h \mathbf{z})$; making that as small as possible gives
$\mathbf{z} = -(\mathbf{J}_s \mathbf{N}_h)^+(\mathbf{s} + \mathbf{J}_s\boldsymbol{\Delta}_h)$, so

$$\boldsymbol{\Delta} = \boldsymbol{\Delta}_h - \mathbf{N}_h\,(\mathbf{J}_s \mathbf{N}_h)^+\,(\mathbf{s} + \mathbf{J}_s\boldsymbol{\Delta}_h). \tag{8.1}$$

The step is applied with (1.2) and repeated until $\mathbf{g}_h = 0$. Only the frames move; the modal
coordinates stay zero. Elastic transmissions are kept unstretched (as a hard equation), unless all
their coordinates are given. Large positions are approached in several increments, so that each
Newton solve starts close to its solution.

### 8.2 Equilibrium
At rest, $\boldsymbol{\nu} = \dot{\boldsymbol{\nu}} = 0$ and (4.6) reduces to a balance between the forces and the constraint
forces, at constant inputs $\mathbf{u}_0$:

$$\mathbf{r}(\mathbf{q}, \boldsymbol{\lambda}) = \mathbf{F}(\mathbf{q}, 0, \mathbf{u}_0) - \mathbf{G}^\top\boldsymbol{\lambda} = 0, \qquad \mathbf{g}(\mathbf{q}) = 0.$$

The unknowns are $\mathbf{q}$ and $\boldsymbol{\lambda}$. Linearising both equations for a step $(\boldsymbol{\Delta}, \Delta\boldsymbol{\lambda})$,
with the tangent stiffness $\mathbf{K}_t = -\partial \mathbf{r}/\partial \mathbf{q}$,

$$\mathbf{K}_t\,\boldsymbol{\Delta} + \mathbf{G}^\top\Delta\boldsymbol{\lambda} = \mathbf{r}, \qquad \mathbf{G}\,\boldsymbol{\Delta} = -\mathbf{g}.$$

Solving these as one system fails when $\mathbf{K}_t$ is singular, which is common: a free mechanism, such
as a hinge without a spring, has no stiffness. So the step is split in three, each a
least-squares solve:

1. *Correct the constraints.* The smallest step with $\mathbf{G}\boldsymbol{\Delta}_c = -\mathbf{g}$ is $\boldsymbol{\Delta}_c = -\mathbf{G}^+ \mathbf{g}$.
2. *Balance the forces.* The rest of the step must not disturb the constraints, so it lies in the
   null space $\mathbf{N}$ of $\mathbf{G}$: $\boldsymbol{\Delta} = \boldsymbol{\Delta}_c + \mathbf{N} \mathbf{z}$. Multiplying the first equation by $\mathbf{N}^\top$
   removes $\Delta\boldsymbol{\lambda}$, because $\mathbf{N}^\top \mathbf{G}^\top = 0$, and leaves
   $\mathbf{N}^\top \mathbf{K}_t \mathbf{N}\,\mathbf{z} = \mathbf{N}^\top(\mathbf{r} - \mathbf{K}_t\boldsymbol{\Delta}_c)$. With the pseudo-inverse, a direction without
   stiffness gets no step, so free mechanisms stay where they are.
3. *Update the multipliers.* What remains of the first equation, $\mathbf{G}^\top\Delta\boldsymbol{\lambda} = \mathbf{r} - \mathbf{K}_t\boldsymbol{\Delta}$.

In one line,

$$\boldsymbol{\Delta}_c = -\mathbf{G}^+ \mathbf{g}, \qquad
\boldsymbol{\Delta} = \boldsymbol{\Delta}_c + \mathbf{N}\,(\mathbf{N}^\top \mathbf{K}_t \mathbf{N})^+ \mathbf{N}^\top(\mathbf{r} - \mathbf{K}_t\boldsymbol{\Delta}_c), \qquad
\Delta\boldsymbol{\lambda} = (\mathbf{G}^\top)^+(\mathbf{r} - \mathbf{K}_t\boldsymbol{\Delta}). \tag{8.2}$$

The step is shortened (halved) when it does not reduce the residuals, and the iteration stops when
both the force residual and the constraint residual are below the tolerance.

$\mathbf{K}_t$ is computed by finite differences, which is the expensive part. Two things keep it cheap. It
is differenced over the frame coordinates only; for the modal coordinates the known $\mathbf{K}_c$ is used.
And it is reused over several iterations, until progress slows. As in Section 7, an approximate
tangent changes the path of the iteration, not its result. If a free coordinate carries a net
force, no equilibrium exists; the iteration then stalls and fmbd reports that it did not converge.

## 9 Linearisation

For frequencies and control design we want a linear version of (4.6) about an equilibrium
$(\mathbf{q}_0, \boldsymbol{\lambda}_0, \mathbf{u}_0)$ from Section 8.

**Linear equations.** Consider a small motion about the equilibrium: $\mathbf{q} = \mathbf{q}_0\oplus\Delta \mathbf{q}$,
$\boldsymbol{\nu} = \Delta\dot {\mathbf{q}}$, $\boldsymbol{\lambda} = \boldsymbol{\lambda}_0 + \Delta\boldsymbol{\lambda}$ and $\mathbf{u} = \mathbf{u}_0 + \Delta \mathbf{u}$. Putting this
into (4.6) and keeping the first-order terms gives

$$\mathbf{M}\,\Delta\ddot {\mathbf{q}} + \mathbf{C}\,\Delta\dot {\mathbf{q}} + \mathbf{K}\,\Delta \mathbf{q} + \mathbf{G}^\top\Delta\boldsymbol{\lambda} = \mathbf{B}\,\Delta \mathbf{u}, \qquad \mathbf{G}\,\Delta \mathbf{q} = 0,$$

with

$$\mathbf{K} = -\frac{\partial (\mathbf{F} - \mathbf{G}^\top\boldsymbol{\lambda})}{\partial \mathbf{q}}, \qquad
\mathbf{C} = -\frac{\partial \mathbf{F}}{\partial \boldsymbol{\nu}}, \qquad \mathbf{B} = \frac{\partial \mathbf{F}}{\partial \mathbf{u}}. \tag{9.1}$$

$\mathbf{K}$ contains more than the springs. $\mathbf{F}$ also holds gravity, so $\mathbf{K}$ includes its effect (a pendulum
has stiffness without any spring). And $\boldsymbol{\lambda}_0$ is held fixed inside the derivative, so the
term $\mathbf{G}^\top\boldsymbol{\lambda}_0$ adds the stiffness of loaded constraint forces that change direction when
the system moves. The same goes for preloaded springs.

The derivatives are taken by differences on the same force routine that `simulate` uses, so the
linear model and the simulation always agree. $\mathbf{K}$ uses central differences in the increments
(1.2). $\mathbf{F}$ is quadratic in $\boldsymbol{\nu}$ and linear in $\mathbf{u}$, so the differences for $\mathbf{C}$ and $\mathbf{B}$ are exact.

**Removing the constraints.** $\mathbf{G}\,\Delta \mathbf{q} = 0$ says that only some directions are allowed. Let
the columns of $\mathbf{T}$ be an orthonormal basis of these directions (the null space of $\mathbf{G}$, so
$\mathbf{G}\mathbf{T} = 0$) and write $\Delta \mathbf{q} = \mathbf{T} \mathbf{z}$. The independent coordinates $\mathbf{z}$ are as many as the system
has degrees of freedom. Multiplying the linear equations by $\mathbf{T}^\top$ removes the constraint
forces, because $\mathbf{T}^\top \mathbf{G}^\top = 0$:

$$\mathbf{T}^\top \mathbf{M} \mathbf{T}\,\ddot {\mathbf{z}} + \mathbf{T}^\top \mathbf{C} \mathbf{T}\,\dot {\mathbf{z}} + \mathbf{T}^\top \mathbf{K} \mathbf{T}\, \mathbf{z} = \mathbf{T}^\top \mathbf{B} \mathbf{u}, \tag{9.2}$$

where $\mathbf{u}$ now stands for the change $\Delta \mathbf{u}$. This is what `linearize` returns.

**Modes and state space.** From here on $\mathbf{M}$, $\mathbf{C}$, $\mathbf{K}$, $\mathbf{B}$ are the reduced matrices of (9.2). The
undamped modes solve $\mathbf{K} \mathbf{V} = \mathbf{M} \mathbf{V}\,\boldsymbol{\Omega}^2$, with the natural frequencies on the diagonal of
$\boldsymbol{\Omega}$ (`modes`). Normalised so that $\mathbf{V}^\top \mathbf{M} \mathbf{V} = \mathbf{I}$, they also give $\mathbf{V}^\top \mathbf{K} \mathbf{V} = \boldsymbol{\Omega}^2$.
Writing $\mathbf{z} = \mathbf{V}\boldsymbol{\eta}$, with $\boldsymbol{\eta}$ now the mode amplitudes (not those of Section 3), and multiplying
by $\mathbf{V}^\top$:

$$\ddot{\boldsymbol{\eta}} + \left(2\zeta\boldsymbol{\Omega} + \mathbf{V}^\top \mathbf{C} \mathbf{V}\right)\dot{\boldsymbol{\eta}} + \boldsymbol{\Omega}^2\boldsymbol{\eta} = \mathbf{V}^\top \mathbf{B} \mathbf{u}. \tag{9.3}$$

The damping has two parts: that of the dampers in the model, $\mathbf{V}^\top \mathbf{C} \mathbf{V}$, and an optional damping
ratio $\zeta$ per mode. With the state $[\boldsymbol{\eta}; \dot{\boldsymbol{\eta}}]$ this is the state-space model of
`statespace`. Outputs that measure accelerations use (9.3) to replace $\ddot{\boldsymbol{\eta}}$, which gives a
direct term from $\mathbf{u}$ to the output.

Often only the first few modes are kept. The model then responds too little to a constant input,
because the stiffness of the dropped modes is missing. To repair this, the static shapes
$\mathbf{K}^{-1}\mathbf{B}$, the deflection under each input, are added to the kept modes (after making them
$\mathbf{M}$-orthogonal to these), so the static gain stays right.

All of this assumes $(\mathbf{q}_0, \mathbf{u}_0)$ is an equilibrium. If it is not, the constant part of (4.6) does
not cancel and the frequencies mean little, so fmbd warns in that case.

## Concluding remarks

We have presented the formulation of fmbd, from the coordinates of a single unit to a modal
state-space model of a constrained system. Rigid and elastic bodies enter the same equations of
motion (4.6): an elastic unit only adds modal coordinates, and its inertia terms are evaluated
from invariants at a cost independent of the mesh. Constraints, applied forces and all analyses
act on these equations through the same few quantities, the anchor Jacobians and the force
evaluation, which keeps the linear models consistent with the simulation.

The formulation has the limits of its ingredients. The elasticity is linear with small strains,
and the floating frame assumes that the deformation relative to the frame stays small; large
deformations call for a geometrically nonlinear description. The reduced basis captures the
dynamics up to the highest kept fixed-interface frequency, and the static response only to the
interface motion and to the loads whose shapes were added. Springs on the rotation of spherical
joints act on the components of the rotation vector, which deviates from springs on Cardan
angles when such a joint turns far about one axis.

## References

::: {#refs}
:::
