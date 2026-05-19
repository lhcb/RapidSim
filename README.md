# RapidSim [![DOI](https://zenodo.org/badge/50026095.svg)](https://zenodo.org/badge/latestdoi/50026095)

RapidSim is a fast Monte Carlo generator for simulation of heavy-quark hadron decays.
It can be useful for generating background or toy data sets for acceptance studies.

* It uses [`TGenPhaseSpace`][tgenphasespace] to generate b/c quark hadron decays.
* [FONLL][fonll] is used to give the b/c the correct production kinematics for the LHC.
* The daughter particle momenta are smeared correctly using user-defined resolutions.
* Cuts can be made on daughter particle properties (pT, eta, etc).
* Particle mass hypotheses can be swapped to investigate effect of mis-identification.

A more detailed description can be found here https://arxiv.org/abs/1612.07489


## Setup

RapidSim is built with [pixi], which provides ROOT, EvtGen and the necessary
compilers in a self-contained environment on Linux and macOS:

```shell
$ # Install pixi: https://pixi.sh/latest/installation
$ pixi run build
```

To run a decay (validation files live in the source tree):

```shell
$ pixi run -- .build/src/RapidSim.exe validation/Bs2Jpsiphi 10000 1
```
If the max event to save/select is specified and >0, the generation will be stopped early if this number is reached.

For an interactive shell with ROOT, EvtGen and the compiler toolchain on the
PATH:

```shell
$ pixi shell
```

To run the full validation suite:

```shell
$ pixi run test
```

Optionally set `RAPIDSIM_CONFIG` to a directory containing additional
configuration files (smearing, PID, FONLL, particles); this is searched in
addition to the in-tree `config/` and `rootfiles/` directories.

[pixi]: https://pixi.sh

## Decays

To generate a new decay mode you must write a `.decay` file using the following syntax:

* particles are named according to `$RAPIDSIM_ROOT/config/particles.dat` or `$RAPIDSIM_CONFIG/config/particles.dat` and separated by 
  spaces
* a decay is denoted by `->`
* subdecays are demarcated by braces `{}`, e.g.

```
Bs0 -> { Jpsi -> mu+ mu- } { phi -> K+ K- }
```

When RapidSim is first run for a new decay mode it will generate a simple `.config` file.
These files allow properties of the particles and global settings to be configured for each decay.
Run RapidSim with 0 events to generate if you just want to produce the `.config` file.
The default values listed below are those obtained when using the automatically generated `.config` file. 
Note that when lines are removed from the `.config` file this behaviour may change.

## Re-decay

It is possible to specify a command-line option to re-decay a particular heavy-quark hadron.
This means that for a particular hadron, the same parent kinematics are retained but the 
kinematics of the decay products (and their various detector-level smearings) are recomputed.

## Configuration

Global settings should be defined at the start of the file using the syntax:

```
<setting> : <value>
``````

Particle settings should be defined after the corresponding `@#` tag using the same syntax.

### Global settings

* `seed`:
  * Sets the seed for the random number generator
  * Default: 0

* `acceptance`:
  * Sets the type of geometric acceptance to apply to the decay
  * Options are `Any`, `ParentIn`, `AllIn` and `AllDownstream`
  * Default: `Any`

* `geometry`:
  * Sets the detector geometry to apply acceptance for
  * Options are `4pi` and `LHCb`
  * Default: `LHCb`

* `energy`:
  * Sets the pp collision energy (in TeV) used to get the correct parent 
    kinematics
  * Supported: `7`, `8`, `13`, `13.6` and `14`
  * more types may be added in $RAPIDSIM_ROOT/rootfiles/fonll or $RAPIDSIM_CONFIG/rootfiles/fonll

* `parent`:
  * Overrides to determined flavour of the parent particle used to get the 
    correct parent kinematics
  * Supported: `b` and `c`
  * More types may be added in [`rootfiles/fonll`][fonll_dir]

* `ptRange`:
  * Define the pT range to generate the parent in.
  * Syntax: `ptRange : <min> <max>`
  * Defaults: 0–300 (for `4pi` geometry) or 0–100 (for `LHCb` geometry)

* `etaRange`:
  * Define the pseudorapidity range to generate the parent in.
  * Syntax: `etaRange : <min> <max>`
  * Defaults: -8–8 for 4pi geometry and 1–6 for LHCb geometry

* `minWidth`:
  * Sets the minimum resonance width (in GeV) to be generated
  * Narrower resonances will be generated with a fixed mass
  * Default: 0.001 GeV

* `maxAttempts`:
  * Sets the maximum number of attempts allowed to generate each decay
  * Default: 1000

* `paramsStable`:
  * Defines the set of parameters to be added to the histograms/tree
    for each stable particle in the decay.
  * Syntax: `paramsStable : <params>`, where `<params>` is a space separated 
    list of the types defined in [the Parameters section](#parameters).
  * Example: `paramsStable : P PT`

* `paramsDecaying`:
  * Defines the set of parameters to be added to the histograms/tree
    for each decaying particle in the decay.
  * Syntax is `paramsDecaying : <params>` where `<params>` is a space separated 
    list of the types defined in [the Parameters section](#parameters).
  * Example: `paramsDecaying : M P PT`

* `paramsTwoBody`:
  * Defines the set of parameters to be added to the histograms/tree
    for each two-body combination of particles in a 3+-body decay.
  * Syntax is `paramsTwoBody : <params>` where `<params>` is a space separated 
    list of the types defined in [the Parameters section](#parameters).
  * Example: `paramsTwoBody : M M2`

* `paramsThreeBody`:
  * Defines the set of parameters to be added to the histograms/tree
    for each three-body combination of particles in a 4+-body decay.
  * Syntax is `paramsThreeBody : <params>` where `<params>` is a space 
    separated list of the types defined in [the Parameters 
    section](#parameters).
  * Example: `paramsThreeBody : M M2`

* `param`:
  * Defines a new parameter to be added to the histograms/tree
  * Syntax is `param : <name> <type> <particles> [TRUE]`, where:
    * `<type>` must be one of the types defined in [the Parameters section](#parameters),
    * `<particles>` is a space separated list of particle indices
    * The optional `TRUE` argument means the parameter will be calculated using 
      the true unsmeared momenta of the particles
  * Example: `param : mSq12 M2 1 2 TRUE`

* `cut`:
  * Applies a cut on a given parameter
  * Syntax is `cut : <param> <type> <min/max> [<max>]`, where
    * `<param>` is the name of a parameter (must be defined using `param`)
    * `<type>` is one of "min", "max", "range" or "veto"
    * `<min>` and/or `<max>` define(s) the cut value(s)

* `shape`:
  * Sets a 1D or 2D PDF to generate events according to
  * Syntax is `shape : <file> <hist> <paramX> [<paramY>]`, where:
    * `<file>` is the file containing the histogram (path must be absolute or relative to run directory),
    * `<hist>` is the name of the histogram,
    * `<paramX>` is the name of the parameter on the X-axis (must be defined 
      using `param`)
  * For 2D histograms `<paramY>` is the name of the parameter on the Y-axis
  * RapidSim supports TH1F, TH1D, TH2F and TH2D histograms,
  * The dimensionality of the histogram will be inferred from the number of 
    parameters given
  * Defaults to phase-space distribution

* `useEvtGen` :
  * Perform decays using the external EvtGen generator
  * Syntax is `useEvtGen : TRUE`
  * Note any value for this parameter will turn EvtGen ON (even FALSE)
  * To turn EvtGen OFF do not define this parameter

* `evtGenUsePHOTOS` :
  * If EvtGen is used to generate decays then PHOTOS will be used for final state radiation
  * Syntax is `evtGenUsePHOTOS : TRUE`
  * Note any value for this parameter will turn PHOTOS ON (even FALSE)
  * To turn PHOTOS OFF do not define this parameter
  * Note also that if EvtGen is not in use then this parameter will have no effect

* `pid` :
  * Sets the particle identification scheme to use
  * Syntax is `pid : <scheme>`, where 
    * `<scheme>` is the name of the file that defines the scheme (default `LHCbGenericPID`)
  * More types may be defined in $RAPIDSIM_ROOT/config/pid or $RAPIDSIM_CONFIG/config/pid

### Particle settings

* `name`:
  * A user-friendly name for this particle to be used in variable names
  * Default: Automatically generated unique name of particle

* `smear`:
  * The type of momentum or IP smearing to apply to this particle
    * May be used multiple times for the same particle to define different types of smearing
  * Supported momentum types: `LHCbGeneric`, `LHCbElectron`, `AtlasMuon`, or `AtlasHadron`
  * Supported IP types: `LHCbGenericIP`
    * More types may be defined in $RAPIDSIM_ROOT/config/smear or $RAPIDSIM_CONFIG/config/smear
  * Default: `LHCbElectron` (for electrons/positrons), otherwise `LHCbGeneric`

* `invisible`:
  * Whether the particle should be treated as invisible
  * Invisible particles are not included when determining non-truth parameters
  * Invisible particles are also ignored when applying the `AllIn` or `AllDownstream` 
    acceptance requirement for the default LHCb geometry
  * A corrected mass parameter will be added automatically if any particles are invisible
  * Default: `true` (for neutrinos), otherwise `false`

* `altMass`:
  * Adds alternative mass hypotheses for this particle
  * Syntax is `altMass : <particles>` where `<particles>` is a space separated 
    list of particle type names as listed in [`config/particles.dat`][particles_file]
  * If alternative hypotheses are defined for multiple particles then 
    parameters will be calculated for all combinations of one or two 
    mis-identifications

* `evtGenModel` :
  * Sets the decay model to be used if this particle is decayed using EvtGen
  * Syntax is `evtGenModel : <model> [<params>]` where `<model>` is the name of 
    the model as used in DECAY.DEC and `<params>` is a space separated list of
    parameters to be passed to the model
  * Default: PHSP

## Parameters

* `M`: The invariant mass of the combination of the given particles
* `M2`: The squared invariant mass of the combination of the given particles
* `MT`: The transverse mass of the combination of the given particles
* `E`: The energy of the combination of the given particles
* `ET`: The transverse energy of the combination of the given particles
* `P`: The total momentum of the combination of the given particles
* `PX`: The X momentum of the combination of the given particles
* `PY`: The Y momentum of the combination of the given particles
* `PZ`: The Z momentum of the combination of the given particles
* `PT`: The transverse momentum of the combination of the given particles
* `vtxX`: The X coordinate of the decay vertex of the given particles
* `vtxY`: The Y coordinate of the decay vertex of the given particles
* `vtxZ`: The Z coordinate of the decay vertex of the given particles
* `origX`: The X coordinate of the origin vertex of the given particles
* `origY`: The Y coordinate of the origin vertex of the given particles
* `origZ`: The Z coordinate of the origin vertex of the given particles
* `IP` : Impact parameter to own primary vertex
* `SIGMAIP` : Error on impact parameter to own primary vertex
* `MINIP` : Minimum impact parameter to any primary vertex
* `SIGMAMINIP` : Error on minimum impact parameter to a primary vertex
* `FD` : Flight distance
* `eta`: The pseudorapidity of the combination
* `phi`: The azimuthal angle of the combination
* `y`: The rapidity of the combination
* `gamma`: The relativistic gamma factor of the combination
* `beta`: The velocity of the combination
* `theta`: The angle between the first two particles in the rest frame of 
  the combination of the remaining particles
  e.g. 1 2 1 3 would give the angle between 1 and 2 in the rest frame 
  of a resonance in `m_13` (i.e. `theta_13` in a 3-body decay)
* `costheta`: The cosine of theta
* `Mcorr`: The corrected mass of the combination of the given particles
  correcting for any invisible particles
* `ProbNNmu` : Particle ID feature for the muon hypothesis
* `ProbNNe` : Particle ID feature for the electron hypothesis
* `ProbNNpi` : Particle ID feature for the pion hypothesis
* `ProbNNk` : Particle ID feature for the kaon hypothesis
* `ProbNNp` : Particle ID feature for the proton hypothesis

## External Generators

### EvtGen

EvtGen may be used to generate decays allowing for non-phasespace decay models.

* The pixi environment provides EvtGen (with Pythia, PHOTOS and Tauola
  externals) automatically — no manual setup is required.
* Enable it for a specific decay by setting `useEvtGen : TRUE` in the
  `.config` file (see the global settings above). PHOTOS final-state radiation
  is enabled with `evtGenUsePHOTOS : TRUE`.
* At present, all identical particles must decay into the same mode.

## TODO

* Add more particles to [`particles.dat`][particles_file]
* Add derived class for particle specific parameters to simplify RapidParam interface
* Combine `.decay` and `.config` files to simplify interface

## Authors

 * Greig Cowan 2015
 * Dan Craik 2016
 * Matt Needham 2015

## Additional contributions

 * Léa Dreyfus 2025

[tgenphasespace]: https://root.cern.ch/doc/master/classTGenPhaseSpace.html
[fonll]: http://cacciari.web.cern.ch/cacciari/fonll/
[fonll_dir]: rootfiles/fonll
[particles_file]: config/particles.dat
