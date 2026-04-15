# libcloudph++ documentation

🌦️ Welcome to the libcloudph++ documentation. This is a work in progress.

Libcloudph++ is described in detail in the following peer-reviewed publications:
- [libcloudph++ 1.0: a single-moment bulk, double-moment bulk, and particle-based warm-rain microphysics library in C++](https://gmd.copernicus.org/articles/8/1677/2015/)
- [libcloudph++ 2.0: aqueous-phase chemistry extension of the particle-based cloud microphysics scheme](https://gmd.copernicus.org/articles/11/3623/2018/)


---
### Microphysical Schemes:


####  💧 **Lagrangian Scheme**
Particle-based microphysics using the Super-Droplet Method, based on Shima et al. (2009) and Shima et al. (2020).
- **Available processes**
  - Detailed description of droplet growth / evaporation
  - Aerosol processes
  - Collision-coalescence
  - Advection and sedimentation
  - Ice nucleation (singular or time-dependent formulation) / melting
  - Depositional growth of ice
- **Type**: Compiled library (`src/`, `include/libcloudph++/lagrangian/`)

####  ☁️ **1-Moment Bulk Scheme**
Single-moment bulk scheme based on Kessler (1995), with ice microphysics parametrization based on Grabowski (1999).
Only the total mass of water per category (cloud / rain / iceA / iceB) is considered, in addition to heat and moisture content.

- **Available processes**:
    - Condensation and evaporation
    - Autoconversion and collection 
    - Sedimentation
    - Ice processes including homogeneous and heterogeneous nucleation, growth by deposition and riming
- **Type**: Header-only library (`include/libcloudph++/blk_1m/`)

####  🌧️ **2-Moment Bulk Scheme**
Double-moment scheme based on Morrison and Grabowski (2007). 
  Condensed water is divided into two categories:
  cloud water and rain water. In addition to the total mass of
  water in both categories, concentrations
  are also predicted. 
- **Available processes**:
    - Condensation and evaporation
    - Autoconversion and collection
    - Sedimentation
    - Ice processess not implemented (yet)
- **Type**: Header-only library (`include/libcloudph++/blk_2m/`)
---


###  Library Components

#### **⚙️ Shared Utilities** (`include/libcloudph++/common/`)
- Equations for thermodynamics, droplet growth, terminal velocity, etc.
- Physical constants
- Mathematical utilities, numerical methods and helpers

#### ☁️ Use cases (`models/`)

- **Kinematic 2D model** - A simplified two-dimensional atmospheric model implementation for testing and demonstrating microphysical schemes


#### 🔨 Build System

The project uses **CMake** as its build system:
- `CMakeLists.txt`: Main build configuration
- `cmake/`: Additional CMake modules and utilities
- Header-only schemes require no compilation
- Lagrangian scheme produces linkable libraries

#### 🔗 Language Bindings (`bindings/`)

Interfaces for **Python** with NumPy.

#### ✅ Testing (`tests/`)

- **Unit tests**: individual component validation
- **Benchmark cases**: performance and accuracy validation
- **Inter-scheme comparisons**

---