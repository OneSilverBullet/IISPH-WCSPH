# 🌊 **IISPH-WCSPH: SPH Fluid Simulation Framework**

### 📜 **0. Project Motivation**

Welcome to the **IISPH-WCSPH** project! This project serves as a **basic framework** for the **SPH** (Smoothed Particle Hydrodynamics) field, created to fill the gap of limited open-source SPH simulation resources. The goal is to provide a reliable platform for those in the SPH domain to build upon and experiment with.

🔍 **Key Objectives:**

1. **Lab Requirement:** Build a solid SPH fluid simulation platform as per the lab's needs.
2. **Open Source:** Share this framework with the broader SPH community, hoping to assist researchers and developers in fluid simulation.

### 📚 **1. Project Details**

- **Programming Language:** C++
- **Framework:** Windows
- **Libraries Used:** QT
- **Development Environment:** Visual Studio 2017
- **Note:** This framework mainly uses **C++** (about 90%), with minimal code related to **QT** or rendering. It leverages **QT's OpenGL API** for visualization. If you're familiar with C++, you'll find the code easy to read and extend.

### ⚙️ **2. Key Features**

This project provides a comprehensive set of features to work with SPH simulations:

- **Basic Math Library:** Includes essentials like **MATRIX**, **VECTOR**, and more for fluid computations.
- **Spatial Partitioning Framework:** The backbone for implementing particle-based fluid simulations in SPH.
- **WCSPH Fluid Simulation:** One of the classic methods used for fluid dynamics in SPH.
- **IISPH Fluid Simulation:** Another classical SPH method for simulating fluid motion.
- **Distance Field-Based Boundary Detection:** Foundation for **solid-liquid coupling** in fluid simulation.
- **Marching Cubes Algorithm:** Converts particle-based fluid data into a **meshed surface**, allowing for realistic visualization.
- **QT Visualization Interface:** A user-friendly interface for visualizing the fluid simulation results, making it easier to modify and experiment with different methods.

### 🔥 **3. Current Progress**

Here's what the framework currently offers:

- Fluid simulation with realistic water dynamics 🏞️💧
- Visual feedback for different fluid behaviors 🎮
- Easy to integrate and modify methods for further research and experimentation 🔧


![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/1.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/2.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/3.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/4.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/5.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/6.png)
