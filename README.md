# DEEP (District Energy & Environmental Planning) Framework

This project provides high-performance Fortran tools for the analysis and simulation of building energy systems at district scale, including thermal networks and distributed renewables. It is intended for both the academic community and professional practitioners (engineers, energy consultants, and building designers).

---

## Project Goals

- provide **scientifically sound** numerical tools
- ensure **high computational efficiency**
- support **reuse and integration** in professional workflows
- maintain a **clean, well-documented, and verifiable** code base

---

## Project Status

- ✔ stable core - under active development
- ✔ supported platforms: Linux / Windows  
- ✔ tested compilers: Intel ifx/ifort

---

## Requirements

- Fortran 2008 (or later)
- Third-party dependencies (see below)

---

## Third-Party Software

This project uses the following open-source third-party software:

- **json-fortran**  
  License: BSD 3-Clause  
  Copyright (c) Jacob Williams and contributors  
  https://github.com/jacobwilliams/json-fortran

The corresponding license texts are included in the respective directories or in the project distribution.

---

## License

This project is released under the **Apache License, Version 2.0**.

- free use in academic and professional contexts
- commercial use permitted
- no warranty is provided
- no liability is assumed for the use of the results

See the `LICENSE` file for the full license text.

---

## Disclaimer

This software is provided **“as is”**, without warranty of any kind, express or implied.  
The author(s) assume **no responsibility** for the use of this software in engineering design, regulatory compliance, or contractual applications.

Users are responsible for validating the correctness and suitability of the software for their intended purpose.

---

## 📖 How to Cite

If you use DEEP, please cite the software as:

DEEP Development Team (2026).  
DEEP – District Energy & Environmental Planning Framework.  
GitHub repository: https://github.com/marcelloa/deep-district-energy

## 📚 Publications

The following publications describe the methodology and applications of DEEP:

1. Aprile, M., Villa, G., Scoccia, R., Angelotti, A. (2026). A Tool for the Optimal Design of Climate Resilient Thermal Energy Networks. Part 1: the District Building Energy Performance Module. In: Proceedings of the 15th REHVA HVAC World Congress - CLIMA 2025
   
2. Villa, G., Angelotti, A., Aprile, M. (2026). A Tool for the Optimal Design of Climate Resilient Thermal Energy Networks. Part 2: District Heating and Cooling Network Model. In: Proceedings of the 15th REHVA HVAC World Congress - CLIMA 2025. Lecture Notes in Civil Engineering, vol 762. Springer, Cham. https://doi.org/10.1007/978-3-032-06806-4_30

3. Angelotti, A., Villa, G., Aprile, M. (2026). A Tool for the Optimal Design of Climate Resilient Thermal Energy Networks. Part 3: the Ground Loop Model. In: Proceedings of the 15th REHVA HVAC World Congress - CLIMA 2025.

---

## Contributing

Contributions, bug reports, and suggestions are welcome.

- open an **issue** to report problems or request features
- submit a **pull request** for code contributions
- please follow the coding style and contribution guidelines of the project

---

## 👥 Project Team

### Project Lead
**Marcello Aprile**  
Politecnico di Milano  
Role: Concept, Architecture, Core Development, Building & HVAC module

### Core Contributors
**Adriana Angelotti**  
Politecnico di Milano  
Role: Ground loop module

**Giorgio Villa**  
Politecnico di Milano  
Role: Hydraulic network module


## Contact

For questions, collaborations, or academic inquiries:

**Marcello Aprile**  
marcello.aprile@polimi.it

---

## Acknowledgements

The present research was carried out within the Project “Network 4 Energy Sustainable Transition – NEST” funded by the National Recovery and Resilience Plan (NRRP), Mission 4 Component 2 Investment 1.3 - Call for tender No. 1561 of 11.10.2022 of Ministero dell’Università e della Ricerca (MUR); funded by the European Union – NextGenera-tionEU. Project code PE0000021, Concession Decree No. 1561 of 11.10.2022 adopted by Ministero dell’Università e della Ricerca (MUR), CUP - D43C22003090001, according to attachment E of Decree No. 1561/2022.
