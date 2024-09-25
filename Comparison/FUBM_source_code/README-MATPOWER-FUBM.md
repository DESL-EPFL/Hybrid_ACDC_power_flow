![MATPOWER-FUBM][logo-fubm]

A Power System Simulation Package for MATLAB
--------------------------------------------

- **MATPOWER Website**          - https://matpower.org
- **MATPOWER GitHub Project**   - https://github.com/MATPOWER/matpower
- **MATPOWER-FUBM GitHub Project**   - https://github.com/AbrahamAlvarezB/matpower-fubm

MATPOWER is a package of M-files for solving power flow, continuation
power flow and optimal power flow problems using MATLAB or Octave. It
is intended as a simulation tool for researchers and educators that is
easy to use and modify. MATPOWER is designed to give the best
performance possible while keeping the code simple to understand and
modify.

MATPOWER-FUBM is an extended and modified version of MATPOWER. It combines
The Flexible Universal Branch Model (FUBM) with the MATPOWER tool to solve
power flow and optimal power flow problems for AC, DC and AC/DC power grids.
The FUBM formulation provides a direct link between the AC and DC parts of 
the grid allowing for solving the entire network within a unified frame of
reference (not sequentially). It is also capable of realistically model any
element within the AC/DC powergrid, ranging from conventional AC transmission
lines to multiple types of AC/DC interface devices such as Voltage Source 
Converters (VSC). The FUBM model introduces optional control variables for 
voltage, power, and voltage droop control. 

MATPOWER and MATPOWER-FUBM are intended as a simulation tool for researchers
and educators that is easy to use and modify. Both of them are designed to 
give the best performance possible while keeping the code simple to understand 
and modify.

The FUBM does not create any conflicts with the traditional usage of MATPOWER.
Even though the MATPOWER-FUBM tool is designed to be fully compatible with 
MATPOWER, it has only been fully tested for Power Flow and Optimal Power Flow.
The MATPOWER-FUBM tool is not compatible with Octave. However, this feature will
be added in the future.

MATPOWER releases can be downloaded from the [MATPOWER website][1],
and the latest stable and work-in-progress versions can always be
downloaded or cloned from the [MATPOWER GitHub project][2]. The
`master` branch should always contain a stable version.

MATPOWER-FUBM releases can be downloaded from the [MATPOWER-FUBM 
GitHub project][40].

System Requirements
-------------------
*   [MATLAB][3] version 7.3 (R2006b) or later, or
*   [GNU Octave][4] version 4 or later (Only MATPOWER)


Getting MATPOWER-FUBM
---------------------
#### Current Development Version

There are two options for obtaining the most recent development version
of MATPOWER-FUBM from the `master` branch on GitHub.
**Note:** This does _not_ include the [MATPOWER Extras][7d].

1. Clone the [MATPOWER-FUBM repository from GitHub][40].
   *Use this option if you want to be able to easily update to the current
   development release, with the latest bug fixes and new features, using a
   simple `git pull` command, or if you want to help with testing or
   or development. This requires that you have a [Git client][5] (GUI
   or command-line) installed.*
    - From the command line:
        - `git clone https://github.com/AbrahamAlvarezB/matpower-fubm.git`
    - Or, from the [MATPOWER-FUBM GitHub repository page][40]:
        - Click the green **Clone or download** button, then **Open in Desktop**.

2. Download a ZIP file of the MATPOWER-FUBM repository from GitHub.
   *Use this option if you need features or fixes introduced since
   the latest versioned release, but you do not have access to or
   are not ready to begin using Git (but don't be afraid to
   [give Git a try][6]).*
    - Go to the [MATPOWER-FUBM GitHub repository page][40].
    - Click the green **Clone or download** button, then **Download ZIP**.

#### Contributing to MATPOWER-FUBM
**For Contributing to MATPOWER-FUBM**
Please contact:
1. Abraham Alvarez-Bustos at abraham.alvarez-bustos@durham.ac.uk 
                          or snoop_and@hotmail.com
2. Behzad Kazemtabrizi at behzad.kazemtabrizi@durham.ac.uk
3. Mahmoud Shahbazi at mahmoud.shahbazi@durham.ac.uk


Installation
------------

Installation and use of MATPOWER-FUBM requires familiarity with the
basic operation of MATLAB. Make sure you follow the installation 
instructions for the version of MATPOWER you are installing. The 
installation process for MATPOWER-FUBM has been included to the 
original installation script developed by MATPOWER.

1.  **Get a copy of MATPOWER-FUBM** as described above. Clone the 
    repository or download and extract the ZIP file of the MATPOWER
    -FUBM distribution and place the resulting directory in the 
    location of your choice and call it anything you like. We will 
    use `<MATPOWER-FUBM>` as a placeholder to denote the path to this
    directory (the one containing `install_matpower.m`). The files in
    `<MATPOWER-FUBM>` should not need to be modified, so it is 
    recommended that they be kept separate from your own code.

2.  **Run the installer.**
    - Open MATLAB and change to the `<MATPOWER-FUBM>` directory.
    - Run the installer and follow the directions to add the
      required directories to your MATLAB path, by typing:

            install_matpower

3.  **That's it.** There is no step 3.
    - But, if you chose not to have the installer run the test suite for
      you in step 2, you can run it now to verify that MATPOWER is
      installed and functioning properly, by typing:

            test_matpower


Running MATPOWER-FUBM
---------------------
The FUBM requires extra data to identify AC and DC branches, 
Transformers, Phase Shifter Transformers (PST), Controlled 
Tap Transformers (CTT), Voltage Source Converters (VSC), 
STATCOMS, and their controls. This data is added as extra 
columns in the branch matrix of each case in the data folder.
The MATPOWER-FUBM code includes a **Quick guide** on the **docs**
folder, for details on how to simulate different elements and
their voltage and power controls.

To run a controlled AC/DC Newton Power Flow on the modified 30-bus
system specified in the file `fubm_case_30_2MTDC_ctrls_vt2_pf_dp.m`,
with the default algorithm options, at the MATLAB prompt, type:

```matlab
runpf('fubm_case_30_2MTDC_ctrls_vt2_pf_dp')
```

To run a controlled AC/DC Optimal Power Flow on the modified 30-bus
system specified in the file `fubm_case_30_2MTDC_ctrls_vt2_pf_dp.m`,
The algorithm options should be adjusted as follows:
At the MATLAB prompt, type:

```matlab
mpopt = mpoption('opf.ac.solver', 'MIPS', 'mips.max_it',5000,'opf.violation',1e-6);
runopf('fubm_case_30_2MTDC_ctrls_vt2_pf_dp',mpopt)
```

By default, the results of the simulation are pretty-printed to the
screen, but the solution can also be optionally returned in a `results`
struct. The following example shows how simple it is, after running a DC
OPF on the 118-bus system in `case118.m`, to access the final objective
function value, the real power output of generator 6 and the power flow
in branch 51.

```matlab
results = rundcopf('case118');
final_objective = results.f;
gen6_output     = results.gen(6, PG);
branch51_flow   = results.branch(51, PF);
```

It works similarly for the FUBM extra variables. For example, running
the opf using KNITRO on the case 'fubm_case_30_2MTDC_ctrls_vt2_pf', 
the optimal theta_sh variable controlling an active power flow of 7.5MW
on a PST can be accessed as:

```matlab
PF = 14; SHIFT = 10;
mpopt = mpoption('opf.ac.solver', 'KNITRO','knitro.tol_x',1e-10,'knitro.tol_f',1e-4,'opf.violation',1e-6);
results = runopf('fubm_case_30_2MTDC_ctrls_vt2_pf_dp',mpopt)
final_objective = results.f;
PST_P_flow = results.branch(25, PF);
PST_Theta_sh = results.branch(25, SHIFT);
```

For additional MATPOWER info, see the [MATPOWER User's Manual][8], the
[on-line function reference][9], or the built-in help documentation for
the various MATPOWER and MATPOWER-FUBM functions. For example:

    help runpf
    help runopf
    help mpoption
    help caseformat
    help bustypes_fubm


Documentation
-------------

There are five primary sources of documentation for MATPOWER-FUBM.
  - [MATPOWER User's Manual][8]
  - [MOST User's Manual][10]
  - [MATPOWER Online Function Reference][9]
  - [MATPOWER-FUBM Quick Guide][41]
  - MATLAB's `help` command

#### Manuals

The MATPOWER and MOST User's Manuals are included in the distribution
([`docs/MATPOWER-manual.pdf`][8] and [`most/docs/MOST-manual.pdf`][10]) and
the latest released versions are always available online, respectively, at:
  - https://matpower.org/MATPOWER-manual.pdf
  - https://matpower.org/MOST-manual.pdf.

Previous versions are also available at
  - https://matpower.org/doc/manuals

#### Built-in Help

Each M-file has its own documentation which can be accessed by typing at
the MATLAB prompt:

    help <name of M-file>

Documentation for the case data file format can be found by typing:

    help caseformat

If something is still unclear after checking the manual and the help,
the source code *is* the documentation. :wink:

#### Changes

Changes to MATPOWER in each released version are summarized in the
[release notes](docs/relnotes), found in `docs/relnotes` and in
Appendix H of the [MATPOWER User's Manual][8]. A complete, detailed
change log, even for unreleased versions, is available in the
[`CHANGES.md`][11] file.


Contributing
------------

**For Contributing to MATPOWER-FUBM**
Please contact:
1. Abraham Alvarez-Bustos at abraham.alvarez-bustos@durham.ac.uk 
                          or snoop_and@hotmail.com
2. Behzad Kazemtabrizi at behzad.kazemtabrizi@durham.ac.uk
3. Mahmoud Shahbazi at mahmoud.shahbazi@durham.ac.uk


Sponsoring the FUBM Project
-------------------------------

If you have found the FUBM to be valuable, please consider supporting
the project by becoming a sponsor. Please contact:
1. Abraham Alvarez-Bustos at abraham.alvarez-bustos@durham.ac.uk 
                          or snoop_and@hotmail.com
2. Behzad Kazemtabrizi at behzad.kazemtabrizi@durham.ac.uk
3. Mahmoud Shahbazi at mahmoud.shahbazi@durham.ac.uk

Any contributions from the community or other sponsors free us to focus on
that support and the development of valuable new features.

-----------------------------------------
MATPOWER-FUBM Publications and Tech Notes
-----------------------------------------

10. A. Alvarez-Bustos and B. Kazemtabrizi, ["Flexible general branch 
    model unified power flow algorithm for future flexible AC/DC 
    networks,"][42] *2018 IEEE International Conference on Environment
    and Electrical Engineering (EEEIC / I&CPS Europe)*, Palermo, Italy,
    Jun. 2018.
    doi: [10.1109/EEEIC.2018.8493705][43].

11. *Under Review* A. Alvarez-Bustos, B. Kazemtabrizi, M. Shahbazi, and
    E. Acha-Daza, "Universal Branch Model for the Solution of Optimal 
    Power Flows in Hybrid AC/DC Grids," *International Journal of 
    Electrical Power and Energy Systems*, vol. XX, no. X, pp. XX–XX,
    Month. 20XX.  
    doi: XX.XXXX/XXXXX.20XX.XXXXXXX XX.


12. *To Be Published* A. Alvarez-Bustos, B. Kazemtabrizi, M. Shahbazi, and
    R. D. Zimmerman, "MATPOWER-FUBM: Flexible Universal Branch Model for 
    Matpower’s Optimal Power Flow and Power Flow Tools for Hybrid AC/DC 
    Power Systems Research and Education," *Power Systems, IEEE Transactions
    on*, vol. XX, no. X, pp. XX–XX, Month. 20XX.  
    doi: XX.XXXX/XXXXX.20XX.XXXXXXX XX.

13. A. Alvarez-Bustos, "AC/DC Optimal Power Flow and Power flow Equations and
    their Derivatives in Complex Matrix Notation using FUBM for MATPOWER," 
    *MATPOWER-FUBM Technical Note*, Ago 2020.  
    Available:
    https://matpower.org/docs/TN-FUBM-Derivatives.pdf  
    doi: XX.XXXX/XXXXX.20XX.XXXXXXX XX.

14. A. Alvarez-Bustos, "MATPOWER-FUBM Quick Guide," *MATPOWER-FUBM Quick Guide*,
    Ago 2020. Available:
    https://matpower.org/docs/MATPOWER-FUBM-Quick-Guide.pdf  
    doi: XX.XXXX/XXXXX.20XX.XXXXXXX XX.

[Citing MATPOWER-FUBM][31]
---------------------

We request that publications derived from the use of MATPOWER-FUBM,
the included data files, explicitly acknowledge that fact by citing the
appropriate paper(s) and the softwares itself. Please notice that every
publication derived from the use of MATPOWER-FUBM must cite, the FUBM, 
the MATPOWER-FUBM and also MATPOWER. 

#### Papers

All publications derived from the use of *MATPOWER-FUBM*, the FUBM model or
the included FUBM data files, should cite the following papers:

>   *To Be Published* A. Alvarez-Bustos, B. Kazemtabrizi, M. Shahbazi, and
    R. D. Zimmerman, "MATPOWER-FUBM: Flexible Universal Branch Model for 
    Matpower’s Optimal Power Flow and Power Flow Tools for Hybrid AC/DC 
    Power Systems Research and Education," *Power Systems, IEEE Transactions
    on*, vol. XX, no. X, pp. XX–XX, Month. 20XX.  
    doi: XX.XXXX/XXXXX.20XX.XXXXXXX XX.

>   *Under Review* A. Alvarez-Bustos, B. Kazemtabrizi, M. Shahbazi, and
    E. Acha-Daza, "Universal Branch Model for the Solution of Optimal 
    Power Flows in Hybrid AC/DC Grids," *International Journal of 
    Electrical Power and Energy Systems*, vol. XX, no. X, pp. XX–XX,
    Month. 20XX.  
    doi: XX.XXXX/XXXXX.20XX.XXXXXXX XX.

All publications derived from the use of *MATPOWER*, or the included data
files, should cite the 2011 MATPOWER paper:

>   R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, "MATPOWER:
    Steady-State Operations, Planning and Analysis Tools for Power Systems
    Research and Education," *Power Systems, IEEE Transactions on*, vol. 26,
    no. 1, pp. 12-19, Feb. 2011.  
    doi: [10.1109/TPWRS.2010.2051168][13]


Work making specific reference to the [MATPOWER Interior Point Solver
(MIPS)][32] should also cite:

>   H. Wang, C. E. Murillo-Sánchez, R. D. Zimmerman, R. J. Thomas, "On
    Computational Issues of Market-Based Optimal Power Flow," *Power Systems,
    IEEE Transactions on*, vol. 22, no. 3, pp. 1185-1193, Aug. 2007.  
    doi: [10.1109/TPWRS.2007.901301][17]

NOTE: Some of the case files included with MATPOWER request the citation
of additional publications. This includes the ACTIVSg, PEGASE, and RTE
cases. Details are available in the help text at the top of the
corresponding case files.

#### Software
To cite the MATPOWER-FUBM software, should cite the following and MATPOWER:

>   *To Be Published* A. Alvarez-Bustos, B. Kazemtabrizi, M. Shahbazi, and
    R. D. Zimmerman, "MATPOWER-FUBM: Flexible Universal Branch Model for 
    Matpower’s Optimal Power Flow and Power Flow Tools for Hybrid AC/DC 
    Power Systems Research and Education," *Power Systems, IEEE Transactions
    on*, vol. XX, no. X, pp. XX–XX, Month. 20XX.  
    doi: XX.XXXX/XXXXX.20XX.XXXXXXX XX.

>   A. Alvarez-Bustos (2020). MATPOWER-FUBM (Version 1.0)
    [Software]. Available: https://github.com/AbrahamAlvarezB/matpower-fubm
    doi: XX.XXXX/XXXXX.20XX.XXXXXXX 

To cite the MATPOWER software generally, without reference to a specific
version, use the following citation and DOI, with *\<YEAR\>* replaced by the
year of the most recent release:

>   R. D. Zimmerman, C. E. Murillo-Sanchez (*\<YEAR\>*). *MATPOWER*
    [Software]. Available: https://matpower.org  
    doi: [10.5281/zenodo.3236535][33]

A list of versions with release dates and version-specific DOI's can be
found via the general DOI at https://doi.org/10.5281/zenodo.3236535.

For the sake of reproducibility of research results, it is best to cite
the specific commit version of the software used and date.


#### User's Manuals

The MATPOWER, MIPS and MOST User's Manuals should also be cited
explicitly in work that refers to or is derived from their content. As
with the software, the citation and DOI can be version-specific or
general, as appropriate. For version 7.0 of the [MATPOWER User's Manual][8],
use:

>   R. D. Zimmerman, C. E. Murillo-Sanchez. *MATPOWER User's Manual,
    Version 7.0.* 2019.  
    [Online]. Available: https://matpower.org/docs/MATPOWER-manual-7.0.pdf  
    doi: [10.5281/zenodo.3251118](https://doi.org/10.5281/zenodo.3251118)

For a version non-specific citation, use the following citation and DOI,
with *\<YEAR\>* replaced by the year of the most recent release:

>   R. D. Zimmerman, C. E. Murillo-Sanchez. *MATPOWER User's Manual.* *\<YEAR\>*.  
    [Online]. Available: https://matpower.org/docs/MATPOWER-manual.pdf  
    doi: [10.5281/zenodo.3236519][34]

A list of versions of the User's Manual with release dates and
version-specific DOI's can be found via the general DOI at
https://doi.org/10.5281/zenodo.3236519.

For information on citing the MIPS or MOST User's Manuals, please see
the [`mips/CITATION`][35] and [`most/CITATION`][36] files, respectively.

#### Recommendation

In the interest of facilitating research reproducibility and thereby
increasing the value of your MATPOWER-related research publications, we
strongly encourage you to also publish, whenever possible, all of the
code and data required to generate the results you are publishing.
[Zenodo/GitHub][37] and [IEEE DataPort][38] are two of [many available
options][39].


E-mail Lists
------------

There are two e-mail lists available to serve the MATPOWER community:

- [**Discussion List**][26] ([MATPOWER-L][26]) – to facilitate discussion
  among MATPOWER users and provide a forum for help with MATPOWER
  related questions

- [**Developer List**][27] ([MATPOWER-DEV-L][27]) – to provide a forum
  for discussion among MATPOWER users and developers related to the
  development of the MATPOWER software or proposed contributions

For details see the [Mailing Lists section][28] of the
[MATPOWER website][1].

Please select the most appropriate list for your post and do *not*
cross-post to both Discussion and Developer lists. Bug reports,
software patches, proposed enhancements, etc. should be submitted to
the [issue tracker on GitHub][29].


Optional Packages
-----------------

There are numerous optional packages to enhance the performance of
MATPOWER that must be installed separately. The terms of use and
license agreements vary. Some are free of charge for all to use,
others are only free for academic use, and others may require a
commercial license. Please see Appendix G of the [MATPOWER User's
Manual][8] for details.


License and Terms of Use
------------------------

MATPOWER is distributed as open-source under the [3-clause BSD license][30].

---

[1]: https://matpower.org
[2]: https://github.com/MATPOWER/matpower
[3]: https://www.mathworks.com/
[4]: https://www.gnu.org/software/octave/
[5]: https://git-scm.com/downloads
[6]: https://git-scm.com
[7]: CONTRIBUTING.md
[7a]: https://hub.docker.com/
[7b]: https://www.docker.com
[7c]: https://hub.docker.com/r/matpower/matpower-desktop
[7d]: https://github.com/MATPOWER/matpower-extras
[7e]: docker/MATPOWER-Docker.md
[8]: docs/MATPOWER-manual.pdf
[9]: https://matpower.org/docs/ref/
[10]: most/docs/MOST-manual.pdf
[11]: CHANGES.md
[12]: https://matpower.org/docs/MATPOWER-paper.pdf
[13]: https://doi.org/10.1109/TPWRS.2010.2051168
[14]: https://matpower.org/docs/MATPOWER-OPF.pdf
[15]: https://doi.org/10.1109/PES.2009.5275967
[16]: https://matpower.org/docs/MATPOWER-OPF-slides.pdf
[17]: https://doi.org/10.1109/TPWRS.2007.901301
[18]: https://doi.org/10.1109/TSG.2013.2281001
[19]: https://doi.org/10.1109/TSTE.2018.2865454
[20]: https://matpower.org/docs/TN1-OPF-Auctions.pdf
[21]: https://matpower.org/docs/TN2-OPF-Derivatives.pdf
[22]: https://matpower.org/docs/TN3-More-OPF-Derivatives.pdf
[23]: https://matpower.org/docs/TN4-OPF-Derivatives-Cartesian.pdf
[24]: https://github.com/MATPOWER/most
[26]: https://matpower.org/mailing-lists/#discusslist
[27]: https://matpower.org/mailing-lists/#devlist
[28]: https://matpower.org/mailing-lists
[29]: https://github.com/MATPOWER/matpower/issues
[30]: LICENSE
[31]: CITATION
[32]: https://github.com/MATPOWER/mips
[33]: https://doi.org/10.5281/zenodo.3236535
[34]: https://doi.org/10.5281/zenodo.3236519
[35]: mips/CITATION
[36]: most/CITATION
[37]: https://guides.github.com/activities/citable-code/
[38]: https://ieee-dataport.org
[39]: https://www.re3data.org
[40]: https://github.com/AbrahamAlvarezB/matpower-fubm
[41]: docs/MATPOWER-FUBM-Quick-Guide.pdf
[42]: https://ieeexplore.ieee.org/document/8493705
[43]: https://doi.org/10.1109/EEEIC.2018.8493705

[logo]: docs/src/images/MATPOWER-md.png
[logo-fubm]: docs/src/images/MATPOWER-FUBM-md.png
