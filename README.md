# Power_Flow_in_Hybrid_ACDC_networks
The source code of the paper titled: 'General and unified model of the power flow problem in hybrid AC/DC Networks' (https://ieeexplore.ieee.org/document/10475554)

The power flow model includes the AC network, DC network and the interfacing converters in a unified way. The interfacing converters are modelled generically and can operate in different control modes: voltage mode or power mode.

The following examples are given:
- Run the file **Main_Microgrid.m** for the power flow of the 26-node hybrid AC/DC microgrid
- Run the file **Main_AsyncHVDC.m** for the power flow of the HVDC system that interfaces two asynchronous HVAC networks

A comparison with a Matpower-based model is given in the folder named 'Comparison'. The computational time of our method is around 10x faster with the same accuracy.
The additional results presented in the paper appendix are given in the folder named 'Appendix'

In further research, the power flow model is extended with a more detailed losses model and the option to include grid-forming interfacing converters.
You can contact the authors for more information. willem.lambrichts@gmail.com
