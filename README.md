# one-layer-empc
Implementation of one-layer economic model predictive control.

A common approach for model predictive control is to divide real-time optimization (RTO) and model predictive control (MPC) layers. Another approach is to combine both layers into an economic model predictive control (EMPC). This code is an implementation of EMPC without RTO layer.

There are two case examples, which are a cyclic process and distillation column A (plus CSTR). They are located inside the folder models. Main codes for ideal multistage are:
- NMPCCyclic.m
- MultistageNMPCDistColA.m 

For sensitivity-based codes (you need MINOS QP solver to run) are:
- AsMsNmpcCyclic.m
- AsMsNMPCDistColA.m.  
