OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];
creg c[4];

swap q[0],q[1];
s q[3];
h q[0];
swap q[1],q[2];
z q[0];
cswap q[1],q[2],q[3];
cx q[0],q[1];
s q[2];
z q[0];
ccx q[1],q[2],q[3];
h q[1];
cx q[0],q[1];
s q[2];
z q[3];
h q[0];
