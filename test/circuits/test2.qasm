OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];
creg c[5];

u1(pi/2) q[0];
h q[1];
s q[2];
swap q[2],q[3];
sdg q[3];
cx q[3],q[4];
