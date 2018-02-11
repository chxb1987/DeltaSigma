`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    14:18:36 01/07/2018 
// Design Name: 
// Module Name:    mux_16_struct 
// Project Name: 
// Target Devices: 
// Tool versions: 
// Description: 
//
// Dependencies: 
//
// Revision: 
// Revision 0.01 - File Created
// Additional Comments: 
//
//////////////////////////////////////////////////////////////////////////////////
module mux_16_struct(
    input [15:0] a,
    input [15:0] b,
    input en,
    output [15:0] c
    );
wire [15:0] d, e;
wire en_;
not(en_,en);
and a1(d[0],en_,a[0]);
and a2(d[1],en_,a[1]);
and a3(d[2],en_,a[2]);
and a4(d[3],en_,a[3]);
and a5(d[4],en_,a[4]);
and a6(d[5],en_,a[5]);
and a7(d[6],en_,a[6]);
and a8(d[7],en_,a[7]);
and a9(d[8],en_,a[8]);
and a10(d[9],en_,a[9]);
and a11(d[10],en_,a[10]);
and a12(d[11],en_,a[11]);
and a13(d[12],en_,a[12]);
and a14(d[13],en_,a[13]);
and a15(d[14],en_,a[14]);
and a16(d[15],en_,a[15]);

and b1(e[0],en,b[0]);
and b2(e[1],en,b[1]);
and b3(e[2],en,b[2]);
and b4(e[3],en_,b[3]);
and b5(e[4],en_,b[4]);
and b6(e[5],en_,b[5]);
and b7(e[6],en_,b[6]);
and b8(e[7],en_,b[7]);
and b9(e[8],en_,b[8]);
and b10(e[9],en_,b[9]);
and b11(e[10],en_,b[10]);
and b12(e[11],en_,b[11]);
and b13(e[12],en_,b[12]);
and b14(e[13],en_,b[13]);
and b15(e[14],en_,b[14]);
and b16(e[15],en_,b[15]);

or o1(c[0],d[0],e[0]);
or o2(c[1],d[1],e[1]);
or o3(c[2],d[2],e[2]);
or o4(c[3],d[3],e[3]);
or o5(c[4],d[4],e[4]);
or o6(c[5],d[5],e[5]);
or o7(c[6],d[6],e[6]);
or o8(c[7],d[7],e[7]);
or o9(c[8],d[8],e[8]);
or o10(c[9],d[9],e[9]);
or o11(c[10],d[10],e[10]);
or o12(c[11],d[11],e[11]);
or o13(c[12],d[12],e[12]);
or o14(c[13],d[13],e[13]);
or o15(c[14],d[14],e[14]);
or o16(c[15],d[15],e[15]);

endmodule
