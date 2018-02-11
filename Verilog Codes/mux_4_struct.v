`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    14:51:01 01/06/2018 
// Design Name: 
// Module Name:    mux_4_struct 
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
module mux_4_struct(
    input [3:0] a,
    input [3:0] b,
    input en,
    output [3:0] c
    );
wire [3:0] d, e;
wire en_;
not(en_,en);
and a1(d[0],en_,a[0]);
and a2(d[1],en_,a[1]);
and a3(d[2],en_,a[2]);
and a4(d[3],en_,a[3]);

and b1(e[0],en,b[0]);
and b2(e[1],en,b[1]);
and b3(e[2],en,b[2]);
and b4(e[3],en,b[3]);

or o1(c[0],d[0],e[0]);
or o2(c[1],d[1],e[1]);
or o3(c[2],d[2],e[2]);
or o4(c[3],d[3],e[3]);

endmodule
