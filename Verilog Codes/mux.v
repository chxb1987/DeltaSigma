`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    17:08:09 01/04/2018 
// Design Name: 
// Module Name:    mux 
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
module mux(
    input a,
    input b,
    input en,
    output out
    );
wire d,e,en_;
not n1(en_,en);
and a1(d,en_,a);
and b1(e,en,b);
or o1(out,d,e);

endmodule
