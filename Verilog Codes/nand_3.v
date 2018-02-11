`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    17:20:07 01/03/2018 
// Design Name: 
// Module Name:    nand_3 
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
module nand_3(
    input a,
    input b,
    input c,
    output out
    );
wire d;
nand n1(d,a,b);
or(out,~c,d);

endmodule
