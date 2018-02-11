`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    16:44:09 01/03/2018 
// Design Name: 
// Module Name:    dff_struct 
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
module dff_struct(
    input din,
    input CLK,
    input set,
    input reset,
    output q,
    output q_
    );
wire a1,a2,a3,a4,set_,reset_; 
not m1(set_,set);
not m2(reset_,reset);
nand_3 n1(.a(set_),.b(a1),.c(a2),.out(a3));
nand_3 n2(.a(reset_),.b(a3),.c(CLK),.out(a2));
nand_3 n3(.a(a2),.b(CLK),.c(a1),.out(a4));
nand_3 n4(.a(a4),.b(reset_),.c(din),.out(a1));
nand_3 n5(.a(a2),.b(set_),.c(q_),.out(q));
nand_3 n6(.a(a4),.b(reset_),.c(q),.out(q_));

endmodule
