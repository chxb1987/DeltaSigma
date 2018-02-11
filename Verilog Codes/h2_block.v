`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    17:06:05 01/03/2018 
// Design Name: 
// Module Name:    h2_block 
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
module h2_block(
    input [15:0] in,
    input CLK,
	 input [15:0] reset,
    output [15:0] out
    );
wire CLK_;
wire [15:0] in_s,in_sd,in_sdnot,in_sdnotd,out_tmp;
not n2(CLK_,CLK);
shifter_1 s1(in,1'b1,in_s);
reg_16 r1(in_s,CLK,16'b0000000000000000,reset,in_sd);
two_cpl_16 n1(in_sd,in_sdnot);
reg_16 r2(in_sdnot,CLK_,16'b0000000000000000,reset,in_sdnotd);
add_16 a1(in_sdnot,in_s,out_tmp);
reg_16 r3(out_tmp,CLK,16'b0000000000000000,reset,out);

endmodule
