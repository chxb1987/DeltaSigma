`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    14:26:11 01/05/2018 
// Design Name: 
// Module Name:    h1_block_matlab_imperfect 
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
module h1_block_matlab_imperfect(
    input [15:0] in,
    input CLK,
	 input [15:0] reset,
    output [15:0] out
    );
wire [15:0] in_d,in_ds,in_dss,in_dsnet,in_dd,in_dds,in_ddsnot;
reg_16 r1(in,CLK,16'b0000000000000000,reset,in_d);
reg_16 r2(in_d,CLK,16'b0000000000000000,reset,in_dd);
shifter_1 s1(in_d,1'b1,in_ds);
shifter_1 s2(in_ds,1'b1,in_dss);
shifter_1 s3(in_dd,1'b1,in_dds);
two_cpl_16 n1(in_dds,in_ddsnot);
add_16 a1(in_ds,in_dss,in_dsnet);
add_16 a2(in_dsnet,in_ddsnot,out);

endmodule
