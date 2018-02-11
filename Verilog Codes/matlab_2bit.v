`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    14:53:16 01/06/2018 
// Design Name: 
// Module Name:    matlab_2bit 
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
/*module matlab_2bit(
    input v_tmp,
	 input [15:0] v,
    input CLK,
	 input [15:0] reset,
    output [15:0] v_lsli
    );
wire [15:0] v_tmp12,v_12,v_tmp_1,v_1,v_2;
wire v0_;

mux_16_struct m1(16'b1100000000000000,16'b0100000000000000,v_tmp,v_tmp12[15:13]);
not n1(v_12[15],v[1]);
not n2(v0_,v[0]);
not n3(v_12[14],v0_);
//assign v_12[14] = v[0];
not n4(v_12[13],1'b0);
//assign v_12[13] = 1'b1;
not n5(v_12[12],1'b1);
not n6(v_12[11],1'b1);
not n7(v_12[10],1'b1);
not n8(v_12[9],1'b1);
not n9(v_12[8],1'b1);
not n10(v_12[7],1'b1);
not n11(v_12[6],1'b1);
not n12(v_12[5],1'b1);
not n13(v_12[4],1'b1);
not n14(v_12[3],1'b1);
not n15(v_12[2],1'b1);
not n16(v_12[1],1'b1);
not n17(v_12[0],1'b1);
//assign v_12[12:0] = 13'b0000000000000;

h1_block_matlab_imperfect h1(v_tmp12,CLK,reset,v_tmp_1);
h2_block h21(v_12,CLK,reset,v_1);
h2_block h22(v_1,CLK,reset,v_2);

add_16 a1(v_tmp_1,v_2,v_lsli);

endmodule
*/
module matlab_2bit(
    input [15:0] v_tmp12,
	 input [15:0] v_12,
    input CLK,
	 input [15:0] reset,
    output [15:0] v_lsli
    );
wire [15:0] v_tmp_1,v_tmp_1f1,v_tmp_1f2,v_tmp_1f3,v_tmp_1f4,v_tmp_1fnet;
wire [15:0] v_1,v_1f1,v_1f2,v_2;
wire CLK_;
not n1(CLK_,CLK);

h1_block_matlab_imperfect h1(v_tmp12,CLK,reset,v_tmp_1);
h2_block h21(v_12,CLK,reset,v_1);

reg_16 r1(v_1,CLK,16'b0000000000000000,reset,v_1f1);
//reg_16 r2(v_1f1,CLK_,16'b0000000000000000,reset,v_1f2);

h2_block h22(v_1f1,CLK_,reset,v_2);

reg_16 r3(v_tmp_1,CLK_,16'b0000000000000000,reset,v_tmp_1f1);
reg_16 r4(v_tmp_1f1,CLK,16'b0000000000000000,reset,v_tmp_1f2);
reg_16 r5(v_tmp_1f2,CLK_,16'b0000000000000000,reset,v_tmp_1f3);
reg_16 r6(v_tmp_1f3,CLK,16'b0000000000000000,reset,v_tmp_1f4);
reg_16 r7(v_tmp_1f4,CLK_,16'b0000000000000000,reset,v_tmp_1fnet);

add_16 a1(v_tmp_1fnet,v_2,v_lsli);

endmodule