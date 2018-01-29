`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    15:23:46 01/09/2018 
// Design Name: 
// Module Name:    reg_16 
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
module reg_16(
    input [15:0] a,
    input CLK,
    input [15:0] set,
    input [15:0] reset,
    output [15:0] b
    );
wire [15:0] b_;
dff_struct d1(a[0],CLK,set[0],reset[0],b[0],b_[0]);
dff_struct d2(a[1],CLK,set[1],reset[1],b[1],b_[1]);
dff_struct d3(a[2],CLK,set[2],reset[2],b[2],b_[2]);
dff_struct d4(a[3],CLK,set[3],reset[3],b[3],b_[3]);
dff_struct d5(a[4],CLK,set[4],reset[4],b[4],b_[4]);
dff_struct d6(a[5],CLK,set[5],reset[5],b[5],b_[5]);
dff_struct d7(a[6],CLK,set[6],reset[6],b[6],b_[6]);
dff_struct d8(a[7],CLK,set[7],reset[7],b[7],b_[7]);
dff_struct d9(a[8],CLK,set[8],reset[8],b[8],b_[8]);
dff_struct d10(a[9],CLK,set[9],reset[9],b[9],b_[9]);
dff_struct d11(a[10],CLK,set[10],reset[10],b[10],b_[10]);
dff_struct d12(a[11],CLK,set[11],reset[11],b[11],b_[11]);
dff_struct d13(a[12],CLK,set[12],reset[12],b[12],b_[12]);
dff_struct d14(a[13],CLK,set[13],reset[13],b[13],b_[13]);
dff_struct d15(a[14],CLK,set[14],reset[14],b[14],b_[14]);
dff_struct d16(a[15],CLK,set[15],reset[15],b[15],b_[15]);

endmodule
