`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    16:48:41 01/03/2018 
// Design Name: 
// Module Name:    shifter_1 
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
module shifter_1(
    input [15:0] xin,
    input en,
    output [15:0] xout
    );
mux m1(xin[0],xin[1],en,xout[0]);
mux m2(xin[1],xin[2],en,xout[1]);
mux m3(xin[2],xin[3],en,xout[2]);
mux m4(xin[3],xin[4],en,xout[3]);
mux m5(xin[4],xin[5],en,xout[4]);
mux m6(xin[5],xin[6],en,xout[5]);
mux m7(xin[6],xin[7],en,xout[6]);
mux m8(xin[7],xin[8],en,xout[7]);
mux m9(xin[8],xin[9],en,xout[8]);
mux m10(xin[9],xin[10],en,xout[9]);
mux m11(xin[10],xin[11],en,xout[10]);
mux m12(xin[11],xin[12],en,xout[11]);
mux m13(xin[12],xin[13],en,xout[12]);
mux m14(xin[13],xin[14],en,xout[13]);
mux m15(xin[14],xin[15],en,xout[14]);
mux m16(xin[15],xin[15],en,xout[15]);

endmodule
