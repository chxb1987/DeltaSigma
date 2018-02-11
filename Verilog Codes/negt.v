`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    17:02:40 01/03/2018 
// Design Name: 
// Module Name:    negt 
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
module negt(
    input [15:0] A,
    output [15:0] A_
    );
xor x1(A_[0],A[0],1'b1);
xor x2(A_[1],A[1],1'b1);
xor x3(A_[2],A[2],1'b1);
xor x4(A_[3],A[3],1'b1);
xor x5(A_[4],A[4],1'b1);
xor x6(A_[5],A[5],1'b1);
xor x7(A_[6],A[6],1'b1);
xor x8(A_[7],A[7],1'b1);
xor x9(A_[8],A[8],1'b1);
xor x10(A_[9],A[9],1'b1);
xor x11(A_[10],A[10],1'b1);
xor x12(A_[11],A[11],1'b1);
xor x13(A_[12],A[12],1'b1);
xor x14(A_[13],A[13],1'b1);
xor x15(A_[14],A[14],1'b1);
xor x16(A_[15],A[15],1'b1);

endmodule
