`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    15:49:10 01/09/2018 
// Design Name: 
// Module Name:    two_cpl_16 
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
module two_cpl_16(
    input [15:0] A,
    output [15:0] B
    );
wire [15:0] A_;
negt n1(A,A_);
add_16 a100(A_,16'b0000000000000001,B);

endmodule
