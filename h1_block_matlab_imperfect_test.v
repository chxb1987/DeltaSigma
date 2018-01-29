`timescale 1ns / 1ps

////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:
//
// Create Date:   15:14:42 01/09/2018
// Design Name:   h1_block_matlab_imperfect
// Module Name:   D:/LEARN/BTP/h1_block_matlab_imperfect_test.v
// Project Name:  BTP
// Target Device:  
// Tool versions:  
// Description: 
//
// Verilog Test Fixture created by ISE for module: h1_block_matlab_imperfect
//
// Dependencies:
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
////////////////////////////////////////////////////////////////////////////////

module h1_block_matlab_imperfect_test;

	// Inputs
	reg [15:0] in;
	reg CLK;
	reg [15:0] reset;

	// Outputs
	wire [15:0] out;

	// Instantiate the Unit Under Test (UUT)
	h1_block_matlab_imperfect uut (
		.in(in), 
		.CLK(CLK), 
		.reset(reset),
		.out(out)
	);

	initial begin
		// Initialize Inputs
		in = 16'b1000010101100000;
		CLK = 0;
		reset = 16'b1111111111111111;

		// Wait 100 ns for global reset to finish
		#10;
		reset = 16'b0000000000000000;
		#10;
      in = 16'b0110110101100000;
		#20;
      in = 16'b0010110110010000;
		#20;
      in = 16'b1111100100000000;
		#20;
      in = 16'b1111100111100000;
		// Add stimulus here

	end
always begin #10; CLK <= ~CLK; end      
endmodule

