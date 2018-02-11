`timescale 1ns / 1ps

////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:
//
// Create Date:   15:13:06 01/17/2018
// Design Name:   FA
// Module Name:   D:/LEARN/BTP/FA_test.v
// Project Name:  BTP
// Target Device:  
// Tool versions:  
// Description: 
//
// Verilog Test Fixture created by ISE for module: FA
//
// Dependencies:
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
////////////////////////////////////////////////////////////////////////////////

module FA_test;

	// Inputs
	reg A;
	reg B;
	reg C0;

	// Outputs
	wire S;
	wire C1;

	// Instantiate the Unit Under Test (UUT)
	FA uut (
		.A(A), 
		.B(B), 
		.C0(C0), 
		.S(S), 
		.C1(C1)
	);

	initial begin
		// Initialize Inputs
		A = 0;
		B = 0;
		C0 = 0;

		// Wait 100 ns for global reset to finish
		#20;
      C0 = 1;
 		#20;
      B = 1;
		C0 = 0;
		#20;
      C0 = 1;
		#20;
      A = 1;
		B = 0;
		C0 = 0;
		#20;
      C0 = 1;
		#20;
		B = 1;
      C0 = 0;
		#20;
      C0 = 1;
		// Add stimulus here

	end
      
endmodule

