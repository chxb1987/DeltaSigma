`timescale 1ns / 1ps

////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:
//
// Create Date:   16:10:16 01/10/2018
// Design Name:   h2_block
// Module Name:   D:/LEARN/BTP/h2_block_test.v
// Project Name:  BTP
// Target Device:  
// Tool versions:  
// Description: 
//
// Verilog Test Fixture created by ISE for module: h2_block
//
// Dependencies:
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
////////////////////////////////////////////////////////////////////////////////

module h2_block_test;

	// Inputs
	reg [15:0] in;
	reg CLK;
	reg [15:0] reset;

	// Outputs
	wire [15:0] out;

	// Instantiate the Unit Under Test (UUT)
	h2_block uut (
		.in(in), 
		.CLK(CLK), 
		.reset(reset), 
		.out(out)
	);

	initial begin
		// Initialize Inputs
		in = 16'b0000000000000000;
		CLK = 0;
		reset = 16'b1111111111111111;

		// Wait 100 ns for global reset to finish
		#10;
		reset = 16'b0000000000000000;
		#10;
      in = 16'b0111010000000000;
		#20;
      in = 16'b1110000000011000;
		#20;
      in = 16'b1110010110111000;
		#20;
      in = 16'b0000000001110000;
      #20;
      in = 16'b0000000000000000;  
		// Add stimulus here

	end
always begin #10; CLK <= ~CLK; end      
endmodule

