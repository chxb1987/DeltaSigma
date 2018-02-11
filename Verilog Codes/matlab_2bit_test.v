`timescale 1ns / 1ps

////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:
//
// Create Date:   15:05:17 01/11/2018
// Design Name:   matlab_2bit
// Module Name:   D:/LEARN/BTP/matlab_2bit_test.v
// Project Name:  BTP
// Target Device:  
// Tool versions:  
// Description: 
//
// Verilog Test Fixture created by ISE for module: matlab_2bit
//
// Dependencies:
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
////////////////////////////////////////////////////////////////////////////////

module matlab_2bit_test;

	// Inputs
	reg [15:0] v_tmp12;
	reg [15:0] v_12;
	reg CLK;
	reg [15:0] reset;

	// Outputs
	wire [15:0] v_lsli;

	// Instantiate the Unit Under Test (UUT)
	matlab_2bit uut (
		.v_tmp12(v_tmp12), 
		.v_12(v_12), 
		.CLK(CLK), 
		.reset(reset), 
		.v_lsli(v_lsli)
	);

	initial begin
		// Initialize Inputs
		v_tmp12 = 16'b1000010101100000;
		v_12 = 16'b1000010101100000;
		CLK = 0;
		reset = 16'b1111111111111111;

		// Wait 100 ns for global reset to finish
		#10;
		reset = 16'b0000000000000000;
		#10;
      v_tmp12 = 16'b0110110101100000;
		v_12 = 16'b0110110101100000;
		#20;
      v_tmp12 = 16'b0010110110010000;
		v_12 = 16'b0010110110010000;
		#20;
      v_tmp12 = 16'b1111100100000000;
		v_12 = 16'b1111100100000000;
		#20;
      v_tmp12 = 16'b1111100111100000;
		v_12 = 16'b1111100111100000;
        
		// Add stimulus here

	end
always begin #10; CLK <= ~CLK; end      
endmodule

